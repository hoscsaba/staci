#include "epanet_water_age.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace {
constexpr double kFlowEpsilon = 1.0e-12;
constexpr double kVolumeEpsilon = 1.0e-12;
}

EpanetWaterAgeModel::EpanetWaterAgeModel(
    std::vector<EpanetWaterAgeNode> nodes,
    std::vector<EpanetWaterAgeLink> links,
    double quality_timestep_s)
    : nodes_(std::move(nodes)), links_(std::move(links)),
      segments_(links_.size()), node_age_s_(nodes_.size(), 0.0),
      quality_timestep_s_(quality_timestep_s) {
    if (!(quality_timestep_s_ > 0.0) || !std::isfinite(quality_timestep_s_))
        throw std::invalid_argument("Water-age quality timestep must be positive and finite.");
    for (std::size_t index = 0; index < links_.size(); ++index) {
        const EpanetWaterAgeLink &link = links_[index];
        if (link.from_node >= nodes_.size() || link.to_node >= nodes_.size())
            throw std::invalid_argument("Water-age link has an invalid node index.");
        if (link.volume_m3 > kVolumeEpsilon)
            segments_[index].push_back(Segment{link.volume_m3, 0.0});
    }
}

void EpanetWaterAgeModel::advance(
    double duration_s, const std::vector<double> &link_flows_m3s) {
    if (duration_s < 0.0 || !std::isfinite(duration_s))
        throw std::invalid_argument("Water-age transport duration must be finite and non-negative.");
    if (link_flows_m3s.size() != links_.size())
        throw std::invalid_argument("Water-age flow vector does not match the link index.");
    double remaining = duration_s;
    while (remaining > 1.0e-9) {
        // A modest fixed cap limits numerical dispersion without tying the
        // global cost to the smallest/fastest pipe in a large network.
        double dt = std::min(remaining, std::min(quality_timestep_s_, 60.0));
        if (!(dt > 0.0) || !std::isfinite(dt))
            throw std::runtime_error("Water-age transport selected an invalid timestep.");
        step(dt, link_flows_m3s);
        remaining -= dt;
    }
    mix_nodes(0.0, link_flows_m3s);
}

void EpanetWaterAgeModel::step(
    double dt_s, const std::vector<double> &link_flows_m3s) {
    for (std::deque<Segment> &pipe : segments_)
        for (Segment &segment : pipe)
            segment.age_s += dt_s;
    mix_nodes(dt_s, link_flows_m3s);
    for (std::size_t index = 0; index < links_.size(); ++index) {
        const double flow = link_flows_m3s[index];
        if (std::abs(flow) <= kFlowEpsilon || links_[index].volume_m3 <= kVolumeEpsilon)
            continue;
        const std::size_t upstream = flow > 0.0
            ? links_[index].from_node : links_[index].to_node;
        displace(index, flow, dt_s, node_age_s_[upstream]);
    }
}

void EpanetWaterAgeModel::mix_nodes(
    double elapsed_s, const std::vector<double> &link_flows_m3s) {
    std::vector<double> fallback = node_age_s_;
    for (double &age : fallback)
        age += elapsed_s;
    for (std::size_t index = 0; index < nodes_.size(); ++index)
        if (nodes_[index].fixed_zero_age)
            fallback[index] = 0.0;

    std::vector<double> base_mass(nodes_.size(), 0.0);
    std::vector<double> base_flow(nodes_.size(), 0.0);
    for (std::size_t index = 0; index < links_.size(); ++index) {
        const double flow = link_flows_m3s[index];
        if (std::abs(flow) <= kFlowEpsilon || links_[index].volume_m3 <= kVolumeEpsilon)
            continue;
        const std::size_t downstream = flow > 0.0
            ? links_[index].to_node : links_[index].from_node;
        const double magnitude = std::abs(flow);
        base_mass[downstream] += magnitude * outlet_age(index, flow);
        base_flow[downstream] += magnitude;
    }

    std::vector<double> ages = fallback;
    const std::size_t iteration_limit = std::max<std::size_t>(4, nodes_.size() + 1);
    for (std::size_t iteration = 0; iteration < iteration_limit; ++iteration) {
        std::vector<double> mass = base_mass;
        std::vector<double> inflow = base_flow;
        for (std::size_t index = 0; index < links_.size(); ++index) {
            const double flow = link_flows_m3s[index];
            if (std::abs(flow) <= kFlowEpsilon || links_[index].volume_m3 > kVolumeEpsilon)
                continue;
            const std::size_t upstream = flow > 0.0
                ? links_[index].from_node : links_[index].to_node;
            const std::size_t downstream = flow > 0.0
                ? links_[index].to_node : links_[index].from_node;
            const double magnitude = std::abs(flow);
            mass[downstream] += magnitude * ages[upstream];
            inflow[downstream] += magnitude;
        }
        std::vector<double> next = fallback;
        double maximum_change = 0.0;
        for (std::size_t index = 0; index < nodes_.size(); ++index) {
            if (nodes_[index].fixed_zero_age)
                next[index] = 0.0;
            else if (inflow[index] > kFlowEpsilon)
                next[index] = mass[index] / inflow[index];
            maximum_change = std::max(maximum_change,
                                      std::abs(next[index] - ages[index]));
        }
        ages.swap(next);
        if (maximum_change < 1.0e-9)
            break;
    }
    node_age_s_.swap(ages);
}

double EpanetWaterAgeModel::outlet_age(
    std::size_t link_index, double flow_m3s) const {
    const std::deque<Segment> &pipe = segments_[link_index];
    if (pipe.empty()) {
        const std::size_t upstream = flow_m3s >= 0.0
            ? links_[link_index].from_node : links_[link_index].to_node;
        return node_age_s_[upstream];
    }
    return flow_m3s >= 0.0 ? pipe.back().age_s : pipe.front().age_s;
}

void EpanetWaterAgeModel::displace(
    std::size_t link_index, double flow_m3s, double dt_s,
    double inlet_age_s) {
    std::deque<Segment> &pipe = segments_[link_index];
    const double volume = std::min(std::abs(flow_m3s) * dt_s,
                                   links_[link_index].volume_m3);
    double remaining = volume;
    while (remaining > kVolumeEpsilon && !pipe.empty()) {
        Segment &outlet = flow_m3s > 0.0 ? pipe.back() : pipe.front();
        const double removed = std::min(remaining, outlet.volume_m3);
        outlet.volume_m3 -= removed;
        remaining -= removed;
        if (outlet.volume_m3 <= kVolumeEpsilon) {
            if (flow_m3s > 0.0) pipe.pop_back(); else pipe.pop_front();
        }
    }
    if (volume <= kVolumeEpsilon)
        return;
    const Segment inlet{volume, inlet_age_s};
    if (flow_m3s > 0.0) pipe.push_front(inlet); else pipe.push_back(inlet);
}

std::vector<double> EpanetWaterAgeModel::link_average_age_s(
    const std::vector<double> &link_flows_m3s) const {
    if (link_flows_m3s.size() != links_.size())
        throw std::invalid_argument("Water-age flow vector does not match the link index.");
    std::vector<double> result;
    result.reserve(links_.size());
    for (std::size_t index = 0; index < links_.size(); ++index) {
        if (segments_[index].empty()) {
            const std::size_t upstream = link_flows_m3s[index] >= 0.0
                ? links_[index].from_node : links_[index].to_node;
            result.push_back(node_age_s_[upstream]);
            continue;
        }
        double weighted_age = 0.0;
        double volume = 0.0;
        for (const Segment &segment : segments_[index]) {
            weighted_age += segment.volume_m3 * segment.age_s;
            volume += segment.volume_m3;
        }
        result.push_back(volume > kVolumeEpsilon ? weighted_age / volume : 0.0);
    }
    return result;
}

EpanetChemicalModel::EpanetChemicalModel(
    std::vector<EpanetChemicalNode> nodes,
    std::vector<EpanetChemicalLink> links,
    double quality_timestep_s)
    : nodes_(std::move(nodes)), links_(std::move(links)),
      segments_(links_.size()), node_concentration_kgm3_(nodes_.size(), 0.0),
      quality_timestep_s_(quality_timestep_s) {
    if (!(quality_timestep_s_ > 0.0) || !std::isfinite(quality_timestep_s_))
        throw std::invalid_argument("Chemical quality timestep must be positive and finite.");
    for (std::size_t index = 0; index < nodes_.size(); ++index) {
        const double value = nodes_[index].initial_concentration_kgm3;
        if (value < 0.0 || !std::isfinite(value))
            throw std::invalid_argument("Initial chemical concentration must be finite and non-negative.");
        node_concentration_kgm3_[index] = value;
    }
    for (std::size_t index = 0; index < links_.size(); ++index) {
        const EpanetChemicalLink &link = links_[index];
        if (link.from_node >= nodes_.size() || link.to_node >= nodes_.size())
            throw std::invalid_argument("Chemical link has an invalid node index.");
        if (link.volume_m3 < 0.0 || !std::isfinite(link.volume_m3) ||
            link.initial_concentration_kgm3 < 0.0 ||
            !std::isfinite(link.initial_concentration_kgm3) ||
            !std::isfinite(link.reaction_coefficient_per_s))
            throw std::invalid_argument("Chemical link properties must be finite and non-negative where applicable.");
        if (link.volume_m3 > kVolumeEpsilon)
            segments_[index].push_back(Segment{
                link.volume_m3, link.initial_concentration_kgm3});
    }
}

void EpanetChemicalModel::advance(
    double duration_s, const std::vector<double> &link_flows_m3s,
    const std::vector<double> &external_inflows_m3s,
    const std::vector<EpanetChemicalSource> &sources) {
    if (duration_s < 0.0 || !std::isfinite(duration_s))
        throw std::invalid_argument("Chemical transport duration must be finite and non-negative.");
    if (link_flows_m3s.size() != links_.size() ||
        external_inflows_m3s.size() != nodes_.size() || sources.size() != nodes_.size())
        throw std::invalid_argument("Chemical transport vectors do not match the network index.");
    double remaining = duration_s;
    while (remaining > 1.0e-9) {
        double dt = std::min(remaining, quality_timestep_s_);
        if (!(dt > 0.0) || !std::isfinite(dt))
            throw std::runtime_error("Chemical transport selected an invalid timestep.");
        step(dt, link_flows_m3s, external_inflows_m3s, sources);
        remaining -= dt;
    }
}

void EpanetChemicalModel::step(
    double dt_s, const std::vector<double> &link_flows_m3s,
    const std::vector<double> &external_inflows_m3s,
    const std::vector<EpanetChemicalSource> &sources) {
    for (std::size_t index = 0; index < segments_.size(); ++index) {
        const double multiplier = std::exp(
            links_[index].reaction_coefficient_per_s * dt_s);
        for (Segment &segment : segments_[index])
            segment.concentration_kgm3 = std::max(
                0.0, segment.concentration_kgm3 * multiplier);
    }
    // EPANET mixes the complete volume arriving during a quality timestep,
    // not merely the instantaneous concentration at a pipe outlet. This is
    // essential when a concentration front reaches a node part-way through a
    // timestep.
    std::vector<std::vector<std::size_t> > incoming(nodes_.size());
    std::vector<std::vector<std::size_t> > outgoing(nodes_.size());
    std::vector<std::size_t> indegree(nodes_.size(), 0);
    for (std::size_t index = 0; index < links_.size(); ++index) {
        const double flow = link_flows_m3s[index];
        if (std::abs(flow) <= kFlowEpsilon)
            continue;
        const std::size_t upstream = flow > 0.0
            ? links_[index].from_node : links_[index].to_node;
        const std::size_t downstream = flow > 0.0
            ? links_[index].to_node : links_[index].from_node;
        outgoing[upstream].push_back(index);
        incoming[downstream].push_back(index);
        ++indegree[downstream];
    }
    std::deque<std::size_t> ready;
    for (std::size_t node = 0; node < nodes_.size(); ++node)
        if (indegree[node] == 0)
            ready.push_back(node);
    std::vector<std::size_t> order;
    std::vector<bool> listed(nodes_.size(), false);
    while (!ready.empty()) {
        const std::size_t node = ready.front();
        ready.pop_front();
        if (listed[node]) continue;
        listed[node] = true;
        order.push_back(node);
        for (std::size_t link_index : outgoing[node]) {
            const double flow = link_flows_m3s[link_index];
            const std::size_t downstream = flow > 0.0
                ? links_[link_index].to_node : links_[link_index].from_node;
            if (indegree[downstream] > 0 && --indegree[downstream] == 0)
                ready.push_back(downstream);
        }
    }
    for (std::size_t node = 0; node < nodes_.size(); ++node)
        if (!listed[node]) order.push_back(node);

    const auto take_outflow = [&](std::size_t link_index, double flow,
                                  double requested, double &volume,
                                  double &mass) {
        std::deque<Segment> &pipe = segments_[link_index];
        double remaining = requested;
        while (remaining > kVolumeEpsilon && !pipe.empty()) {
            Segment &segment = flow > 0.0 ? pipe.back() : pipe.front();
            const double amount = std::min(remaining, segment.volume_m3);
            volume += amount;
            mass += amount * segment.concentration_kgm3;
            segment.volume_m3 -= amount;
            remaining -= amount;
            if (segment.volume_m3 <= kVolumeEpsilon) {
                if (flow > 0.0) pipe.pop_back(); else pipe.pop_front();
            }
        }
    };
    const auto add_inflow = [&](std::size_t link_index, double flow,
                                double volume, double concentration) {
        if (links_[link_index].volume_m3 <= kVolumeEpsilon)
            return;
        std::deque<Segment> &pipe = segments_[link_index];
        Segment *last = nullptr;
        if (!pipe.empty()) last = flow > 0.0 ? &pipe.front() : &pipe.back();
        if (last != nullptr && std::abs(last->concentration_kgm3 - concentration) < 1.0e-12)
            last->volume_m3 += volume;
        else if (flow > 0.0)
            pipe.push_front(Segment{volume, concentration});
        else
            pipe.push_back(Segment{volume, concentration});
    };

    for (std::size_t node : order) {
        double inflow_volume = 0.0;
        double inflow_mass = 0.0;
        for (std::size_t link_index : incoming[node]) {
            const double flow = link_flows_m3s[link_index];
            const double volume = std::abs(flow) * dt_s;
            if (links_[link_index].volume_m3 <= kVolumeEpsilon) {
                const std::size_t upstream = flow > 0.0
                    ? links_[link_index].from_node : links_[link_index].to_node;
                inflow_volume += volume;
                inflow_mass += volume * node_concentration_kgm3_[upstream];
            } else {
                take_outflow(link_index, flow, volume, inflow_volume, inflow_mass);
            }
        }
        const double external_volume = std::max(0.0, external_inflows_m3s[node]) * dt_s;
        inflow_volume += external_volume;
        if (!nodes_[node].fixed_external_concentration && inflow_volume > kVolumeEpsilon)
            node_concentration_kgm3_[node] = inflow_mass / inflow_volume;
        if (nodes_[node].fixed_external_concentration)
            node_concentration_kgm3_[node] = nodes_[node].initial_concentration_kgm3;

        double outflow_volume = 0.0;
        for (std::size_t link_index : outgoing[node])
            outflow_volume += std::abs(link_flows_m3s[link_index]) * dt_s;
        const EpanetChemicalSource &source = sources[node];
        if (outflow_volume > kVolumeEpsilon) {
            if (nodes_[node].fixed_external_concentration &&
                source.type == EpanetChemicalSourceType::Concentration)
                node_concentration_kgm3_[node] = source.strength_si;
            else if (nodes_[node].fixed_external_concentration &&
                     source.type == EpanetChemicalSourceType::Mass)
                node_concentration_kgm3_[node] =
                    source.strength_si * dt_s / outflow_volume;
            else if (nodes_[node].fixed_external_concentration &&
                     source.type == EpanetChemicalSourceType::Setpoint)
                node_concentration_kgm3_[node] = std::max(
                    source.strength_si - node_concentration_kgm3_[node], 0.0);
            else if (nodes_[node].fixed_external_concentration &&
                     source.type == EpanetChemicalSourceType::FlowPaced)
                node_concentration_kgm3_[node] = source.strength_si;
            else if (source.type == EpanetChemicalSourceType::Concentration &&
                external_volume > kVolumeEpsilon)
                node_concentration_kgm3_[node] +=
                    source.strength_si * external_volume / outflow_volume;
            else if (source.type == EpanetChemicalSourceType::Mass)
                node_concentration_kgm3_[node] += source.strength_si * dt_s / outflow_volume;
            else if (source.type == EpanetChemicalSourceType::Setpoint)
                node_concentration_kgm3_[node] = std::max(
                    node_concentration_kgm3_[node], source.strength_si);
            else if (source.type == EpanetChemicalSourceType::FlowPaced)
                node_concentration_kgm3_[node] += source.strength_si;
        }
        node_concentration_kgm3_[node] = std::max(0.0, node_concentration_kgm3_[node]);
        for (std::size_t link_index : outgoing[node]) {
            const double flow = link_flows_m3s[link_index];
            add_inflow(link_index, flow, std::abs(flow) * dt_s,
                       node_concentration_kgm3_[node]);
        }
    }
}

void EpanetChemicalModel::mix_nodes(
    const std::vector<double> &link_flows_m3s,
    const std::vector<double> &external_inflows_m3s,
    const std::vector<EpanetChemicalSource> &sources) {
    std::vector<double> base_mass(nodes_.size(), 0.0);
    std::vector<double> base_flow(nodes_.size(), 0.0);
    for (std::size_t index = 0; index < links_.size(); ++index) {
        const double flow = link_flows_m3s[index];
        if (std::abs(flow) <= kFlowEpsilon || links_[index].volume_m3 <= kVolumeEpsilon)
            continue;
        const std::size_t downstream = flow > 0.0
            ? links_[index].to_node : links_[index].from_node;
        const double magnitude = std::abs(flow);
        base_mass[downstream] += magnitude * outlet_concentration(index, flow);
        base_flow[downstream] += magnitude;
    }

    std::vector<double> values = node_concentration_kgm3_;
    const std::size_t iteration_limit = std::max<std::size_t>(4, nodes_.size() + 1);
    for (std::size_t iteration = 0; iteration < iteration_limit; ++iteration) {
        std::vector<double> mass = base_mass;
        std::vector<double> inflow = base_flow;
        for (std::size_t index = 0; index < links_.size(); ++index) {
            const double flow = link_flows_m3s[index];
            if (std::abs(flow) <= kFlowEpsilon || links_[index].volume_m3 > kVolumeEpsilon)
                continue;
            const std::size_t upstream = flow > 0.0
                ? links_[index].from_node : links_[index].to_node;
            const std::size_t downstream = flow > 0.0
                ? links_[index].to_node : links_[index].from_node;
            const double magnitude = std::abs(flow);
            mass[downstream] += magnitude * values[upstream];
            inflow[downstream] += magnitude;
        }
        std::vector<double> next = values;
        double maximum_change = 0.0;
        for (std::size_t index = 0; index < nodes_.size(); ++index) {
            const EpanetChemicalSource &source = sources[index];
            if (nodes_[index].fixed_external_concentration) {
                next[index] = nodes_[index].initial_concentration_kgm3;
            } else {
                const double external = std::max(0.0, external_inflows_m3s[index]);
                if (source.type == EpanetChemicalSourceType::Concentration && external > 0.0) {
                    mass[index] += external * source.strength_si;
                    inflow[index] += external;
                } else if (source.type == EpanetChemicalSourceType::Mass) {
                    mass[index] += source.strength_si;
                }
                if (inflow[index] > kFlowEpsilon)
                    next[index] = mass[index] / inflow[index];
                if (source.type == EpanetChemicalSourceType::Setpoint)
                    next[index] = source.strength_si;
                else if (source.type == EpanetChemicalSourceType::FlowPaced)
                    next[index] += source.strength_si;
            }
            next[index] = std::max(0.0, next[index]);
            maximum_change = std::max(maximum_change,
                                      std::abs(next[index] - values[index]));
        }
        values.swap(next);
        if (maximum_change < 1.0e-12)
            break;
    }
    node_concentration_kgm3_.swap(values);
}

double EpanetChemicalModel::outlet_concentration(
    std::size_t link_index, double flow_m3s) const {
    const std::deque<Segment> &pipe = segments_[link_index];
    if (pipe.empty()) {
        const std::size_t upstream = flow_m3s >= 0.0
            ? links_[link_index].from_node : links_[link_index].to_node;
        return node_concentration_kgm3_[upstream];
    }
    return flow_m3s >= 0.0 ? pipe.back().concentration_kgm3
                           : pipe.front().concentration_kgm3;
}

void EpanetChemicalModel::displace(
    std::size_t link_index, double flow_m3s, double dt_s,
    double inlet_concentration_kgm3) {
    std::deque<Segment> &pipe = segments_[link_index];
    const double volume = std::min(std::abs(flow_m3s) * dt_s,
                                   links_[link_index].volume_m3);
    double remaining = volume;
    while (remaining > kVolumeEpsilon && !pipe.empty()) {
        Segment &outlet = flow_m3s > 0.0 ? pipe.back() : pipe.front();
        const double removed = std::min(remaining, outlet.volume_m3);
        outlet.volume_m3 -= removed;
        remaining -= removed;
        if (outlet.volume_m3 <= kVolumeEpsilon) {
            if (flow_m3s > 0.0) pipe.pop_back(); else pipe.pop_front();
        }
    }
    if (volume <= kVolumeEpsilon)
        return;
    const Segment inlet{volume, inlet_concentration_kgm3};
    if (flow_m3s > 0.0) pipe.push_front(inlet); else pipe.push_back(inlet);
}

std::vector<double> EpanetChemicalModel::link_average_concentration_kgm3(
    const std::vector<double> &link_flows_m3s) const {
    if (link_flows_m3s.size() != links_.size())
        throw std::invalid_argument("Chemical flow vector does not match the link index.");
    std::vector<double> result;
    result.reserve(links_.size());
    for (std::size_t index = 0; index < links_.size(); ++index) {
        if (segments_[index].empty()) {
            const std::size_t upstream = link_flows_m3s[index] >= 0.0
                ? links_[index].from_node : links_[index].to_node;
            result.push_back(node_concentration_kgm3_[upstream]);
            continue;
        }
        double mass = 0.0;
        double volume = 0.0;
        for (const Segment &segment : segments_[index]) {
            mass += segment.volume_m3 * segment.concentration_kgm3;
            volume += segment.volume_m3;
        }
        result.push_back(volume > kVolumeEpsilon ? mass / volume : 0.0);
    }
    return result;
}
