#include "StaciSession.h"

#include "epanet_steady_quality.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kFlowEpsilon = 1.0e-12;

bool is_pipe(const Agelem &link) { return link.GetType() == "Cso"; }

bool is_fixed_quality_boundary(Agelem &link) {
    return link.Get_Csp_db() == 1 &&
           (link.GetType() == "Vegakna" || link.GetType() == "KonstNyomas");
}

std::string unsupported(const char *kind, const std::string &property,
                        const std::string &id,
                        const std::string &type = std::string()) {
    std::string message = "Unsupported ";
    message += kind;
    message += " property '" + property + "' for '" + id + "'";
    if (!type.empty())
        message += " (STACI type " + type + ")";
    return message + ".";
}

} // namespace

StaciSession::StaciSession(const std::string &network_file)
    : system_(std::make_unique<Staci>(network_file, true)) {
    system_->Set_debug_level(0);
    system_->build_system();
    for (Csomopont *item : system_->cspok) {
        const std::string id = item->Get_nev();
        if (!nodes_.emplace(id, item).second)
            throw std::runtime_error("Duplicate STACI node id: " + id);
        info_.node_ids.push_back(id);
    }
    for (Agelem *item : system_->agelemek) {
        const std::string id = item->Get_nev();
        if (!links_.emplace(id, item).second)
            throw std::runtime_error("Duplicate STACI link id: " + id);
        info_.link_ids.push_back(id);
        info_.link_types.emplace_back(item->GetType());
        info_.link_from_nodes.push_back(item->Get_Cspe_Nev());
        info_.link_to_nodes.push_back(
            item->Get_Csp_db() == 2 ? item->Get_Cspv_Nev() : std::string());
    }
}

Csomopont &StaciSession::node(const std::string &id) const {
    const auto found = nodes_.find(id);
    if (found == nodes_.end())
        throw std::invalid_argument("Unknown STACI node id: " + id);
    return *found->second;
}

Agelem &StaciSession::link(const std::string &id) const {
    const auto found = links_.find(id);
    if (found == links_.end())
        throw std::invalid_argument("Unknown STACI link id: " + id);
    return *found->second;
}

double StaciSession::value_for(const std::vector<double> &values,
                               std::size_t index, std::size_t expected) {
    if (values.size() != 1 && values.size() != expected)
        throw std::invalid_argument(
            "Property values must be scalar or match the number of element ids.");
    const double value = values.size() == 1 ? values.front() : values[index];
    if (!std::isfinite(value))
        throw std::invalid_argument("STACI property values must be finite.");
    return value;
}

void StaciSession::reset_hydraulic_state() {
    system_->ini();
    solved_ = false;
}

StaciHydraulicStatus StaciSession::solve_hydraulics() {
    bool converged = system_->solve_system();
    if (!converged) {
        system_->ini();
        converged = system_->solve_system();
    }
    solved_ = converged;
    return StaciHydraulicStatus{converged};
}

void StaciSession::require_solved() const {
    if (!solved_)
        throw std::runtime_error(
            "No converged hydraulic state is available. Call solveHydraulics first.");
}

std::vector<double> StaciSession::get_node_property(
    const std::vector<std::string> &ids, const std::string &property) const {
    std::vector<double> result;
    result.reserve(ids.size());
    for (const std::string &id : ids) {
        Csomopont &item = node(id);
        if (property == "elevation_m")
            result.push_back(item.Get_h());
        else if (property == "pressure_head_m")
            result.push_back(item.Get_p());
        else if (property == "total_head_m")
            result.push_back(item.Get_h() + item.Get_p());
        else if (property == "demand_m3s")
            result.push_back(item.Get_dprop("demand") / 3600.0);
        else if (property == "water_age_s")
            result.push_back(item.Get_dprop("tt"));
        else if (property == "concentration_kgm3")
            result.push_back(item.Get_dprop("concentration"));
        else if (property == "source_concentration_kgm3")
            result.push_back(item.Get_dprop("cl_be"));
        else
            throw std::invalid_argument(unsupported("node", property, id));
    }
    return result;
}

void StaciSession::set_node_property(const std::vector<std::string> &ids,
                                     const std::string &property,
                                     const std::vector<double> &values) {
    for (std::size_t index = 0; index < ids.size(); ++index) {
        Csomopont &item = node(ids[index]);
        const double value = value_for(values, index, ids.size());
        if (property == "demand_m3s")
            item.Set_dprop("demand", value * 3600.0);
        else if (property == "pressure_head_m")
            item.Set_p(value);
        else if (property == "water_age_s")
            item.Set_dprop("tt", value);
        else if (property == "concentration_kgm3")
            item.Set_dprop("concentration", value);
        else if (property == "source_concentration_kgm3")
            item.Set_dprop("cl_be", value);
        else
            throw std::invalid_argument(unsupported("node", property, ids[index]));
    }
    solved_ = false;
}

std::vector<double> StaciSession::get_link_property(
    const std::vector<std::string> &ids, const std::string &property) const {
    std::vector<double> result;
    result.reserve(ids.size());
    for (const std::string &id : ids) {
        Agelem &item = link(id);
        const std::string type(item.GetType());
        if (property == "flow_rate_m3s")
            result.push_back(item.Get_Q());
        else if (property == "mass_flow_rate_kgs")
            result.push_back(item.Get_mp());
        else if (property == "velocity_ms")
            result.push_back(item.Get_v());
        else if (property == "status")
            result.push_back(item.Is_enabled() ? 1.0 : 0.0);
        else if (property == "water_age_start_s")
            result.push_back(item.Get_tt_start());
        else if (property == "water_age_end_s")
            result.push_back(item.Get_tt_end());
        else if (property == "diameter_m" && is_pipe(item))
            result.push_back(item.Get_dprop("diameter"));
        else if (property == "length_m" && is_pipe(item))
            result.push_back(item.Get_dprop("length"));
        else if (property == "friction_coeff" && is_pipe(item))
            result.push_back(item.Get_dprop("friction_coeff"));
        else if (property == "minor_loss" && is_pipe(item))
            result.push_back(item.Get_dprop("minor_loss"));
        else if (property == "bulk_reaction_per_s" && is_pipe(item))
            result.push_back(item.Get_dprop("cl_k"));
        else if (property == "wall_reaction_m_s" && is_pipe(item))
            result.push_back(item.Get_dprop("cl_w"));
        else if (property == "speed" &&
                 (type == "Szivattyu" || type == "EpanetPowerPump"))
            result.push_back(item.Get_dprop("speed"));
        else if (property == "power_w" && type == "EpanetPowerPump")
            result.push_back(item.Get_dprop("power"));
        else if (property == "position" && type == "JelleggorbesFojtas")
            result.push_back(item.Get_dprop("position"));
        else if (property == "tcv_setting" && type == "JelleggorbesFojtas")
            result.push_back(item.Get_dprop("tcv_setting"));
        else if (property == "water_level_m" && type == "Vegakna")
            result.push_back(item.Get_dprop("water_level"));
        else if (property == "bottom_level_m" && type == "Vegakna")
            result.push_back(item.Get_dprop("bottom_level"));
        else if (property == "boundary_head_m" && type == "KonstNyomas")
            result.push_back(item.Get_dprop("head"));
        else
            throw std::invalid_argument(unsupported("link", property, id, type));
    }
    return result;
}

void StaciSession::set_link_property(const std::vector<std::string> &ids,
                                     const std::string &property,
                                     const std::vector<double> &values) {
    for (std::size_t index = 0; index < ids.size(); ++index) {
        Agelem &item = link(ids[index]);
        const std::string type(item.GetType());
        const double value = value_for(values, index, ids.size());
        if (property == "status")
            item.Set_enabled(value != 0.0);
        else if (property == "diameter_m" && is_pipe(item)) {
            if (value <= 0.0)
                throw std::invalid_argument("Pipe diameter must be positive.");
            item.Set_dprop("diameter", value);
        } else if (property == "friction_coeff" && is_pipe(item))
            item.Set_dprop("friction_coeff", value);
        else if (property == "minor_loss" && is_pipe(item))
            item.Set_dprop("minor_loss", value);
        else if (property == "speed" &&
                 (type == "Szivattyu" || type == "EpanetPowerPump"))
            item.Set_dprop("speed", value);
        else if (property == "power_w" && type == "EpanetPowerPump")
            item.Set_dprop("power", value);
        else if (property == "position" && type == "JelleggorbesFojtas")
            item.Set_dprop("position", value);
        else if (property == "tcv_setting" && type == "JelleggorbesFojtas")
            item.Set_dprop("tcv_setting", value);
        else if (property == "water_level_m" && type == "Vegakna")
            item.Set_dprop("water_level", value);
        else if (property == "bottom_level_m" && type == "Vegakna")
            item.Set_dprop("bottom_level", value);
        else if (property == "boundary_head_m" && type == "KonstNyomas")
            item.Set_dprop("head", value);
        else
            throw std::invalid_argument(unsupported("link", property, ids[index], type));
    }
    solved_ = false;
}

StaciNodeResults StaciSession::node_results() const {
    require_solved();
    StaciNodeResults result;
    result.ids = info_.node_ids;
    for (const std::string &id : info_.node_ids) {
        Csomopont &item = node(id);
        result.elevation_m.push_back(item.Get_h());
        result.pressure_head_m.push_back(item.Get_p());
        result.total_head_m.push_back(item.Get_h() + item.Get_p());
        result.demand_m3s.push_back(item.Get_dprop("demand") / 3600.0);
    }
    return result;
}

StaciLinkResults StaciSession::link_results() const {
    require_solved();
    StaciLinkResults result;
    result.ids = info_.link_ids;
    result.types = info_.link_types;
    result.from_nodes = info_.link_from_nodes;
    result.to_nodes = info_.link_to_nodes;
    for (std::size_t index = 0; index < info_.link_ids.size(); ++index) {
        Agelem &item = link(info_.link_ids[index]);
        result.flow_rate_m3s.push_back(item.Get_Q());
        result.velocity_ms.push_back(item.Get_v());
        result.status.push_back(item.Is_enabled() ? 1.0 : 0.0);
        double headloss = 0.0;
        if (item.Get_Csp_db() == 2) {
            Csomopont &from = node(item.Get_Cspe_Nev());
            Csomopont &to = node(item.Get_Cspv_Nev());
            headloss = from.Get_h() + from.Get_p() - to.Get_h() - to.Get_p();
        }
        result.headloss_m.push_back(headloss);
    }
    return result;
}

StaciSteadyQualityResults StaciSession::solve_steady_quality(
    bool solve_age, bool solve_chemical) const {
    require_solved();
    if (!solve_age && !solve_chemical)
        throw std::invalid_argument("Select steady water age and/or chemical quality.");

    std::unordered_map<std::string, std::size_t> node_index;
    std::vector<EpanetSteadyQualityNode> quality_nodes(info_.node_ids.size());
    std::vector<double> incoming(info_.node_ids.size(), 0.0);
    for (std::size_t index = 0; index < info_.node_ids.size(); ++index) {
        node_index[info_.node_ids[index]] = index;
        Csomopont &item = node(info_.node_ids[index]);
        quality_nodes[index].id = info_.node_ids[index];
        quality_nodes[index].external_inflow_m3s =
            std::max(0.0, -item.Get_dprop("demand") / 3600.0);
        quality_nodes[index].external_concentration_kgm3 =
            item.Get_dprop("cl_be");
    }

    for (const std::string &id : info_.link_ids) {
        Agelem &item = link(id);
        if (!is_fixed_quality_boundary(item))
            continue;
        const auto found = node_index.find(item.Get_Cspe_Nev());
        if (found == node_index.end())
            continue;
        EpanetSteadyQualityNode &boundary = quality_nodes[found->second];
        boundary.fixed_age = true;
        boundary.fixed_chemical = true;
        boundary.fixed_concentration_kgm3 =
            node(item.Get_Cspe_Nev()).Get_dprop("cl_be");
    }

    std::vector<EpanetSteadyQualityLink> quality_links;
    std::vector<std::string> quality_link_ids;
    for (const std::string &id : info_.link_ids) {
        Agelem &item = link(id);
        if (item.Get_Csp_db() != 2)
            continue;
        const auto from = node_index.find(item.Get_Cspe_Nev());
        const auto to = node_index.find(item.Get_Cspv_Nev());
        if (from == node_index.end() || to == node_index.end())
            continue;
        double volume = 0.0;
        double reaction = 0.0;
        if (is_pipe(item)) {
            const double diameter = item.Get_dprop("diameter");
            const double length = item.Get_dprop("length");
            volume = kPi * diameter * diameter * length / 4.0;
            reaction = item.Get_dprop("cl_k") +
                       4.0 * item.Get_dprop("cl_w") / diameter;
        }
        quality_links.push_back(EpanetSteadyQualityLink{
            id, from->second, to->second, item.Get_Q(), volume, reaction});
        quality_link_ids.push_back(id);
        if (item.Get_Q() > kFlowEpsilon)
            incoming[to->second] += item.Get_Q();
        else if (item.Get_Q() < -kFlowEpsilon)
            incoming[from->second] -= item.Get_Q();
    }

    for (std::size_t index = 0; index < quality_nodes.size(); ++index) {
        Csomopont &item = node(info_.node_ids[index]);
        if (!quality_nodes[index].fixed_age &&
            std::abs(item.Get_dprop("demand")) <= kFlowEpsilon &&
            incoming[index] <= kFlowEpsilon) {
            quality_nodes[index].fixed_age = true;
            quality_nodes[index].fixed_chemical = true;
            quality_nodes[index].fixed_concentration_kgm3 =
                item.Get_dprop("cl_be");
        }
    }

    EpanetSteadyQualityModel quality(std::move(quality_nodes),
                                     std::move(quality_links));
    const EpanetSteadyQualityResult solved =
        quality.solve(solve_age, solve_chemical);
    StaciSteadyQualityResults result;
    result.node_ids = info_.node_ids;
    result.link_ids = std::move(quality_link_ids);
    result.node_water_age_s = solved.node_age_s;
    result.link_average_water_age_s = solved.link_average_age_s;
    result.link_travel_time_s = solved.link_travel_time_s;
    result.node_concentration_kgm3 = solved.node_concentration_kgm3;
    result.link_average_concentration_kgm3 =
        solved.link_average_concentration_kgm3;
    return result;
}

StaciSensitivityResults StaciSession::hydraulic_sensitivity(
    const std::string &element_id, const std::string &property) {
    require_solved();
    std::string legacy_property;
    bool node_parameter = false;
    if (property == "diameter_m")
        legacy_property = "diameter";
    else if (property == "friction_coeff")
        legacy_property = "friction_coeff";
    else if (property == "demand_m3s") {
        legacy_property = "demand";
        node_parameter = true;
    } else {
        throw std::invalid_argument(
            "Hydraulic sensitivity supports diameter_m, friction_coeff and demand_m3s.");
    }
    if (node_parameter)
        (void)node(element_id);
    else {
        Agelem &selected = link(element_id);
        if (!is_pipe(selected))
            throw std::invalid_argument(
                "Diameter and friction sensitivity require a Cso element.");
    }

    system_->element_ID = element_id;
    system_->property_ID = legacy_property;
    system_->Compute_dxdmu();
    const std::vector<double> &derivative = system_->Get_dxdmu();
    if (derivative.size() != info_.link_ids.size() + info_.node_ids.size())
        throw std::runtime_error("STACI returned an invalid sensitivity vector size.");

    StaciSensitivityResults result;
    result.parameter_element_id = element_id;
    result.parameter_property = property;
    result.link_ids = info_.link_ids;
    result.node_ids = info_.node_ids;
    for (std::size_t index = 0; index < info_.link_ids.size(); ++index) {
        Agelem &item = link(info_.link_ids[index]);
        double value = derivative[index] / item.Get_ro();
        if (node_parameter)
            value *= node(element_id).Get_dprop("ro");
        result.link_flow_rate_derivative.push_back(value);
    }
    for (std::size_t index = 0; index < info_.node_ids.size(); ++index) {
        double value = derivative[info_.link_ids.size() + index];
        if (node_parameter)
            value *= node(element_id).Get_dprop("ro");
        result.node_pressure_head_derivative.push_back(value);
    }
    return result;
}
