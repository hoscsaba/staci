#include "epanet_extended_simulation.h"

#include "Agelem.h"
#include "Csomopont.h"
#include "Staci.h"
#include "eps_result_writer.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {

struct Record {
    std::size_t line_number;
    std::vector<std::string> fields;
};

struct DemandComponent {
    double base_m3_per_hour;
    std::string pattern_id;
};

struct TankState {
    std::string id;
    double level_m;
    double min_level_m;
    double max_level_m;
    double area_m2;
    double min_volume_m3;
    Agelem *boundary;
};

struct ReservoirState {
    std::string id;
    double base_head_m;
    std::string pattern_id;
    Agelem *boundary;
};

struct LevelControl {
    std::string link_id;
    bool enable;
    std::string tank_id;
    bool below;
    double threshold_m;
    std::size_t line_number;
};

std::string trim(const std::string &value) {
    const std::string whitespace = " \t\r\n";
    const std::size_t begin = value.find_first_not_of(whitespace);
    if (begin == std::string::npos)
        return "";
    const std::size_t end = value.find_last_not_of(whitespace);
    return value.substr(begin, end - begin + 1);
}

std::string upper(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return value;
}

std::vector<std::string> split_fields(const std::string &line) {
    std::istringstream input(line);
    std::vector<std::string> result;
    std::string field;
    while (input >> field)
        result.push_back(field);
    return result;
}

double parse_number(const std::string &value, const std::string &context) {
    char *end = nullptr;
    const double result = std::strtod(value.c_str(), &end);
    if (end == value.c_str() || *end != '\0')
        throw std::runtime_error("Invalid number '" + value + "' in " + context + ".");
    return result;
}

long long parse_time_seconds(const std::string &value) {
    const std::size_t colon = value.find(':');
    if (colon == std::string::npos)
        return static_cast<long long>(std::llround(parse_number(value, "[TIMES]") * 3600.0));
    const double hours = parse_number(value.substr(0, colon), "[TIMES]");
    const double minutes = parse_number(value.substr(colon + 1), "[TIMES]");
    return static_cast<long long>(std::llround(hours * 3600.0 + minutes * 60.0));
}

bool us_units(const std::string &flow_units) {
    const std::string unit = upper(flow_units);
    return unit == "CFS" || unit == "GPM" || unit == "MGD" ||
           unit == "IMGD" || unit == "AFD";
}

double flow_to_m3_per_hour(double value, const std::string &unit) {
    const std::string normalized = upper(unit);
    if (normalized == "CFS") return value * 101.9406477312;
    if (normalized == "GPM") return value * 0.22712470704;
    if (normalized == "MGD") return value * 157.725491;
    if (normalized == "IMGD") return value * 189.420142;
    if (normalized == "AFD") return value * 51.3950766;
    if (normalized == "LPS") return value * 3.6;
    if (normalized == "LPM") return value * 0.06;
    if (normalized == "MLD") return value * 41.6666666667;
    if (normalized == "CMH") return value;
    if (normalized == "CMD") return value / 24.0;
    return value * 3.6;
}

double length_to_metres(double value, bool is_us) {
    return is_us ? value * 0.3048 : value;
}

double diameter_to_metres(double value, bool is_us) {
    return is_us ? value * 0.3048 : value;
}

double volume_to_m3(double value, bool is_us) {
    return is_us ? value * 0.028316846592 : value;
}

std::string csv(const std::string &value) {
    if (value.find_first_of(",\"\r\n") == std::string::npos)
        return value;
    std::string escaped = "\"";
    for (char character : value) {
        escaped += character;
        if (character == '"')
            escaped += '"';
    }
    escaped += '"';
    return escaped;
}

class SimulationModel {
public:
    explicit SimulationModel(const std::string &filename)
        : filename_(filename), flow_units_("LPS"), demand_multiplier_(1.0),
          duration_s_(0), hydraulic_step_s_(3600), pattern_step_s_(3600),
          pattern_start_s_(0), report_step_s_(3600), report_start_s_(0) {
        parse_file();
        parse_options();
        parse_times();
        parse_patterns();
        parse_demands();
        parse_boundaries();
        parse_controls();
        parse_status();
    }

    void run(Staci &system, const std::string &prefix) {
        if (duration_s_ < 0 || hydraulic_step_s_ <= 0 || pattern_step_s_ <= 0 ||
            report_step_s_ <= 0)
            throw std::runtime_error("EPANET duration and time steps must be positive.");

        std::map<std::string, Csomopont *> nodes;
        for (Csomopont *node : system.cspok)
            nodes[node->Get_nev()] = node;
        std::map<std::string, Agelem *> links;
        for (Agelem *link : system.agelemek)
            links[link->Get_nev()] = link;

        std::vector<TankState> tanks;
        for (TankState tank : tank_definitions_) {
            const std::string boundary_id = "EPANET_TANK_" + tank.id;
            const auto boundary = links.find(boundary_id);
            if (boundary == links.end()) {
                warn("TANKS", tank.id, 0, "Imported STACI tank boundary was not found.");
                continue;
            }
            tank.boundary = boundary->second;
            tanks.push_back(tank);
        }

        std::vector<ReservoirState> reservoirs;
        for (ReservoirState reservoir : reservoir_definitions_) {
            const std::string boundary_id = "EPANET_RESERVOIR_" + reservoir.id;
            const auto boundary = links.find(boundary_id);
            if (boundary == links.end()) {
                warn("RESERVOIRS", reservoir.id, 0,
                     "Imported STACI reservoir boundary was not found.");
                continue;
            }
            reservoir.boundary = boundary->second;
            reservoirs.push_back(reservoir);
        }

        for (const auto &status : initial_status_) {
            const auto link = links.find(status.first);
            if (link != links.end() && is_pump(link->second))
                link->second->Set_enabled(status.second);
            else if (link != links.end())
                warn("STATUS", status.first, 0,
                     "EPS status switching is currently implemented only for pumps.");
        }

        std::vector<EpsNodeInfo> result_nodes;
        std::map<std::string, std::uint32_t> node_indices;
        std::map<std::string, std::string> node_types;
        std::map<std::string, double> node_elevations;
        std::map<std::string, std::pair<double, double> > coordinates;
        for (const Record &record : sections_["JUNCTIONS"])
            if (record.fields.size() >= 2) {
                node_types[record.fields[0]] = "JUNCTION";
                node_elevations[record.fields[0]] = length_to_metres(
                    parse_number(record.fields[1], "[JUNCTIONS] elevation"), us_units(flow_units_));
            }
        for (const Record &record : sections_["RESERVOIRS"])
            if (!record.fields.empty()) {
                node_types[record.fields[0]] = "RESERVOIR";
                node_elevations[record.fields[0]] = 0.0;
            }
        for (const Record &record : sections_["TANKS"])
            if (record.fields.size() >= 2) {
                node_types[record.fields[0]] = "TANK";
                node_elevations[record.fields[0]] = length_to_metres(
                    parse_number(record.fields[1], "[TANKS] elevation"), us_units(flow_units_));
            }
        for (const Record &record : sections_["COORDINATES"])
            if (record.fields.size() >= 3)
                coordinates[record.fields[0]] = std::make_pair(
                    length_to_metres(parse_number(record.fields[1], "[COORDINATES] x"),
                                     us_units(flow_units_)),
                    length_to_metres(parse_number(record.fields[2], "[COORDINATES] y"),
                                     us_units(flow_units_)));
        const double missing = std::numeric_limits<double>::quiet_NaN();
        for (std::size_t index = 0; index < system.cspok.size(); ++index) {
            Csomopont *node = system.cspok[index];
            const std::string id = node->Get_nev();
            node_indices[id] = static_cast<std::uint32_t>(index);
            const auto type = node_types.find(id);
            const auto elevation = node_elevations.find(id);
            const auto coordinate = coordinates.find(id);
            result_nodes.push_back(EpsNodeInfo{
                id,
                type == node_types.end() ? "JUNCTION" : type->second,
                elevation == node_elevations.end() ? node->Get_h() : elevation->second,
                coordinate == coordinates.end() ? missing : coordinate->second.first,
                coordinate == coordinates.end() ? missing : coordinate->second.second});
        }

        std::vector<Agelem *> result_link_objects;
        std::vector<EpsLinkInfo> result_links;
        for (Agelem *link : system.agelemek) {
            if (link->Get_nev().rfind("EPANET_TANK_", 0) == 0 ||
                link->Get_nev().rfind("EPANET_RESERVOIR_", 0) == 0)
                continue;
            const auto from = node_indices.find(link->Get_Cspe_Nev());
            const auto to = node_indices.find(link->Get_Cspv_Nev());
            if (from == node_indices.end() || to == node_indices.end())
                continue;
            const bool pipe = link->GetType() == "Cso";
            std::string output_type = upper(link->GetType());
            if (pipe)
                output_type = "PIPE";
            else if (is_pump(link))
                output_type = "PUMP";
            result_link_objects.push_back(link);
            result_links.push_back(EpsLinkInfo{
                link->Get_nev(), output_type, from->second, to->second,
                pipe ? link->Get_dprop("length") : missing,
                pipe ? link->Get_dprop("diameter") : missing});
        }

        std::vector<EpsTankInfo> result_tanks;
        for (const TankState &tank : tanks) {
            const auto node = node_indices.find(tank.id);
            if (node != node_indices.end())
                result_tanks.push_back(EpsTankInfo{
                    tank.id, node->second, tank.min_level_m, tank.max_level_m,
                    std::sqrt(4.0 * tank.area_m2 / 3.14159265358979323846)});
        }

        std::ofstream node_output((prefix + "-nodes.csv").c_str(), std::ios::trunc);
        std::ofstream link_output((prefix + "-links.csv").c_str(), std::ios::trunc);
        std::ofstream tank_output((prefix + "-tanks.csv").c_str(), std::ios::trunc);
        std::ofstream summary_output((prefix + "-summary.csv").c_str(), std::ios::trunc);
        if (!node_output || !link_output || !tank_output || !summary_output)
            throw std::runtime_error("Cannot create EPS CSV output files with prefix: " + prefix);

        node_output << "time_seconds,node_id,elevation_m,pressure_head_m,total_head_m,demand_m3s,converged\n";
        link_output << "time_seconds,link_id,type,node_from,node_to,flow_m3s,velocity_mps,headloss_m,status,converged\n";
        tank_output << "time_seconds,tank_id,level_m,volume_m3,inflow_m3s,min_level_m,max_level_m,converged\n";
        node_output << std::setprecision(15);
        link_output << std::setprecision(15);
        tank_output << std::setprecision(15);

        const int previous_debug_level = system.Get_debug_level();
        system.Set_debug_level(0);
        system.ini();

        long long simulation_step_s = std::gcd(hydraulic_step_s_, pattern_step_s_);
        simulation_step_s = std::gcd(simulation_step_s, report_step_s_);
        if (pattern_start_s_ > 0)
            simulation_step_s = std::gcd(simulation_step_s, pattern_start_s_);
        if (report_start_s_ > 0)
            simulation_step_s = std::gcd(simulation_step_s, report_start_s_);
        if (simulation_step_s <= 0)
            simulation_step_s = hydraulic_step_s_;

        EpsOutputMetadata output_metadata{
            filename_, duration_s_, hydraulic_step_s_, simulation_step_s,
            pattern_step_s_, report_step_s_};
        EpsResultWriter result_writer(prefix, result_nodes, result_links,
                                      result_tanks, output_metadata, 16);
        if (!result_writer.hdf5_enabled())
            warn("OUTPUT", "HDF5", 0,
                 "HDF5 support was not available when STACI was built; only SI CSV and JSON files are written.");

        std::size_t state_count = 0;
        std::size_t failed_count = 0;
        long long time_s = 0;
        while (time_s <= duration_s_) {
            apply_demands(nodes, time_s);
            apply_reservoir_heads(reservoirs, time_s);
            for (TankState &tank : tanks) {
                tank.boundary->Set_dprop("water_level", tank.level_m);
                tank.boundary->Set_enabled(true);
            }
            apply_controls(links, tanks);

            bool converged = system.solve_system();
            bool constrained_tank = false;
            for (TankState &tank : tanks) {
                const double flow = tank.boundary->Get_Q();
                if ((tank.level_m <= tank.min_level_m + 1.0e-10 && flow < 0.0) ||
                    (tank.level_m >= tank.max_level_m - 1.0e-10 && flow > 0.0)) {
                    tank.boundary->Set_enabled(false);
                    constrained_tank = true;
                }
            }
            if (constrained_tank)
                converged = system.solve_system() && converged;
            ++state_count;
            if (!converged)
                ++failed_count;

            if (should_report(time_s)) {
                const EpsResultFrame frame = collect_frame(
                    system, result_nodes, result_links, result_link_objects,
                    tanks, time_s, converged);
                result_writer.append(frame);
                write_state(node_output, link_output, tank_output, result_nodes,
                            result_links, result_tanks, frame);
            }

            if (time_s == duration_s_)
                break;
            const long long step_s = std::min(simulation_step_s, duration_s_ - time_s);
            update_tanks(tanks, static_cast<double>(step_s));
            time_s += step_s;
        }
        system.Set_debug_level(previous_debug_level);

        result_writer.finish(state_count, failed_count, warnings_);

        summary_output << "property,value\n"
                       << "input_file," << csv(filename_) << "\n"
                       << "duration_seconds," << duration_s_ << "\n"
                       << "hydraulic_timestep_seconds," << hydraulic_step_s_ << "\n"
                       << "effective_simulation_timestep_seconds," << simulation_step_s << "\n"
                       << "pattern_timestep_seconds," << pattern_step_s_ << "\n"
                       << "report_timestep_seconds," << report_step_s_ << "\n"
                       << "hydraulic_states," << state_count << "\n"
                       << "failed_states," << failed_count << "\n"
                       << "warnings," << warnings_.size() << "\n";

        for (const std::string &warning : warnings_)
            std::cerr << warning << '\n';
        std::cout << "\nEPANET EPS complete: " << state_count << " hydraulic states, "
                  << failed_count << " failed states.\n";
        if (result_writer.hdf5_enabled())
            std::cout << "Results: " << prefix << ".h5, " << prefix
                      << ".meta.json and SI CSV files.\n";
        else
            std::cout << "Results: " << prefix
                      << ".meta.json and SI CSV files (HDF5 unavailable).\n";
    }

private:
    std::string filename_;
    std::map<std::string, std::vector<Record> > sections_;
    std::map<std::string, std::vector<double> > patterns_;
    std::map<std::string, std::vector<DemandComponent> > demands_;
    std::vector<TankState> tank_definitions_;
    std::vector<ReservoirState> reservoir_definitions_;
    std::vector<LevelControl> controls_;
    std::map<std::string, bool> initial_status_;
    std::vector<std::string> warnings_;
    std::set<std::string> warning_keys_;
    std::set<std::string> missing_patterns_;
    std::string flow_units_;
    std::string default_pattern_;
    double demand_multiplier_;
    long long duration_s_;
    long long hydraulic_step_s_;
    long long pattern_step_s_;
    long long pattern_start_s_;
    long long report_step_s_;
    long long report_start_s_;

    void parse_file() {
        std::ifstream input(filename_.c_str());
        if (!input)
            throw std::runtime_error("Cannot open EPANET input file: " + filename_);
        std::string section;
        std::string line;
        std::size_t line_number = 0;
        while (std::getline(input, line)) {
            ++line_number;
            const std::size_t comment = line.find(';');
            const std::string data = trim(line.substr(0, comment));
            if (data.empty())
                continue;
            if (data.front() == '[' && data.back() == ']') {
                section = upper(trim(data.substr(1, data.size() - 2)));
                continue;
            }
            sections_[section].push_back(Record{line_number, split_fields(data)});
        }
    }

    const std::vector<Record> &records(const std::string &section) const {
        static const std::vector<Record> empty;
        const auto found = sections_.find(section);
        return found == sections_.end() ? empty : found->second;
    }

    void parse_options() {
        for (const Record &record : records("OPTIONS")) {
            if (record.fields.size() < 2)
                continue;
            const std::string key = upper(record.fields[0]);
            if (key == "UNITS")
                flow_units_ = upper(record.fields[1]);
            else if (key == "PATTERN")
                default_pattern_ = record.fields[1];
            else if (record.fields.size() > 2 && key == "DEMAND" &&
                     upper(record.fields[1]) == "MULTIPLIER")
                demand_multiplier_ = parse_number(record.fields[2], "DEMAND MULTIPLIER");
        }
    }

    void parse_times() {
        for (const Record &record : records("TIMES")) {
            if (record.fields.size() < 2)
                continue;
            const std::string first = upper(record.fields[0]);
            if (first == "DURATION")
                duration_s_ = parse_time_seconds(record.fields[1]);
            else if (record.fields.size() > 2 && first == "HYDRAULIC" &&
                     upper(record.fields[1]) == "TIMESTEP")
                hydraulic_step_s_ = parse_time_seconds(record.fields[2]);
            else if (record.fields.size() > 2 && first == "PATTERN" &&
                     upper(record.fields[1]) == "TIMESTEP")
                pattern_step_s_ = parse_time_seconds(record.fields[2]);
            else if (record.fields.size() > 2 && first == "PATTERN" &&
                     upper(record.fields[1]) == "START")
                pattern_start_s_ = parse_time_seconds(record.fields[2]);
            else if (record.fields.size() > 2 && first == "REPORT" &&
                     upper(record.fields[1]) == "TIMESTEP")
                report_step_s_ = parse_time_seconds(record.fields[2]);
            else if (record.fields.size() > 2 && first == "REPORT" &&
                     upper(record.fields[1]) == "START")
                report_start_s_ = parse_time_seconds(record.fields[2]);
        }
    }

    void parse_patterns() {
        for (const Record &record : records("PATTERNS")) {
            if (record.fields.size() < 2)
                continue;
            std::vector<double> &pattern = patterns_[record.fields[0]];
            for (std::size_t i = 1; i < record.fields.size(); ++i)
                pattern.push_back(parse_number(record.fields[i], "[PATTERNS]"));
        }
    }

    void parse_demands() {
        for (const Record &record : records("JUNCTIONS")) {
            if (record.fields.size() < 2)
                continue;
            const double base = record.fields.size() > 2
                ? parse_number(record.fields[2], "[JUNCTIONS]") : 0.0;
            const std::string pattern = record.fields.size() > 3
                ? record.fields[3] : default_pattern_;
            demands_[record.fields[0]].push_back(
                DemandComponent{flow_to_m3_per_hour(base, flow_units_), pattern});
        }
        for (const Record &record : records("DEMANDS")) {
            if (record.fields.size() < 2)
                continue;
            const std::string pattern = record.fields.size() > 2
                ? record.fields[2] : default_pattern_;
            demands_[record.fields[0]].push_back(
                DemandComponent{flow_to_m3_per_hour(
                    parse_number(record.fields[1], "[DEMANDS]"), flow_units_), pattern});
        }
    }

    void parse_boundaries() {
        const bool is_us = us_units(flow_units_);
        for (const Record &record : records("TANKS")) {
            if (record.fields.size() < 7)
                continue;
            const double diameter = diameter_to_metres(
                parse_number(record.fields[5], "[TANKS] diameter"), is_us);
            tank_definitions_.push_back(TankState{
                record.fields[0],
                length_to_metres(parse_number(record.fields[2], "[TANKS] initial level"), is_us),
                length_to_metres(parse_number(record.fields[3], "[TANKS] minimum level"), is_us),
                length_to_metres(parse_number(record.fields[4], "[TANKS] maximum level"), is_us),
                3.14159265358979323846 * diameter * diameter / 4.0,
                volume_to_m3(parse_number(record.fields[6], "[TANKS] minimum volume"), is_us),
                nullptr});
            if (record.fields.size() > 7 && !record.fields[7].empty())
                warn("TANKS", record.fields[0], record.line_number,
                     "Tank volume curves are approximated with the cylindrical diameter.");
        }
        for (const Record &record : records("RESERVOIRS")) {
            if (record.fields.size() < 2)
                continue;
            reservoir_definitions_.push_back(ReservoirState{
                record.fields[0],
                length_to_metres(parse_number(record.fields[1], "[RESERVOIRS]"), is_us),
                record.fields.size() > 2 ? record.fields[2] : "",
                nullptr});
        }
    }

    void parse_controls() {
        const bool is_us = us_units(flow_units_);
        for (const Record &record : records("CONTROLS")) {
            if (record.fields.size() >= 8 && upper(record.fields[0]) == "LINK" &&
                (upper(record.fields[2]) == "OPEN" || upper(record.fields[2]) == "CLOSED") &&
                upper(record.fields[3]) == "IF" && upper(record.fields[4]) == "NODE" &&
                (upper(record.fields[6]) == "BELOW" || upper(record.fields[6]) == "ABOVE")) {
                controls_.push_back(LevelControl{
                    record.fields[1], upper(record.fields[2]) == "OPEN", record.fields[5],
                    upper(record.fields[6]) == "BELOW",
                    length_to_metres(parse_number(record.fields[7], "[CONTROLS]"), is_us),
                    record.line_number});
            } else {
                warn("CONTROLS", record.fields.size() > 1 ? record.fields[1] : "",
                     record.line_number,
                     "Only OPEN/CLOSED pump controls based on tank levels are executed in EPS mode.");
            }
        }
    }

    void parse_status() {
        for (const Record &record : records("STATUS")) {
            if (record.fields.size() < 2)
                continue;
            const std::string status = upper(record.fields[1]);
            if (status == "OPEN" || status == "CLOSED")
                initial_status_[record.fields[0]] = status == "OPEN";
            else
                warn("STATUS", record.fields[0], record.line_number,
                     "Numeric link settings are not executed in EPS mode.");
        }
    }

    double pattern_value(const std::string &pattern_id, long long time_s) {
        if (pattern_id.empty())
            return 1.0;
        const auto found = patterns_.find(pattern_id);
        if (found == patterns_.end() || found->second.empty()) {
            if (missing_patterns_.insert(pattern_id).second)
                warn("PATTERNS", pattern_id, 0, "Referenced pattern was not found; multiplier 1.0 was used.");
            return 1.0;
        }
        const long long shifted = std::max(0LL, time_s + pattern_start_s_);
        const std::size_t index = static_cast<std::size_t>(shifted / pattern_step_s_) %
                                  found->second.size();
        return found->second[index];
    }

    void apply_demands(const std::map<std::string, Csomopont *> &nodes,
                       long long time_s) {
        for (const auto &node_demands : demands_) {
            const auto node = nodes.find(node_demands.first);
            if (node == nodes.end())
                continue;
            double demand = 0.0;
            for (const DemandComponent &component : node_demands.second)
                demand += component.base_m3_per_hour *
                          pattern_value(component.pattern_id, time_s);
            node->second->Set_dprop("demand", demand * demand_multiplier_);
        }
    }

    void apply_reservoir_heads(std::vector<ReservoirState> &reservoirs,
                               long long time_s) {
        for (ReservoirState &reservoir : reservoirs)
            reservoir.boundary->Set_dprop(
                "head", reservoir.base_head_m * pattern_value(reservoir.pattern_id, time_s));
    }

    void apply_controls(const std::map<std::string, Agelem *> &links,
                        const std::vector<TankState> &tanks) {
        std::map<std::string, double> levels;
        for (const TankState &tank : tanks)
            levels[tank.id] = tank.level_m;
        for (const LevelControl &control : controls_) {
            const auto link = links.find(control.link_id);
            const auto level = levels.find(control.tank_id);
            if (link == links.end() || level == levels.end())
                continue;
            if (!is_pump(link->second)) {
                warn("CONTROLS", control.link_id, control.line_number,
                     "Only pump OPEN/CLOSED controls are currently switchable.");
                continue;
            }
            const bool condition = control.below
                ? level->second < control.threshold_m
                : level->second > control.threshold_m;
            if (condition)
                link->second->Set_enabled(control.enable);
        }
    }

    void update_tanks(std::vector<TankState> &tanks, double step_s) {
        for (TankState &tank : tanks) {
            const double next = tank.level_m + tank.boundary->Get_Q() * step_s / tank.area_m2;
            if (next < tank.min_level_m) {
                tank.level_m = tank.min_level_m;
                warn("TANKS", tank.id, 0, "Minimum level reached; level was clamped.");
            } else if (next > tank.max_level_m) {
                tank.level_m = tank.max_level_m;
                warn("TANKS", tank.id, 0, "Maximum level reached; level was clamped.");
            } else {
                tank.level_m = next;
            }
        }
    }

    bool should_report(long long time_s) const {
        return time_s >= report_start_s_ &&
               (time_s - report_start_s_) % report_step_s_ == 0;
    }

    static bool is_pump(Agelem *link) {
        return link->GetType() == "Szivattyu" || link->GetType() == "EpanetPowerPump";
    }

    EpsResultFrame collect_frame(Staci &system,
                                 const std::vector<EpsNodeInfo> &node_info,
                                 const std::vector<EpsLinkInfo> &link_info,
                                 const std::vector<Agelem *> &links,
                                 const std::vector<TankState> &tanks,
                                 long long time_s, bool converged) {
        EpsResultFrame frame;
        frame.time_s = time_s;
        frame.converged = converged;
        frame.iterations = 0; // The legacy solver does not expose this counter yet.
        std::map<std::string, double> total_head;
        for (std::size_t index = 0; index < system.cspok.size(); ++index) {
            Csomopont *node = system.cspok[index];
            const double head = node->Get_p() + node->Get_h();
            total_head[node->Get_nev()] = head;
            frame.node_head_m.push_back(head);
            frame.node_pressure_head_m.push_back(head - node_info[index].elevation_m);
            frame.node_demand_m3s.push_back(node->Get_dprop("demand") / 3600.0);
        }
        for (std::size_t index = 0; index < links.size(); ++index) {
            Agelem *link = links[index];
            const double flow = link->Get_Q();
            const double diameter = link_info[index].diameter_m;
            const double area = std::isfinite(diameter)
                ? 3.14159265358979323846 * diameter * diameter / 4.0 : 0.0;
            frame.link_flow_rate_m3s.push_back(flow);
            frame.link_velocity_ms.push_back(area > 0.0
                ? flow / area : std::numeric_limits<double>::quiet_NaN());
            frame.link_headloss_m.push_back(
                total_head[link->Get_Cspe_Nev()] - total_head[link->Get_Cspv_Nev()]);
            frame.link_status.push_back(link->Is_enabled() ? 1 : 0);
        }
        for (const TankState &tank : tanks) {
            frame.tank_level_m.push_back(tank.level_m);
            frame.tank_volume_m3.push_back(
                tank.min_volume_m3 + tank.area_m2 * (tank.level_m - tank.min_level_m));
            frame.tank_inflow_m3s.push_back(tank.boundary->Get_Q());
        }
        return frame;
    }

    void write_state(std::ofstream &node_output, std::ofstream &link_output,
                     std::ofstream &tank_output,
                     const std::vector<EpsNodeInfo> &nodes,
                     const std::vector<EpsLinkInfo> &links,
                     const std::vector<EpsTankInfo> &tanks,
                     const EpsResultFrame &frame) {
        for (std::size_t index = 0; index < nodes.size(); ++index)
            node_output << frame.time_s << ',' << csv(nodes[index].id) << ','
                        << nodes[index].elevation_m << ','
                        << frame.node_pressure_head_m[index] << ','
                        << frame.node_head_m[index] << ','
                        << frame.node_demand_m3s[index] << ','
                        << (frame.converged ? 1 : 0) << '\n';
        for (std::size_t index = 0; index < links.size(); ++index)
            link_output << frame.time_s << ',' << csv(links[index].id) << ','
                        << csv(links[index].type) << ','
                        << csv(nodes[links[index].from_node].id) << ','
                        << csv(nodes[links[index].to_node].id) << ','
                        << frame.link_flow_rate_m3s[index] << ','
                        << frame.link_velocity_ms[index] << ','
                        << frame.link_headloss_m[index] << ','
                        << static_cast<unsigned int>(frame.link_status[index]) << ','
                        << (frame.converged ? 1 : 0) << '\n';
        for (std::size_t index = 0; index < tanks.size(); ++index)
            tank_output << frame.time_s << ',' << csv(tanks[index].id) << ','
                        << frame.tank_level_m[index] << ','
                        << frame.tank_volume_m3[index] << ','
                        << frame.tank_inflow_m3s[index] << ','
                        << tanks[index].min_level_m << ',' << tanks[index].max_level_m << ','
                        << (frame.converged ? 1 : 0) << '\n';
    }

    void warn(const std::string &section, const std::string &element,
              std::size_t line, const std::string &message) {
        std::ostringstream warning;
        warning << "WARNING [EPANET][EPS][" << section << ']';
        if (!element.empty())
            warning << '[' << element << ']';
        if (line != 0)
            warning << " line " << line;
        warning << ": " << message;
        if (warning_keys_.insert(warning.str()).second)
            warnings_.push_back(warning.str());
    }
};

} // namespace

EpanetExtendedSimulation::EpanetExtendedSimulation(
    const std::string &input_filename, const std::string &output_prefix)
    : input_filename_(input_filename), output_prefix_(output_prefix) {}

void EpanetExtendedSimulation::run(Staci &system) {
    if (output_prefix_.empty())
        throw std::runtime_error("EPANET EPS requires -o <result-prefix>.");
    const std::filesystem::path prefix_path(output_prefix_);
    if (!prefix_path.parent_path().empty())
        std::filesystem::create_directories(prefix_path.parent_path());
    SimulationModel model(input_filename_);
    model.run(system, output_prefix_);
}
