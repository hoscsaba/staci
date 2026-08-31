#include "epanet_extended_simulation.h"

#include "Agelem.h"
#include "Csomopont.h"
#include "EpanetPump.h"
#include "JelleggorbesFojtas.h"
#include "Staci.h"
#include "epanet_water_age.h"
#include "epanet_steady_quality.h"
#include "eps_result_writer.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>
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

struct ChemicalSourceDefinition {
    EpanetChemicalSourceType type = EpanetChemicalSourceType::None;
    double strength_si = 0.0;
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

enum class ControlAction {
    Open,
    Closed,
    Active,
    Setting
};

enum class ControlTrigger {
    NodeAbove,
    NodeBelow,
    ElapsedTime,
    ClockTime
};

struct SimpleControl {
    std::string link_id;
    ControlAction action;
    double setting;
    ControlTrigger trigger;
    std::string node_id;
    double threshold_si;
    long long event_time_s;
    std::size_t line_number;
};

enum class RuleLogic { And, Or };
enum class RuleRelation { Equal, NotEqual, Less, LessEqual, Greater, GreaterEqual };
enum class RuleVariable {
    Demand, Head, Pressure, Level, FillTime, DrainTime,
    Flow, Status, Setting, Time, ClockTime
};
enum class RuleStatus { Open, Closed, Active };

struct RulePremise {
    RuleLogic logic;
    std::string object;
    std::string id;
    RuleVariable variable;
    RuleRelation relation;
    double value;
    double tolerance;
    RuleStatus status;
    bool status_value;
    std::size_t line_number;
};

struct RuleAction {
    std::string link_id;
    ControlAction action;
    double setting;
    std::size_t line_number;
};

struct EpanetRule {
    std::string id;
    std::vector<RulePremise> premises;
    std::vector<RuleAction> then_actions;
    std::vector<RuleAction> else_actions;
    double priority = 0.0;
};

class EpanetRuleEngine {
public:
    using NumericLookup = std::function<bool(const RulePremise &, double &)>;
    using StatusLookup = std::function<bool(const RulePremise &, RuleStatus &)>;

    void add(EpanetRule rule) { rules_.push_back(std::move(rule)); }
    bool empty() const { return rules_.empty(); }
    std::size_t size() const { return rules_.size(); }

    std::vector<RuleAction> select_actions(const NumericLookup &numeric_lookup,
                                           const StatusLookup &status_lookup,
                                           long long interval_start_s,
                                           long long time_s,
                                           long long start_clock_s) const {
        struct Selected { RuleAction action; double priority; };
        std::map<std::string, Selected> selected;
        for (const EpanetRule &rule : rules_) {
            const bool condition = evaluate_rule(
                rule, numeric_lookup, status_lookup,
                interval_start_s, time_s, start_clock_s);
            const std::vector<RuleAction> &actions = condition
                ? rule.then_actions : rule.else_actions;
            for (const RuleAction &action : actions) {
                const auto existing = selected.find(action.link_id);
                if (existing == selected.end())
                    selected.emplace(action.link_id, Selected{action, rule.priority});
                else if (rule.priority > existing->second.priority)
                    existing->second = Selected{action, rule.priority};
            }
        }
        std::vector<RuleAction> result;
        for (const auto &entry : selected)
            result.push_back(entry.second.action);
        return result;
    }

private:
    std::vector<EpanetRule> rules_;

    static bool compare(double actual, RuleRelation relation, double expected,
                        double tolerance = 0.0) {
        switch (relation) {
        case RuleRelation::Equal: return std::abs(actual - expected) <= tolerance;
        case RuleRelation::NotEqual: return std::abs(actual - expected) > tolerance;
        case RuleRelation::Less: return actual < expected;
        case RuleRelation::LessEqual: return actual <= expected;
        case RuleRelation::Greater: return actual > expected;
        case RuleRelation::GreaterEqual: return actual >= expected;
        }
        return false;
    }

    static bool time_compare(const RulePremise &premise,
                             long long interval_start_s, long long time_s,
                             long long start_clock_s) {
        long long first = interval_start_s;
        long long last = time_s;
        const long long expected = static_cast<long long>(std::llround(premise.value));
        if (premise.variable == RuleVariable::ClockTime) {
            const long long day = 24 * 3600;
            first = (first + start_clock_s) % day;
            last = (last + start_clock_s) % day;
        }
        if (premise.relation == RuleRelation::Equal ||
            premise.relation == RuleRelation::NotEqual) {
            const bool crossed = last < first
                ? expected >= first || expected <= last
                : expected >= first && expected <= last;
            return premise.relation == RuleRelation::Equal ? crossed : !crossed;
        }
        return compare(static_cast<double>(last), premise.relation,
                       static_cast<double>(expected), 0.0);
    }

    static bool evaluate_premise(const RulePremise &premise,
                                 const NumericLookup &numeric_lookup,
                                 const StatusLookup &status_lookup,
                                 long long interval_start_s, long long time_s,
                                 long long start_clock_s) {
        if (premise.variable == RuleVariable::Time ||
            premise.variable == RuleVariable::ClockTime)
            return time_compare(premise, interval_start_s, time_s, start_clock_s);
        if (premise.status_value) {
            RuleStatus actual = RuleStatus::Closed;
            if (!status_lookup(premise, actual))
                return false;
            const bool equal = actual == premise.status;
            return premise.relation == RuleRelation::NotEqual ? !equal : equal;
        }
        double actual = 0.0;
        return numeric_lookup(premise, actual) &&
               compare(actual, premise.relation, premise.value, premise.tolerance);
    }

    static bool evaluate_rule(const EpanetRule &rule,
                              const NumericLookup &numeric_lookup,
                              const StatusLookup &status_lookup,
                              long long interval_start_s, long long time_s,
                              long long start_clock_s) {
        bool result = true;
        for (const RulePremise &premise : rule.premises) {
            if (premise.logic == RuleLogic::Or) {
                if (!result)
                    result = evaluate_premise(premise, numeric_lookup, status_lookup,
                                              interval_start_s, time_s, start_clock_s);
            } else {
                if (!result)
                    return false;
                result = evaluate_premise(premise, numeric_lookup, status_lookup,
                                          interval_start_s, time_s, start_clock_s);
            }
        }
        return result;
    }
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

long long parse_duration_seconds(const std::vector<std::string> &fields,
                                 std::size_t value_index,
                                 const std::string &context) {
    if (value_index >= fields.size())
        throw std::runtime_error("Missing time value in " + context + ".");
    if (fields[value_index].find(':') != std::string::npos)
        return parse_time_seconds(fields[value_index]);
    double multiplier = 3600.0;
    if (value_index + 1 < fields.size()) {
        const std::string unit = upper(fields[value_index + 1]);
        if (unit == "SEC" || unit == "SECOND" || unit == "SECONDS") multiplier = 1.0;
        else if (unit == "MIN" || unit == "MINUTE" || unit == "MINUTES") multiplier = 60.0;
        else if (unit == "HOUR" || unit == "HOURS") multiplier = 3600.0;
        else if (unit == "DAY" || unit == "DAYS") multiplier = 86400.0;
    }
    return static_cast<long long>(std::llround(
        parse_number(fields[value_index], context) * multiplier));
}

long long parse_clock_seconds(const std::vector<std::string> &fields,
                              std::size_t value_index,
                              const std::string &context) {
    if (value_index >= fields.size())
        throw std::runtime_error("Missing clock time in " + context + ".");
    long long seconds = parse_time_seconds(fields[value_index]);
    if (value_index + 1 < fields.size()) {
        const std::string suffix = upper(fields[value_index + 1]);
        if (suffix == "AM" || suffix == "PM") {
            const long long minutes_and_seconds = seconds % 3600;
            long long hour = seconds / 3600;
            if (hour < 1 || hour > 12)
                throw std::runtime_error("Invalid 12-hour clock value in " + context + ".");
            if (hour == 12)
                hour = 0;
            if (suffix == "PM")
                hour += 12;
            seconds = hour * 3600 + minutes_and_seconds;
        }
    }
    if (seconds < 0 || seconds >= 24 * 3600)
        throw std::runtime_error("Clock time must be within one day in " + context + ".");
    return seconds;
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
          pressure_units_(""), specific_gravity_(1.0), start_clock_s_(0),
          duration_s_(0), hydraulic_step_s_(3600), pattern_step_s_(3600),
          pattern_start_s_(0), report_step_s_(3600), report_start_s_(0),
          rule_step_s_(0), quality_timestep_s_(300), quality_age_(false),
          quality_chemical_(false), chemical_units_("mg/L"),
          bulk_order_(1.0), wall_order_(1.0),
          global_bulk_per_s_(0.0), global_wall_mps_(0.0) {
        parse_file();
        parse_options();
        parse_times();
        parse_patterns();
        parse_demands();
        parse_quality();
        parse_boundaries();
        parse_controls();
        parse_rules();
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
            if (link == links.end())
                continue;
            auto *valve = dynamic_cast<JelleggorbesFojtas *>(link->second);
            if (valve != nullptr) {
                if (status.second == RuleStatus::Active)
                    valve->SetEpanetTcvStatus(EpanetTcvStatus::Active);
                else if (status.second == RuleStatus::Open)
                    valve->SetEpanetTcvStatus(EpanetTcvStatus::Open);
                else
                    valve->SetEpanetTcvStatus(EpanetTcvStatus::Closed);
            } else {
                link->second->Set_enabled(status.second != RuleStatus::Closed);
            }
        }
        for (const auto &setting : initial_settings_) {
            const auto link = links.find(setting.first);
            if (link == links.end())
                continue;
            auto *pump = dynamic_cast<EpanetPumpConfigurable *>(link->second);
            auto *valve = dynamic_cast<JelleggorbesFojtas *>(link->second);
            if (pump != nullptr) {
                pump->SetOperatingSpeed(setting.second);
                link->second->Set_enabled(setting.second > 0.0);
            } else if (valve != nullptr)
                valve->SetEpanetTcvSetting(setting.second);
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
            const bool tcv = link->GetType() == "JelleggorbesFojtas";
            std::string output_type = upper(std::string(link->GetType()));
            if (pipe)
                output_type = "PIPE";
            else if (is_pump(link))
                output_type = "PUMP";
            else if (tcv)
                output_type = "VALVE";
            result_link_objects.push_back(link);
            result_links.push_back(EpsLinkInfo{
                link->Get_nev(), output_type, from->second, to->second,
                pipe ? link->Get_dprop("length") : missing,
                pipe ? link->Get_dprop("diameter")
                     : (tcv ? std::sqrt(4.0 * link->Get_Aref() /
                                       3.14159265358979323846) : missing)});
        }

        std::vector<EpsTankInfo> result_tanks;
        for (const TankState &tank : tanks) {
            const auto node = node_indices.find(tank.id);
            if (node != node_indices.end())
                result_tanks.push_back(EpsTankInfo{
                    tank.id, node->second, tank.min_level_m, tank.max_level_m,
                    std::sqrt(4.0 * tank.area_m2 / 3.14159265358979323846)});
        }

        std::unique_ptr<EpanetWaterAgeModel> water_age;
        if (quality_age_) {
            std::vector<EpanetWaterAgeNode> age_nodes;
            age_nodes.reserve(result_nodes.size());
            for (const EpsNodeInfo &node : result_nodes)
                age_nodes.push_back(EpanetWaterAgeNode{node.type == "RESERVOIR"});
            std::vector<EpanetWaterAgeLink> age_links;
            age_links.reserve(result_links.size());
            for (const EpsLinkInfo &link : result_links) {
                double volume = 0.0;
                if (link.type == "PIPE" && std::isfinite(link.length_m) &&
                    std::isfinite(link.diameter_m) && link.length_m > 0.0 &&
                    link.diameter_m > 0.0) {
                    const double area = 3.14159265358979323846 *
                        link.diameter_m * link.diameter_m / 4.0;
                    volume = area * link.length_m;
                }
                age_links.push_back(EpanetWaterAgeLink{
                    link.from_node, link.to_node, volume});
            }
            water_age = std::make_unique<EpanetWaterAgeModel>(
                std::move(age_nodes), std::move(age_links),
                static_cast<double>(quality_timestep_s_));
            if (!tanks.empty())
                warn("QUALITY", "AGE", 0,
                     "Tank mixing is not yet represented by the EPS water-age model; tank nodes use junction mixing.");
        }

        std::unique_ptr<EpanetChemicalModel> chemical;
        if (quality_chemical_) {
            std::vector<EpanetChemicalNode> chemical_nodes;
            chemical_nodes.reserve(result_nodes.size());
            for (const EpsNodeInfo &node : result_nodes) {
                const auto initial = initial_quality_kgm3_.find(node.id);
                chemical_nodes.push_back(EpanetChemicalNode{
                    initial == initial_quality_kgm3_.end() ? 0.0 : initial->second,
                    node.type == "RESERVOIR"});
            }
            std::vector<EpanetChemicalLink> chemical_links;
            chemical_links.reserve(result_links.size());
            for (const EpsLinkInfo &link : result_links) {
                double volume = 0.0;
                double reaction = 0.0;
                if (link.type == "PIPE" && std::isfinite(link.length_m) &&
                    std::isfinite(link.diameter_m) && link.length_m > 0.0 &&
                    link.diameter_m > 0.0) {
                    const double area = 3.14159265358979323846 *
                        link.diameter_m * link.diameter_m / 4.0;
                    volume = area * link.length_m;
                    const auto bulk = pipe_bulk_per_s_.find(link.id);
                    const auto wall = pipe_wall_mps_.find(link.id);
                    reaction = bulk == pipe_bulk_per_s_.end()
                        ? global_bulk_per_s_ : bulk->second;
                    const double wall_mps = wall == pipe_wall_mps_.end()
                        ? global_wall_mps_ : wall->second;
                    reaction += 4.0 * wall_mps / link.diameter_m;
                }
                const double from = chemical_nodes[link.from_node].initial_concentration_kgm3;
                const double to = chemical_nodes[link.to_node].initial_concentration_kgm3;
                chemical_links.push_back(EpanetChemicalLink{
                    link.from_node, link.to_node, volume, reaction,
                    0.5 * (from + to)});
            }
            chemical = std::make_unique<EpanetChemicalModel>(
                std::move(chemical_nodes), std::move(chemical_links),
                static_cast<double>(quality_timestep_s_));
            if (!tanks.empty())
                warn("QUALITY", chemical_name_, 0,
                     "Tank chemical storage is not yet represented; tank nodes use instantaneous junction mixing.");
        }

        std::ofstream node_output((prefix + "-nodes.csv").c_str(), std::ios::trunc);
        std::ofstream link_output((prefix + "-links.csv").c_str(), std::ios::trunc);
        std::ofstream tank_output((prefix + "-tanks.csv").c_str(), std::ios::trunc);
        std::ofstream summary_output((prefix + "-summary.csv").c_str(), std::ios::trunc);
        if (!node_output || !link_output || !tank_output || !summary_output)
            throw std::runtime_error("Cannot create EPS CSV output files with prefix: " + prefix);

        node_output << "time_seconds,node_id,elevation_m,pressure_head_m,total_head_m,demand_m3s,converged,water_age_s,chlorine_kgm3\n";
        link_output << "time_seconds,link_id,type,node_from,node_to,flow_m3s,velocity_mps,headloss_m,status,converged,water_age_s,chlorine_kgm3\n";
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
            pattern_step_s_, report_step_s_, quality_timestep_s_, quality_age_,
            quality_chemical_, chemical_name_};
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
            apply_pump_speeds(links, time_s);
            for (TankState &tank : tanks) {
                tank.boundary->Set_dprop("water_level", tank.level_m);
                tank.boundary->Set_enabled(true);
            }
            apply_time_controls(links, time_s);

            bool converged = system.solve_system();
            const std::size_t control_iterations =
                std::max<std::size_t>(1, controls_.size() + 1);
            std::size_t control_iteration = 0;
            for (; control_iteration < control_iterations; ++control_iteration) {
                if (!apply_node_controls(links, nodes, tanks))
                    break;
                converged = system.solve_system() && converged;
            }
            if (control_iteration == control_iterations)
                warn("CONTROLS", "", 0,
                     "Node-based controls did not reach a stable state; the last state was retained.");

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
                std::vector<double> node_age;
                std::vector<double> link_age;
                std::vector<double> node_chlorine;
                std::vector<double> link_chlorine;
                std::vector<double> flows;
                flows.reserve(result_link_objects.size());
                for (Agelem *link : result_link_objects)
                    flows.push_back(link->Get_Q());
                if (water_age) {
                    node_age = water_age->node_age_s();
                    link_age = water_age->link_average_age_s(flows);
                }
                if (chemical) {
                    std::vector<double> external(result_nodes.size(), 0.0);
                    std::vector<EpanetChemicalSource> sources(result_nodes.size());
                    for (std::size_t index = 0; index < system.cspok.size(); ++index)
                        external[index] = std::max(
                            0.0, -system.cspok[index]->Get_dprop("demand") / 3600.0);
                    for (const auto &entry : chemical_sources_) {
                        const auto node = node_indices.find(entry.first);
                        if (node == node_indices.end())
                            continue;
                        sources[node->second] = EpanetChemicalSource{
                            entry.second.type,
                            entry.second.strength_si *
                                pattern_value(entry.second.pattern_id, time_s)};
                    }
                    // duration=0 updates instantaneous zero-volume links and
                    // applies boundary/source concentrations before reporting.
                    chemical->advance(0.0, flows, external, sources);
                    node_chlorine = chemical->node_concentration_kgm3();
                    link_chlorine = chemical->link_average_concentration_kgm3(flows);
                }
                const EpsResultFrame frame = collect_frame(
                    system, result_nodes, result_links, result_link_objects,
                    tanks, time_s, converged, node_age, link_age,
                    node_chlorine, link_chlorine);
                result_writer.append(frame);
                write_state(node_output, link_output, tank_output, result_nodes,
                            result_links, result_tanks, frame);
            }

            if (time_s == duration_s_)
                break;
            const long long next_grid_s =
                (time_s / simulation_step_s + 1) * simulation_step_s;
            long long next_state_s = std::min(next_grid_s, duration_s_);
            const long long control_event_s = next_time_control_event(time_s);
            if (control_event_s > time_s && control_event_s <= duration_s_)
                next_state_s = std::min(next_state_s, control_event_s);
            if (!rule_engine_.empty())
                next_state_s = scan_rules(links, nodes, tanks, time_s, next_state_s);
            const long long step_s = next_state_s - time_s;
            if (water_age) {
                std::vector<double> flows;
                flows.reserve(result_link_objects.size());
                for (Agelem *link : result_link_objects)
                    flows.push_back(link->Get_Q());
                water_age->advance(static_cast<double>(step_s), flows);
            }
            if (chemical) {
                std::vector<double> flows;
                flows.reserve(result_link_objects.size());
                for (Agelem *link : result_link_objects)
                    flows.push_back(link->Get_Q());
                std::vector<double> external(result_nodes.size(), 0.0);
                std::vector<EpanetChemicalSource> sources(result_nodes.size());
                for (std::size_t index = 0; index < system.cspok.size(); ++index)
                    external[index] = std::max(
                        0.0, -system.cspok[index]->Get_dprop("demand") / 3600.0);
                for (const auto &entry : chemical_sources_) {
                    const auto node = node_indices.find(entry.first);
                    if (node == node_indices.end())
                        continue;
                    sources[node->second] = EpanetChemicalSource{
                        entry.second.type,
                        entry.second.strength_si *
                            pattern_value(entry.second.pattern_id, time_s)};
                }
                chemical->advance(static_cast<double>(step_s), flows, external, sources);
            }
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
                       << "quality_mode," << (quality_age_ ? "AGE" :
                           (quality_chemical_ ? "CHEMICAL" : "NONE")) << "\n"
                       << "quality_timestep_seconds," << quality_timestep_s_ << "\n"
                       << "simple_controls," << controls_.size() << "\n"
                       << "rules," << rule_engine_.size() << "\n"
                       << "rule_timestep_seconds," << rule_step_s_ << "\n"
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

    void run_steady_quality(Staci &system, const std::string &prefix,
                            const std::string &requested_mode,
                            bool compute_sensitivity) {
        std::string mode = upper(trim(requested_mode));
        bool solve_age = quality_age_;
        bool solve_chemical = quality_chemical_;
        if (!mode.empty()) {
            solve_age = mode == "AGE" || mode == "BOTH";
            solve_chemical = mode == "CHEMICAL" || mode == "CHLORINE" ||
                             mode == "BOTH";
            if (!solve_age && !solve_chemical)
                throw std::runtime_error(
                    "--quality-mode must be age, chemical, chlorine or both.");
        }
        if (!solve_age && !solve_chemical)
            throw std::runtime_error(
                "No water-quality mode is selected in [OPTIONS] QUALITY or --quality-mode.");
        if (std::abs(bulk_order_ - 1.0) > 1.0e-12 ||
            std::abs(wall_order_ - 1.0) > 1.0e-12)
            throw std::runtime_error(
                "Steady chemical quality requires first-order bulk and wall reactions.");

        std::map<std::string, Csomopont *> nodes;
        for (Csomopont *node : system.cspok)
            nodes[node->Get_nev()] = node;
        std::map<std::string, Agelem *> links;
        for (Agelem *link : system.agelemek)
            links[link->Get_nev()] = link;

        std::vector<ReservoirState> reservoirs;
        for (ReservoirState reservoir : reservoir_definitions_) {
            const auto boundary = links.find("EPANET_RESERVOIR_" + reservoir.id);
            if (boundary != links.end()) {
                reservoir.boundary = boundary->second;
                reservoirs.push_back(reservoir);
            }
        }
        for (const TankState &tank : tank_definitions_) {
            const auto boundary = links.find("EPANET_TANK_" + tank.id);
            if (boundary != links.end()) {
                boundary->second->Set_dprop("water_level", tank.level_m);
                boundary->second->Set_enabled(true);
            }
        }
        if (!tank_definitions_.empty())
            warn("QUALITY", "STEADY", 0,
                 "Tank nodes are fixed quality boundaries at their imported initial values for this hydraulic snapshot; transient storage and mixing require EPS.");

        for (const auto &status : initial_status_) {
            const auto link = links.find(status.first);
            if (link == links.end())
                continue;
            auto *valve = dynamic_cast<JelleggorbesFojtas *>(link->second);
            if (valve != nullptr) {
                valve->SetEpanetTcvStatus(
                    status.second == RuleStatus::Active ? EpanetTcvStatus::Active :
                    (status.second == RuleStatus::Open ? EpanetTcvStatus::Open :
                                                        EpanetTcvStatus::Closed));
            } else {
                link->second->Set_enabled(status.second != RuleStatus::Closed);
            }
        }
        for (const auto &setting : initial_settings_) {
            const auto link = links.find(setting.first);
            if (link == links.end())
                continue;
            if (auto *pump = dynamic_cast<EpanetPumpConfigurable *>(link->second)) {
                pump->SetOperatingSpeed(setting.second);
                link->second->Set_enabled(setting.second > 0.0);
            } else if (auto *valve = dynamic_cast<JelleggorbesFojtas *>(link->second)) {
                valve->SetEpanetTcvSetting(setting.second);
            }
        }

        apply_demands(nodes, 0);
        apply_reservoir_heads(reservoirs, 0);
        apply_pump_speeds(links, 0);
        apply_time_controls(links, 0);
        const int old_debug = system.Get_debug_level();
        system.Set_debug_level(0);
        system.ini();
        const bool converged = system.solve_system();
        system.Set_debug_level(old_debug);
        if (!converged)
            throw std::runtime_error("The fixed hydraulic state did not converge.");

        std::set<std::string> reservoir_ids;
        for (const ReservoirState &reservoir : reservoir_definitions_)
            reservoir_ids.insert(reservoir.id);
        std::set<std::string> tank_ids;
        for (const TankState &tank : tank_definitions_)
            tank_ids.insert(tank.id);
        std::vector<EpanetSteadyQualityNode> quality_nodes;
        std::map<std::string, std::size_t> node_index;
        quality_nodes.reserve(system.cspok.size());
        for (std::size_t index = 0; index < system.cspok.size(); ++index) {
            Csomopont *node = system.cspok[index];
            const std::string id = node->Get_nev();
            node_index[id] = index;
            const auto initial = initial_quality_kgm3_.find(id);
            const bool fixed_boundary = reservoir_ids.count(id) != 0 ||
                                        tank_ids.count(id) != 0;
            EpanetSteadyQualityNode quality_node;
            quality_node.id = id;
            quality_node.fixed_age = fixed_boundary;
            quality_node.fixed_chemical = fixed_boundary;
            quality_node.fixed_concentration_kgm3 =
                initial == initial_quality_kgm3_.end() ? 0.0 : initial->second;
            quality_node.external_inflow_m3s =
                std::max(0.0, -node->Get_dprop("demand") / 3600.0);
            quality_nodes.push_back(quality_node);
        }

        std::vector<Agelem *> quality_link_objects;
        std::vector<EpanetSteadyQualityLink> quality_links;
        std::vector<double> wall_coefficients;
        for (Agelem *link : system.agelemek) {
            if (link->Get_nev().rfind("EPANET_TANK_", 0) == 0 ||
                link->Get_nev().rfind("EPANET_RESERVOIR_", 0) == 0)
                continue;
            const auto from = node_index.find(link->Get_Cspe_Nev());
            const auto to = node_index.find(link->Get_Cspv_Nev());
            if (from == node_index.end() || to == node_index.end())
                continue;
            double volume = 0.0;
            double reaction = 0.0;
            double wall_mps = 0.0;
            if (link->GetType() == "Cso") {
                const double length = link->Get_dprop("length");
                const double diameter = link->Get_dprop("diameter");
                if (length > 0.0 && diameter > 0.0) {
                    volume = 3.14159265358979323846 * diameter * diameter *
                             length / 4.0;
                    const auto bulk = pipe_bulk_per_s_.find(link->Get_nev());
                    const auto wall = pipe_wall_mps_.find(link->Get_nev());
                    reaction = bulk == pipe_bulk_per_s_.end()
                        ? global_bulk_per_s_ : bulk->second;
                    wall_mps = wall == pipe_wall_mps_.end()
                        ? global_wall_mps_ : wall->second;
                    reaction += 4.0 * wall_mps / diameter;
                }
            }
            quality_link_objects.push_back(link);
            quality_links.push_back(EpanetSteadyQualityLink{
                link->Get_nev(), from->second, to->second, link->Get_Q(),
                volume, reaction});
            wall_coefficients.push_back(wall_mps);
        }

        std::vector<double> total_inflow(quality_nodes.size(), 0.0);
        for (std::size_t index = 0; index < quality_nodes.size(); ++index)
            total_inflow[index] = quality_nodes[index].external_inflow_m3s;
        for (const EpanetSteadyQualityLink &link : quality_links) {
            if (link.flow_m3s >= 0.0)
                total_inflow[link.to_node] += link.flow_m3s;
            else
                total_inflow[link.from_node] -= link.flow_m3s;
        }
        for (const auto &entry : chemical_sources_) {
            const auto found = node_index.find(entry.first);
            if (found == node_index.end())
                continue;
            EpanetSteadyQualityNode &node = quality_nodes[found->second];
            const double strength = entry.second.strength_si *
                                    pattern_value(entry.second.pattern_id, 0);
            if (entry.second.type == EpanetChemicalSourceType::Mass)
                node.mass_source_kgs += strength;
            else if (entry.second.type == EpanetChemicalSourceType::Concentration) {
                if (node.fixed_chemical)
                    node.fixed_concentration_kgm3 = strength;
                else
                    node.external_concentration_kgm3 = strength;
            } else if (entry.second.type == EpanetChemicalSourceType::FlowPaced) {
                node.mass_source_kgs += total_inflow[found->second] * strength;
            } else if (entry.second.type == EpanetChemicalSourceType::Setpoint) {
                node.fixed_chemical = true;
                node.fixed_concentration_kgm3 = strength;
                warn("SOURCES", entry.first, 0,
                     "SETPOINT was represented as a fixed steady concentration; the conditional no-removal behavior is nonlinear and is not applied.");
            }
        }

        EpanetSteadyQualitySensitivityInput sensitivity;
        const EpanetSteadyQualitySensitivityInput *sensitivity_pointer = nullptr;
        if (compute_sensitivity) {
            if (system.element_ID.empty() || system.property_ID.empty())
                throw std::runtime_error(
                    "Steady-quality sensitivity requires -e <element> and -p <property>.");
            if (system.property_ID != "diameter" &&
                system.property_ID != "friction_coeff" &&
                system.property_ID != "demand")
                throw std::runtime_error(
                    "Steady-quality sensitivity supports diameter, friction_coeff and demand.");
            system.Compute_dxdmu();
            const std::vector<double> &hydraulic = system.Get_dxdmu();
            sensitivity.parameter_element = system.element_ID;
            sensitivity.parameter_property = system.property_ID;
            sensitivity.link_flow_derivative_m3s.assign(quality_links.size(), 0.0);
            sensitivity.link_volume_derivative_m3.assign(quality_links.size(), 0.0);
            sensitivity.link_reaction_derivative_per_s.assign(quality_links.size(), 0.0);
            sensitivity.external_inflow_derivative_m3s.assign(quality_nodes.size(), 0.0);
            const double parameter_scale = system.property_ID == "demand" ? 3600.0 : 1.0;
            for (std::size_t qindex = 0; qindex < quality_link_objects.size(); ++qindex) {
                Agelem *object = quality_link_objects[qindex];
                const auto original = std::find(system.agelemek.begin(),
                                                system.agelemek.end(), object);
                if (original != system.agelemek.end()) {
                    const std::size_t hydraulic_index = static_cast<std::size_t>(
                        std::distance(system.agelemek.begin(), original));
                    sensitivity.link_flow_derivative_m3s[qindex] =
                        hydraulic[hydraulic_index] / object->Get_ro() * parameter_scale;
                }
                if (system.property_ID == "diameter" &&
                    object->Get_nev() == system.element_ID &&
                    object->GetType() == "Cso") {
                    const double diameter = object->Get_dprop("diameter");
                    sensitivity.link_volume_derivative_m3[qindex] =
                        2.0 * quality_links[qindex].volume_m3 / diameter;
                    sensitivity.link_reaction_derivative_per_s[qindex] =
                        -4.0 * wall_coefficients[qindex] / (diameter * diameter);
                }
            }
            if (system.property_ID == "demand") {
                const auto selected = node_index.find(system.element_ID);
                if (selected != node_index.end() &&
                    system.cspok[selected->second]->Get_dprop("demand") < 0.0)
                    sensitivity.external_inflow_derivative_m3s[selected->second] = -1.0;
            }
            sensitivity_pointer = &sensitivity;
        }

        EpanetSteadyQualityModel model(std::move(quality_nodes), quality_links);
        const EpanetSteadyQualityResult result =
            model.solve(solve_age, solve_chemical, sensitivity_pointer);
        const double missing = std::numeric_limits<double>::quiet_NaN();
        std::ofstream node_output((prefix + "-steady-nodes.csv").c_str(),
                                  std::ios::trunc);
        std::ofstream link_output((prefix + "-steady-links.csv").c_str(),
                                  std::ios::trunc);
        std::ofstream summary_output((prefix + "-steady-summary.csv").c_str(),
                                     std::ios::trunc);
        if (!node_output || !link_output || !summary_output)
            throw std::runtime_error("Cannot create steady-quality SI CSV outputs.");
        node_output << "node_id,water_age_s,chemical_concentration_kgm3\n"
                    << std::setprecision(15);
        for (std::size_t index = 0; index < system.cspok.size(); ++index)
            node_output << csv(system.cspok[index]->Get_nev()) << ','
                        << (solve_age ? result.node_age_s[index] : missing) << ','
                        << (solve_chemical ? result.node_concentration_kgm3[index] : missing)
                        << '\n';
        link_output << "link_id,node_from,node_to,flow_m3s,travel_time_s,water_age_s,chemical_concentration_kgm3,transfer_factor\n"
                    << std::setprecision(15);
        for (std::size_t index = 0; index < quality_links.size(); ++index)
            link_output << csv(quality_links[index].id) << ','
                        << csv(system.cspok[quality_links[index].from_node]->Get_nev()) << ','
                        << csv(system.cspok[quality_links[index].to_node]->Get_nev()) << ','
                        << quality_links[index].flow_m3s << ','
                        << result.link_travel_time_s[index] << ','
                        << (solve_age ? result.link_average_age_s[index] : missing) << ','
                        << (solve_chemical ? result.link_average_concentration_kgm3[index] : missing) << ','
                        << result.link_transfer_factor[index] << '\n';
        summary_output << "property,value\ninput_file," << csv(filename_)
                       << "\nquality_mode," << (solve_age && solve_chemical ? "BOTH" :
                           (solve_age ? "AGE" : "CHEMICAL"))
                       << "\nhydraulic_converged,1\nnodes," << system.cspok.size()
                       << "\nlinks," << quality_links.size()
                       << "\nsensitivity," << (compute_sensitivity ? 1 : 0) << '\n';

        if (compute_sensitivity) {
            std::ofstream output((prefix + "-steady-sensitivity.csv").c_str(),
                                 std::ios::trunc);
            if (!output)
                throw std::runtime_error("Cannot create steady-quality sensitivity CSV.");
            output << "parameter_element,parameter_property,result_quantity,result_id,derivative_per_si_parameter\n"
                   << std::setprecision(15);
            for (std::size_t index = 0; index < system.cspok.size(); ++index) {
                if (solve_age)
                    output << csv(system.element_ID) << ',' << csv(system.property_ID)
                           << ",node_water_age_s," << csv(system.cspok[index]->Get_nev())
                           << ',' << result.node_age_sensitivity[index] << '\n';
                if (solve_chemical)
                    output << csv(system.element_ID) << ',' << csv(system.property_ID)
                           << ",node_chemical_concentration_kgm3,"
                           << csv(system.cspok[index]->Get_nev()) << ','
                           << result.node_concentration_sensitivity[index] << '\n';
            }
            for (std::size_t index = 0; index < quality_links.size(); ++index) {
                output << csv(system.element_ID) << ',' << csv(system.property_ID)
                       << ",link_flow_m3s," << csv(quality_links[index].id)
                       << ',' << sensitivity.link_flow_derivative_m3s[index] << '\n';
                if (solve_age)
                    output << csv(system.element_ID) << ',' << csv(system.property_ID)
                           << ",link_water_age_s," << csv(quality_links[index].id)
                           << ',' << result.link_average_age_sensitivity[index] << '\n';
                if (solve_chemical)
                    output << csv(system.element_ID) << ',' << csv(system.property_ID)
                           << ",link_chemical_concentration_kgm3,"
                           << csv(quality_links[index].id) << ','
                           << result.link_average_concentration_sensitivity[index] << '\n';
            }
        }
        for (const std::string &warning : warnings_)
            std::cerr << warning << '\n';
        std::cout << "\nSteady water-quality solution complete. SI results: "
                  << prefix << "-steady-nodes.csv and " << prefix
                  << "-steady-links.csv.\n";
    }

private:
    std::string filename_;
    std::map<std::string, std::vector<Record> > sections_;
    std::map<std::string, std::vector<double> > patterns_;
    std::map<std::string, std::vector<DemandComponent> > demands_;
    std::map<std::string, double> initial_quality_kgm3_;
    std::map<std::string, ChemicalSourceDefinition> chemical_sources_;
    std::map<std::string, double> pipe_bulk_per_s_;
    std::map<std::string, double> pipe_wall_mps_;
    std::vector<TankState> tank_definitions_;
    std::vector<ReservoirState> reservoir_definitions_;
    std::vector<SimpleControl> controls_;
    EpanetRuleEngine rule_engine_;
    std::map<std::string, RuleStatus> initial_status_;
    std::map<std::string, double> initial_settings_;
    std::map<std::string, double> controlled_pump_settings_;
    std::vector<std::string> warnings_;
    std::set<std::string> warning_keys_;
    std::set<std::string> missing_patterns_;
    std::string flow_units_;
    std::string default_pattern_;
    double demand_multiplier_;
    std::string pressure_units_;
    double specific_gravity_;
    long long start_clock_s_;
    long long duration_s_;
    long long hydraulic_step_s_;
    long long pattern_step_s_;
    long long pattern_start_s_;
    long long report_step_s_;
    long long report_start_s_;
    long long rule_step_s_;
    long long quality_timestep_s_;
    bool quality_age_;
    bool quality_chemical_;
    std::string chemical_name_;
    std::string chemical_units_;
    double bulk_order_;
    double wall_order_;
    double global_bulk_per_s_;
    double global_wall_mps_;

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
            else if (key == "PRESSURE")
                pressure_units_ = upper(record.fields[1]);
            else if (key == "QUALITY") {
                const std::string mode = upper(record.fields[1]);
                quality_age_ = mode == "AGE";
                quality_chemical_ = mode == "CHEMICAL";
                if (quality_chemical_) {
                    chemical_name_ = record.fields.size() > 2
                        ? record.fields[2] : "CHEMICAL";
                    chemical_units_ = record.fields.size() > 3
                        ? record.fields[3] : "mg/L";
                }
            }
            else if (record.fields.size() > 2 && key == "SPECIFIC" &&
                     upper(record.fields[1]) == "GRAVITY")
                specific_gravity_ = parse_number(record.fields[2], "SPECIFIC GRAVITY");
        }
        if (pressure_units_.empty())
            pressure_units_ = us_units(flow_units_) ? "PSI" : "METERS";
        if (specific_gravity_ <= 0.0)
            throw std::runtime_error("EPANET specific gravity must be positive.");
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
            else if (record.fields.size() > 2 && first == "QUALITY" &&
                     upper(record.fields[1]) == "TIMESTEP")
                quality_timestep_s_ = parse_time_seconds(record.fields[2]);
            else if (record.fields.size() > 2 && first == "RULE" &&
                     upper(record.fields[1]) == "TIMESTEP")
                rule_step_s_ = parse_duration_seconds(record.fields, 2, "[TIMES] RULE TIMESTEP");
            else if (record.fields.size() > 2 && first == "START" &&
                     upper(record.fields[1]) == "CLOCKTIME")
                start_clock_s_ = parse_clock_seconds(record.fields, 2, "[TIMES] START CLOCKTIME");
        }
        if (rule_step_s_ <= 0)
            rule_step_s_ = std::max(1LL, hydraulic_step_s_ / 10);
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
        std::set<std::string> overridden_junctions;
        for (const Record &record : records("DEMANDS"))
            if (record.fields.size() >= 2)
                overridden_junctions.insert(record.fields[0]);
        for (const Record &record : records("JUNCTIONS")) {
            if (record.fields.size() < 2)
                continue;
            const double base = record.fields.size() > 2
                ? parse_number(record.fields[2], "[JUNCTIONS]") : 0.0;
            const std::string pattern = record.fields.size() > 3
                ? record.fields[3] : default_pattern_;
            if (overridden_junctions.count(record.fields[0]) == 0)
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

    double concentration_to_si(double value) {
        const std::string unit = upper(chemical_units_);
        if (unit == "MG/L" || unit == "MG/LITER" || unit == "MG/LITRE")
            return value * 1.0e-3;
        if (unit == "UG/L" || unit == "UG/LITER" || unit == "UG/LITRE")
            return value * 1.0e-6;
        if (unit == "KG/M3" || unit == "KG/M^3")
            return value;
        warn("OPTIONS", "QUALITY", 0,
             "Unknown chemical concentration unit '" + chemical_units_ +
             "'; values were interpreted as mg/L and converted to SI kg/m3.");
        return value * 1.0e-3;
    }

    double mass_source_to_si(double value) {
        const std::string unit = upper(chemical_units_);
        if (unit == "UG/L" || unit == "UG/LITER" || unit == "UG/LITRE")
            return value * 1.0e-9 / 60.0;
        if (unit == "KG/M3" || unit == "KG/M^3")
            return value / 60.0;
        // EPANET MASS sources use mass/minute; mg/min is paired with the
        // conventional mg/L chemical concentration declaration.
        return value * 1.0e-6 / 60.0;
    }

    void parse_quality() {
        for (const Record &record : records("QUALITY")) {
            if (record.fields.size() >= 2)
                initial_quality_kgm3_[record.fields[0]] = concentration_to_si(
                    parse_number(record.fields[1], "[QUALITY]"));
        }
        for (const Record &record : records("SOURCES")) {
            if (record.fields.size() < 3) {
                warn("SOURCES", record.fields.empty() ? "" : record.fields[0],
                     record.line_number, "Incomplete chemical source was ignored.");
                continue;
            }
            const std::string type = upper(record.fields[1]);
            ChemicalSourceDefinition source;
            if (type == "CONCEN")
                source.type = EpanetChemicalSourceType::Concentration;
            else if (type == "MASS")
                source.type = EpanetChemicalSourceType::Mass;
            else if (type == "SETPOINT")
                source.type = EpanetChemicalSourceType::Setpoint;
            else if (type == "FLOWPACED")
                source.type = EpanetChemicalSourceType::FlowPaced;
            else {
                warn("SOURCES", record.fields[0], record.line_number,
                     "Unknown source type '" + record.fields[1] + "' was ignored.");
                continue;
            }
            const double raw = parse_number(record.fields[2], "[SOURCES]");
            // EPANET MASS source strength is mass/minute; the other source
            // strengths use the configured concentration unit.
            source.strength_si = source.type == EpanetChemicalSourceType::Mass
                ? mass_source_to_si(raw) : concentration_to_si(raw);
            if (record.fields.size() > 3)
                source.pattern_id = record.fields[3];
            chemical_sources_[record.fields[0]] = source;
        }
        for (const Record &record : records("REACTIONS")) {
            if (record.fields.size() < 3)
                continue;
            const std::string first = upper(record.fields[0]);
            const std::string second = upper(record.fields[1]);
            const double value = parse_number(record.fields[2], "[REACTIONS]");
            if (first == "ORDER" && second == "BULK")
                bulk_order_ = value;
            else if (first == "ORDER" && second == "WALL")
                wall_order_ = value;
            else if (first == "GLOBAL" && second == "BULK")
                global_bulk_per_s_ = value / 86400.0;
            else if (first == "GLOBAL" && second == "WALL")
                global_wall_mps_ = length_to_metres(value, us_units(flow_units_)) / 86400.0;
            else if (first == "BULK")
                pipe_bulk_per_s_[record.fields[1]] = value / 86400.0;
            else if (first == "WALL")
                pipe_wall_mps_[record.fields[1]] =
                    length_to_metres(value, us_units(flow_units_)) / 86400.0;
        }
        if (std::abs(bulk_order_ - 1.0) > 1.0e-12 ||
            std::abs(wall_order_ - 1.0) > 1.0e-12)
            warn("REACTIONS", "ORDER", 0,
                 "STACI chemical EPS currently applies first-order bulk and wall reactions.");
        if (!records("MIXING").empty())
            warn("MIXING", "", records("MIXING").front().line_number,
                 "Tank chemical mixing is not yet supported; junction-style instantaneous mixing is used.");
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
        std::set<std::string> tank_ids;
        for (const TankState &tank : tank_definitions_)
            tank_ids.insert(tank.id);
        for (const Record &record : records("CONTROLS")) {
            if (record.fields.size() < 5 || upper(record.fields[0]) != "LINK") {
                warn("CONTROLS", record.fields.size() > 1 ? record.fields[1] : "",
                     record.line_number, "Invalid simple control syntax; control was ignored.");
                continue;
            }

            ControlAction action = ControlAction::Setting;
            double setting = 0.0;
            const std::string action_text = upper(record.fields[2]);
            if (action_text == "OPEN")
                action = ControlAction::Open;
            else if (action_text == "CLOSED")
                action = ControlAction::Closed;
            else if (action_text == "ACTIVE")
                action = ControlAction::Active;
            else {
                setting = parse_number(record.fields[2], "[CONTROLS] link setting");
                if (setting < 0.0) {
                    warn("CONTROLS", record.fields[1], record.line_number,
                         "Negative link setting was ignored.");
                    continue;
                }
            }

            const std::string trigger_keyword = upper(record.fields[3]);
            if (trigger_keyword == "IF" && record.fields.size() >= 8 &&
                upper(record.fields[4]) == "NODE" &&
                (upper(record.fields[6]) == "ABOVE" || upper(record.fields[6]) == "BELOW")) {
                const std::string node_id = record.fields[5];
                double threshold = parse_number(record.fields[7], "[CONTROLS] node threshold");
                if (tank_ids.count(node_id) != 0) {
                    threshold = length_to_metres(threshold, is_us);
                } else if (pressure_units_ == "KPA") {
                    threshold *= 1000.0 / (1000.0 * specific_gravity_ * 9.81);
                } else if (pressure_units_ == "PSI") {
                    threshold *= 6894.757293168 / (1000.0 * specific_gravity_ * 9.81);
                } else if (pressure_units_ == "FEET") {
                    threshold *= 0.3048;
                } else if (pressure_units_ != "METERS") {
                    warn("CONTROLS", record.fields[1], record.line_number,
                         "Unknown pressure unit; node threshold was interpreted as metres of head.");
                }
                controls_.push_back(SimpleControl{
                    record.fields[1], action, setting,
                    upper(record.fields[6]) == "ABOVE"
                        ? ControlTrigger::NodeAbove : ControlTrigger::NodeBelow,
                    node_id, threshold, 0, record.line_number});
            } else if (trigger_keyword == "AT" && record.fields.size() >= 6 &&
                       upper(record.fields[4]) == "TIME") {
                const long long event_s = parse_time_seconds(record.fields[5]);
                if (event_s < 0) {
                    warn("CONTROLS", record.fields[1], record.line_number,
                         "Negative elapsed control time was ignored.");
                    continue;
                }
                controls_.push_back(SimpleControl{
                    record.fields[1], action, setting, ControlTrigger::ElapsedTime,
                    "", 0.0, event_s, record.line_number});
            } else if (trigger_keyword == "AT" && record.fields.size() >= 6 &&
                       upper(record.fields[4]) == "CLOCKTIME") {
                controls_.push_back(SimpleControl{
                    record.fields[1], action, setting, ControlTrigger::ClockTime,
                    "", 0.0,
                    parse_clock_seconds(record.fields, 5, "[CONTROLS] CLOCKTIME"),
                    record.line_number});
            } else {
                warn("CONTROLS", record.fields[1], record.line_number,
                     "Invalid simple control trigger; control was ignored.");
            }
        }
    }

    double pressure_to_head_m(double value) const {
        if (pressure_units_ == "KPA")
            return value * 1000.0 / (1000.0 * specific_gravity_ * 9.81);
        if (pressure_units_ == "PSI")
            return value * 6894.757293168 / (1000.0 * specific_gravity_ * 9.81);
        if (pressure_units_ == "FEET")
            return value * 0.3048;
        return value;
    }

    static RuleRelation parse_rule_relation(const std::string &text) {
        const std::string value = upper(text);
        if (value == "=" || value == "IS") return RuleRelation::Equal;
        if (value == "<>" || value == "NOT") return RuleRelation::NotEqual;
        if (value == "<" || value == "BELOW") return RuleRelation::Less;
        if (value == "<=") return RuleRelation::LessEqual;
        if (value == ">" || value == "ABOVE") return RuleRelation::Greater;
        if (value == ">=") return RuleRelation::GreaterEqual;
        throw std::runtime_error("Unknown rule relation '" + text + "'.");
    }

    static RuleVariable parse_rule_variable(const std::string &text) {
        const std::string value = upper(text);
        if (value == "DEMAND") return RuleVariable::Demand;
        if (value == "HEAD" || value == "GRADE") return RuleVariable::Head;
        if (value == "PRESSURE") return RuleVariable::Pressure;
        if (value == "LEVEL") return RuleVariable::Level;
        if (value == "FILLTIME") return RuleVariable::FillTime;
        if (value == "DRAINTIME") return RuleVariable::DrainTime;
        if (value == "FLOW") return RuleVariable::Flow;
        if (value == "STATUS") return RuleVariable::Status;
        if (value == "SETTING") return RuleVariable::Setting;
        if (value == "TIME") return RuleVariable::Time;
        if (value == "CLOCKTIME") return RuleVariable::ClockTime;
        throw std::runtime_error("Unknown rule variable '" + text + "'.");
    }

    RulePremise parse_rule_premise(const Record &record, RuleLogic logic) const {
        if (record.fields.size() < 5)
            throw std::runtime_error("Incomplete rule premise.");
        const std::string object = upper(record.fields[1]);
        const bool system = object == "SYSTEM";
        const std::size_t variable_index = system ? 2 : 3;
        const std::size_t relation_index = system ? 3 : 4;
        const std::size_t value_index = system ? 4 : 5;
        if (value_index >= record.fields.size())
            throw std::runtime_error("Incomplete rule premise.");
        static const std::set<std::string> valid_objects = {
            "SYSTEM", "NODE", "JUNCTION", "RESERVOIR", "TANK",
            "LINK", "PIPE", "PUMP", "VALVE"};
        if (valid_objects.count(object) == 0)
            throw std::runtime_error("Unknown rule object '" + object + "'.");

        const RuleVariable variable = parse_rule_variable(record.fields[variable_index]);
        const RuleRelation relation = parse_rule_relation(record.fields[relation_index]);
        const bool node_object = object == "NODE" || object == "JUNCTION" ||
            object == "RESERVOIR" || object == "TANK";
        const bool link_object = object == "LINK" || object == "PIPE" ||
            object == "PUMP" || object == "VALVE";
        if (system && variable != RuleVariable::Demand && variable != RuleVariable::Time &&
            variable != RuleVariable::ClockTime)
            throw std::runtime_error("Invalid SYSTEM rule variable.");
        if (node_object && variable != RuleVariable::Demand && variable != RuleVariable::Head &&
            variable != RuleVariable::Pressure && variable != RuleVariable::Level &&
            variable != RuleVariable::FillTime && variable != RuleVariable::DrainTime)
            throw std::runtime_error("Invalid node rule variable.");
        if (link_object && variable != RuleVariable::Flow && variable != RuleVariable::Status &&
            variable != RuleVariable::Setting)
            throw std::runtime_error("Invalid link rule variable.");

        RulePremise premise{logic, object, system ? "" : record.fields[2], variable,
                            relation, 0.0, 1.0e-3, RuleStatus::Closed, false,
                            record.line_number};
        const std::string raw_value = upper(record.fields[value_index]);
        if (variable == RuleVariable::Status) {
            if (relation != RuleRelation::Equal && relation != RuleRelation::NotEqual)
                throw std::runtime_error("STATUS premises require IS, NOT, =, or <>.");
            premise.status_value = true;
            if (raw_value == "OPEN") premise.status = RuleStatus::Open;
            else if (raw_value == "CLOSED") premise.status = RuleStatus::Closed;
            else if (raw_value == "ACTIVE") premise.status = RuleStatus::Active;
            else throw std::runtime_error("Unknown rule status '" + raw_value + "'.");
            return premise;
        }
        if (variable == RuleVariable::ClockTime) {
            premise.value = static_cast<double>(
                parse_clock_seconds(record.fields, value_index, "[RULES] CLOCKTIME"));
            return premise;
        }
        if (variable == RuleVariable::Time) {
            premise.value = static_cast<double>(
                parse_duration_seconds(record.fields, value_index, "[RULES] TIME"));
            return premise;
        }

        double value = parse_number(record.fields[value_index], "[RULES]");
        if (variable == RuleVariable::Demand) {
            value = flow_to_m3_per_hour(value, flow_units_);
            premise.tolerance = flow_to_m3_per_hour(1.0e-3, flow_units_);
        } else if (variable == RuleVariable::Flow) {
            value = flow_to_m3_per_hour(value, flow_units_) / 3600.0;
            premise.tolerance = flow_to_m3_per_hour(1.0e-3, flow_units_) / 3600.0;
        } else if (variable == RuleVariable::Head || variable == RuleVariable::Level) {
            value = length_to_metres(value, us_units(flow_units_));
            premise.tolerance = length_to_metres(1.0e-3, us_units(flow_units_));
        } else if (variable == RuleVariable::Pressure) {
            value = pressure_to_head_m(value);
            premise.tolerance = pressure_to_head_m(1.0e-3);
        } else if (variable == RuleVariable::FillTime || variable == RuleVariable::DrainTime) {
            value *= 3600.0;
            premise.tolerance = 3.6;
        }
        premise.value = value;
        return premise;
    }

    RuleAction parse_rule_action(const Record &record) const {
        if (record.fields.size() != 6 || upper(record.fields[4]) != "IS")
            throw std::runtime_error("Invalid rule action syntax.");
        const std::string object = upper(record.fields[1]);
        if (object != "LINK" && object != "PIPE" && object != "PUMP" && object != "VALVE")
            throw std::runtime_error("Rule action must target a link.");
        const std::string attribute = upper(record.fields[3]);
        if (attribute != "STATUS" && attribute != "SETTING")
            throw std::runtime_error("Rule action requires STATUS or SETTING.");
        const std::string value = upper(record.fields[5]);
        if (value == "OPEN")
            return RuleAction{record.fields[2], ControlAction::Open, 0.0, record.line_number};
        if (value == "CLOSED")
            return RuleAction{record.fields[2], ControlAction::Closed, 0.0, record.line_number};
        if (value == "ACTIVE")
            return RuleAction{record.fields[2], ControlAction::Active, 0.0, record.line_number};
        const double setting = parse_number(record.fields[5], "[RULES] action setting");
        if (setting < 0.0)
            throw std::runtime_error("Negative rule settings are not supported.");
        return RuleAction{record.fields[2], ControlAction::Setting, setting,
                          record.line_number};
    }

    void parse_rules() {
        enum class Phase { None, Premises, ThenActions, ElseActions };
        EpanetRule current;
        Phase phase = Phase::None;
        std::size_t rule_line = 0;
        auto finish_rule = [&]() {
            if (current.id.empty())
                return;
            if (current.premises.empty() || current.then_actions.empty())
                warn("RULES", current.id, rule_line,
                     "Rule requires at least IF and THEN clauses; it was ignored.");
            else
                rule_engine_.add(std::move(current));
            current = EpanetRule{};
            phase = Phase::None;
        };

        for (const Record &record : records("RULES")) {
            if (record.fields.empty())
                continue;
            const std::string keyword = upper(record.fields[0]);
            try {
                if (keyword == "RULE") {
                    finish_rule();
                    if (record.fields.size() < 2)
                        throw std::runtime_error("RULE requires an ID.");
                    current.id = record.fields[1];
                    rule_line = record.line_number;
                } else if (current.id.empty()) {
                    throw std::runtime_error("Clause appears before RULE.");
                } else if (keyword == "IF") {
                    phase = Phase::Premises;
                    current.premises.push_back(parse_rule_premise(record, RuleLogic::And));
                } else if (keyword == "OR") {
                    if (phase != Phase::Premises)
                        throw std::runtime_error("OR is only valid between premises.");
                    current.premises.push_back(parse_rule_premise(record, RuleLogic::Or));
                } else if (keyword == "THEN") {
                    phase = Phase::ThenActions;
                    current.then_actions.push_back(parse_rule_action(record));
                } else if (keyword == "ELSE") {
                    phase = Phase::ElseActions;
                    current.else_actions.push_back(parse_rule_action(record));
                } else if (keyword == "AND") {
                    if (phase == Phase::Premises)
                        current.premises.push_back(parse_rule_premise(record, RuleLogic::And));
                    else if (phase == Phase::ThenActions)
                        current.then_actions.push_back(parse_rule_action(record));
                    else if (phase == Phase::ElseActions)
                        current.else_actions.push_back(parse_rule_action(record));
                    else
                        throw std::runtime_error("Misplaced AND clause.");
                } else if (keyword == "PRIORITY") {
                    if (record.fields.size() < 2)
                        throw std::runtime_error("PRIORITY requires a value.");
                    current.priority = parse_number(record.fields[1], "[RULES] PRIORITY");
                } else {
                    throw std::runtime_error("Unknown rule clause '" + keyword + "'.");
                }
            } catch (const std::exception &error) {
                warn("RULES", current.id, record.line_number, error.what());
            }
        }
        finish_rule();
    }

    void parse_status() {
        for (const Record &record : records("STATUS")) {
            if (record.fields.size() < 2)
                continue;
            const std::string status = upper(record.fields[1]);
            if (status == "OPEN")
                initial_status_[record.fields[0]] = RuleStatus::Open;
            else if (status == "CLOSED")
                initial_status_[record.fields[0]] = RuleStatus::Closed;
            else if (status == "ACTIVE")
                initial_status_[record.fields[0]] = RuleStatus::Active;
            else {
                const double setting = parse_number(record.fields[1], "[STATUS]");
                if (setting < 0.0)
                    warn("STATUS", record.fields[0], record.line_number,
                         "Negative link settings are not supported.");
                else
                    initial_settings_[record.fields[0]] = setting;
            }
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

    void apply_pump_speeds(const std::map<std::string, Agelem *> &links,
                           long long time_s) {
        for (const auto &link : links) {
            auto *pump = dynamic_cast<EpanetPumpConfigurable *>(link.second);
            if (pump == nullptr)
                continue;
            const EpanetPumpMetadata &metadata = pump->GetEpanetPumpMetadata();
            double speed = !metadata.speed_pattern_id.empty()
                ? pattern_value(metadata.speed_pattern_id, time_s)
                : (metadata.initial_setting_specified
                       ? metadata.initial_setting : metadata.base_speed);
            const auto controlled = controlled_pump_settings_.find(link.first);
            if (controlled != controlled_pump_settings_.end())
                speed = controlled->second;
            pump->SetOperatingSpeed(speed);
        }
    }

    bool apply_control_action(const SimpleControl &control,
                              const std::map<std::string, Agelem *> &links,
                              const std::string &section = "CONTROLS") {
        const auto found = links.find(control.link_id);
        if (found == links.end()) {
            warn(section, control.link_id, control.line_number,
                 "Controlled link was not found in the imported STACI system.");
            return false;
        }
        Agelem *link = found->second;
        const bool old_enabled = link->Is_enabled();
        auto *pump = dynamic_cast<EpanetPumpConfigurable *>(link);
        auto *valve = dynamic_cast<JelleggorbesFojtas *>(link);
        const double old_setting = pump != nullptr ? pump->GetOperatingSpeed()
            : (valve != nullptr ? valve->GetEpanetTcvSetting() : 0.0);
        const EpanetTcvStatus old_valve_status = valve == nullptr
            ? EpanetTcvStatus::Closed : valve->GetEpanetTcvStatus();

        if (control.action == ControlAction::Open) {
            if (valve != nullptr)
                valve->SetEpanetTcvStatus(EpanetTcvStatus::Open);
            else
                link->Set_enabled(true);
        } else if (control.action == ControlAction::Closed) {
            if (valve != nullptr)
                valve->SetEpanetTcvStatus(EpanetTcvStatus::Closed);
            else
                link->Set_enabled(false);
        } else if (control.action == ControlAction::Active) {
            if (valve == nullptr) {
                warn(section, control.link_id, control.line_number,
                     "ACTIVE status requires an imported EPANET TCV; the control was ignored.");
                return false;
            }
            valve->SetEpanetTcvStatus(EpanetTcvStatus::Active);
        } else if (pump != nullptr) {
            pump->SetOperatingSpeed(control.setting);
            controlled_pump_settings_[control.link_id] = control.setting;
            link->Set_enabled(control.setting > 0.0);
        } else if (valve != nullptr) {
            valve->SetEpanetTcvSetting(control.setting);
        } else if (link->GetType() == "Cso") {
            link->Set_enabled(control.setting > 0.0);
        } else {
            warn(section, control.link_id, control.line_number,
                 "Numeric settings currently require a pipe, EPANET pump, or TCV; the control was ignored.");
            return false;
        }

        return link->Is_enabled() != old_enabled ||
               (pump != nullptr &&
                std::abs(pump->GetOperatingSpeed() - old_setting) > 1.0e-12) ||
               (valve != nullptr &&
                (std::abs(valve->GetEpanetTcvSetting() - old_setting) > 1.0e-12 ||
                 valve->GetEpanetTcvStatus() != old_valve_status));
    }

    bool time_control_is_active(const SimpleControl &control, long long time_s) const {
        if (control.trigger == ControlTrigger::ElapsedTime)
            return time_s == control.event_time_s;
        if (control.trigger == ControlTrigger::ClockTime)
            return (start_clock_s_ + time_s) % (24 * 3600) == control.event_time_s;
        return false;
    }

    long long scan_rules(const std::map<std::string, Agelem *> &links,
                         const std::map<std::string, Csomopont *> &nodes,
                         const std::vector<TankState> &tanks,
                         long long time_s, long long candidate_time_s) {
        long long previous_check_s = time_s + 1;
        long long check_s = std::min(
            (time_s / rule_step_s_ + 1) * rule_step_s_, candidate_time_s);
        while (check_s > time_s && check_s <= candidate_time_s) {
            std::map<std::string, double> projected_levels;
            for (const TankState &tank : tanks) {
                const double projected = tank.level_m + tank.boundary->Get_Q() *
                    static_cast<double>(check_s - time_s) / tank.area_m2;
                projected_levels[tank.id] = std::max(
                    tank.min_level_m, std::min(tank.max_level_m, projected));
            }

            const EpanetRuleEngine::NumericLookup numeric_lookup =
                [&](const RulePremise &premise, double &value) -> bool {
                    if (premise.object == "SYSTEM" && premise.variable == RuleVariable::Demand) {
                        value = 0.0;
                        for (const auto &node : nodes)
                            value += node.second->Get_dprop("demand");
                        return true;
                    }
                    if (premise.variable == RuleVariable::Demand ||
                        premise.variable == RuleVariable::Head ||
                        premise.variable == RuleVariable::Pressure ||
                        premise.variable == RuleVariable::Level ||
                        premise.variable == RuleVariable::FillTime ||
                        premise.variable == RuleVariable::DrainTime) {
                        const auto node = nodes.find(premise.id);
                        if (node == nodes.end()) {
                            warn("RULES", premise.id, premise.line_number,
                                 "Referenced rule node was not found.");
                            return false;
                        }
                        const auto level = projected_levels.find(premise.id);
                        if (premise.variable == RuleVariable::Demand)
                            value = node->second->Get_dprop("demand");
                        else if (premise.variable == RuleVariable::Level)
                            value = level == projected_levels.end() ? node->second->Get_p() : level->second;
                        else if (premise.variable == RuleVariable::Head)
                            value = node->second->Get_h() +
                                (level == projected_levels.end() ? node->second->Get_p() : level->second);
                        else if (premise.variable == RuleVariable::Pressure)
                            value = level == projected_levels.end() ? node->second->Get_p() : level->second;
                        else {
                            const auto tank = std::find_if(tanks.begin(), tanks.end(),
                                [&](const TankState &item) { return item.id == premise.id; });
                            if (tank == tanks.end())
                                return false;
                            const double flow = tank->boundary->Get_Q();
                            if (premise.variable == RuleVariable::FillTime) {
                                // EPANET treats FILLTIME as an inapplicable
                                // (false) premise unless the tank is filling.
                                if (flow <= 1.0e-12)
                                    return false;
                                value = (tank->max_level_m - level->second) *
                                    tank->area_m2 / flow;
                            } else {
                                // Likewise DRAINTIME is false while a tank is
                                // stationary or filling; infinity would make
                                // an ABOVE comparison incorrectly true.
                                if (flow >= -1.0e-12)
                                    return false;
                                value = (level->second - tank->min_level_m) *
                                    tank->area_m2 / -flow;
                            }
                        }
                        return true;
                    }
                    const auto link = links.find(premise.id);
                    if (link == links.end()) {
                        warn("RULES", premise.id, premise.line_number,
                             "Referenced rule link was not found.");
                        return false;
                    }
                    if (premise.variable == RuleVariable::Flow) {
                        value = std::abs(link->second->Get_Q());
                        return true;
                    }
                    if (premise.variable == RuleVariable::Setting) {
                        auto *pump = dynamic_cast<EpanetPumpConfigurable *>(link->second);
                        auto *valve = dynamic_cast<JelleggorbesFojtas *>(link->second);
                        if (pump == nullptr && valve == nullptr) {
                            warn("RULES", premise.id, premise.line_number,
                                 "SETTING premises currently require an EPANET pump or TCV.");
                            return false;
                        }
                        value = pump != nullptr ? pump->GetOperatingSpeed()
                                                : valve->GetEpanetTcvSetting();
                        return true;
                    }
                    return false;
                };

            const EpanetRuleEngine::StatusLookup status_lookup =
                [&](const RulePremise &premise, RuleStatus &value) -> bool {
                    const auto link = links.find(premise.id);
                    if (link == links.end()) {
                        warn("RULES", premise.id, premise.line_number,
                             "Referenced rule link was not found.");
                        return false;
                    }
                    auto *valve = dynamic_cast<JelleggorbesFojtas *>(link->second);
                    if (valve != nullptr) {
                        const EpanetTcvStatus status = valve->GetEpanetTcvStatus();
                        value = status == EpanetTcvStatus::Active ? RuleStatus::Active
                            : (status == EpanetTcvStatus::Open ? RuleStatus::Open
                                                               : RuleStatus::Closed);
                    } else {
                        value = link->second->Is_enabled() ? RuleStatus::Open
                                                           : RuleStatus::Closed;
                    }
                    return true;
                };

            const std::vector<RuleAction> actions = rule_engine_.select_actions(
                numeric_lookup, status_lookup, previous_check_s, check_s, start_clock_s_);
            bool changed = false;
            for (const RuleAction &action : actions) {
                const auto target = links.find(action.link_id);
                auto *tcv = target == links.end()
                    ? nullptr : dynamic_cast<JelleggorbesFojtas *>(target->second);
                // EPANET's rule executor applies OPEN only when a link is
                // currently closed. Unlike a simple control, OPEN does not
                // change an already ACTIVE TCV into its fully open state.
                if (action.action == ControlAction::Open && tcv != nullptr &&
                    tcv->GetEpanetTcvStatus() == EpanetTcvStatus::Active)
                    continue;
                const SimpleControl control{
                    action.link_id, action.action, action.setting,
                    ControlTrigger::ElapsedTime, "", 0.0, check_s,
                    action.line_number};
                changed = apply_control_action(control, links, "RULES") || changed;
            }
            if (changed)
                return check_s;
            if (check_s == candidate_time_s)
                break;
            previous_check_s = check_s + 1;
            check_s = std::min(check_s + rule_step_s_, candidate_time_s);
        }
        return candidate_time_s;
    }

    void apply_time_controls(const std::map<std::string, Agelem *> &links,
                             long long time_s) {
        for (const SimpleControl &control : controls_)
            if (time_control_is_active(control, time_s))
                apply_control_action(control, links);
    }

    long long next_time_control_event(long long time_s) const {
        long long next = duration_s_ + 1;
        for (const SimpleControl &control : controls_) {
            if (control.trigger == ControlTrigger::ElapsedTime) {
                if (control.event_time_s > time_s)
                    next = std::min(next, control.event_time_s);
            } else if (control.trigger == ControlTrigger::ClockTime) {
                const long long day_s = 24 * 3600;
                const long long clock_now = (start_clock_s_ + time_s) % day_s;
                long long delta = (control.event_time_s - clock_now + day_s) % day_s;
                if (delta == 0)
                    delta = day_s;
                next = std::min(next, time_s + delta);
            }
        }
        return next;
    }

    bool apply_node_controls(const std::map<std::string, Agelem *> &links,
                             const std::map<std::string, Csomopont *> &nodes,
                             const std::vector<TankState> &tanks) {
        std::map<std::string, double> levels;
        for (const TankState &tank : tanks)
            levels[tank.id] = tank.level_m;

        std::map<std::string, std::tuple<bool, double, int> > state_before;
        for (const SimpleControl &control : controls_) {
            const auto link = links.find(control.link_id);
            if (link == links.end())
                continue;
            auto *pump = dynamic_cast<EpanetPumpConfigurable *>(link->second);
            auto *valve = dynamic_cast<JelleggorbesFojtas *>(link->second);
            const double setting = pump != nullptr ? pump->GetOperatingSpeed()
                : (valve != nullptr ? valve->GetEpanetTcvSetting() : 0.0);
            const int valve_status = valve == nullptr
                ? -1 : static_cast<int>(valve->GetEpanetTcvStatus());
            state_before[control.link_id] = std::make_tuple(
                link->second->Is_enabled(),
                setting, valve_status);
        }

        for (const SimpleControl &control : controls_) {
            if (control.trigger != ControlTrigger::NodeAbove &&
                control.trigger != ControlTrigger::NodeBelow)
                continue;

            double value = 0.0;
            const auto level = levels.find(control.node_id);
            if (level != levels.end()) {
                value = level->second;
            } else {
                const auto node = nodes.find(control.node_id);
                if (node == nodes.end()) {
                    warn("CONTROLS", control.link_id, control.line_number,
                         "Control node '" + control.node_id + "' was not found.");
                    continue;
                }
                value = node->second->Get_p();
            }

            const bool condition = control.trigger == ControlTrigger::NodeBelow
                ? value < control.threshold_si : value > control.threshold_si;
            if (condition)
                apply_control_action(control, links);
        }

        for (const auto &old_state : state_before) {
            const auto link = links.find(old_state.first);
            if (link == links.end())
                continue;
            auto *pump = dynamic_cast<EpanetPumpConfigurable *>(link->second);
            auto *valve = dynamic_cast<JelleggorbesFojtas *>(link->second);
            const double setting = pump != nullptr ? pump->GetOperatingSpeed()
                : (valve != nullptr ? valve->GetEpanetTcvSetting() : 0.0);
            const int valve_status = valve == nullptr
                ? -1 : static_cast<int>(valve->GetEpanetTcvStatus());
            if (link->second->Is_enabled() != std::get<0>(old_state.second) ||
                std::abs(setting - std::get<1>(old_state.second)) > 1.0e-12 ||
                valve_status != std::get<2>(old_state.second))
                return true;
        }
        return false;
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
                                 long long time_s, bool converged,
                                 const std::vector<double> &node_age_s,
                                 const std::vector<double> &link_age_s,
                                 const std::vector<double> &node_chlorine_kgm3,
                                 const std::vector<double> &link_chlorine_kgm3) {
        EpsResultFrame frame;
        frame.time_s = time_s;
        frame.converged = converged;
        frame.iterations = 0; // The legacy solver does not expose this counter yet.
        const double missing = std::numeric_limits<double>::quiet_NaN();
        std::map<std::string, double> total_head;
        for (std::size_t index = 0; index < system.cspok.size(); ++index) {
            Csomopont *node = system.cspok[index];
            const double head = node->Get_p() + node->Get_h();
            total_head[node->Get_nev()] = head;
            frame.node_head_m.push_back(head);
            frame.node_pressure_head_m.push_back(head - node_info[index].elevation_m);
            frame.node_demand_m3s.push_back(node->Get_dprop("demand") / 3600.0);
            frame.node_water_age_s.push_back(
                node_age_s.size() == system.cspok.size() ? node_age_s[index] : missing);
            frame.node_chlorine_kgm3.push_back(
                node_chlorine_kgm3.size() == system.cspok.size()
                    ? node_chlorine_kgm3[index] : missing);
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
            if (is_pump(link))
                frame.link_status.push_back(link->Get_dprop("status") != 0.0 ? 1 : 0);
            else
                frame.link_status.push_back(link->Is_enabled() ? 1 : 0);
            frame.link_water_age_s.push_back(
                link_age_s.size() == links.size() ? link_age_s[index] : missing);
            frame.link_chlorine_kgm3.push_back(
                link_chlorine_kgm3.size() == links.size()
                    ? link_chlorine_kgm3[index] : missing);
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
                        << (frame.converged ? 1 : 0) << ','
                        << frame.node_water_age_s[index] << ','
                        << frame.node_chlorine_kgm3[index] << '\n';
        for (std::size_t index = 0; index < links.size(); ++index)
            link_output << frame.time_s << ',' << csv(links[index].id) << ','
                        << csv(links[index].type) << ','
                        << csv(nodes[links[index].from_node].id) << ','
                        << csv(nodes[links[index].to_node].id) << ','
                        << frame.link_flow_rate_m3s[index] << ','
                        << frame.link_velocity_ms[index] << ','
                        << frame.link_headloss_m[index] << ','
                        << static_cast<unsigned int>(frame.link_status[index]) << ','
                        << (frame.converged ? 1 : 0) << ','
                        << frame.link_water_age_s[index] << ','
                        << frame.link_chlorine_kgm3[index] << '\n';
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

EpanetSteadyQualitySimulation::EpanetSteadyQualitySimulation(
    const std::string &input_filename, const std::string &output_prefix,
    const std::string &quality_mode, bool compute_sensitivity)
    : input_filename_(input_filename), output_prefix_(output_prefix),
      quality_mode_(quality_mode), compute_sensitivity_(compute_sensitivity) {}

void EpanetSteadyQualitySimulation::run(Staci &system) {
    if (output_prefix_.empty())
        throw std::runtime_error(
            "Steady EPANET water quality requires -o <result-prefix>.");
    const std::filesystem::path prefix_path(output_prefix_);
    if (!prefix_path.parent_path().empty())
        std::filesystem::create_directories(prefix_path.parent_path());
    SimulationModel model(input_filename_);
    model.run_steady_quality(system, output_prefix_, quality_mode_,
                             compute_sensitivity_);
}
