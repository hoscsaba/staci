#include "epanet_reader.h"

#include "Cso.h"
#include "Csomopont.h"
#include "EpanetPowerPump.h"
#include "KonstNyomas.h"
#include "Szivattyu.h"
#include "Vegakna.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {

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
    std::vector<std::string> fields;
    std::istringstream input(line);
    std::string field;
    while (input >> field)
        fields.push_back(field);
    return fields;
}

bool starts_with(const std::vector<std::string> &fields,
                 const std::vector<std::string> &prefix) {
    if (fields.size() <= prefix.size())
        return false;
    for (std::size_t i = 0; i < prefix.size(); ++i)
        if (upper(fields[i]) != prefix[i])
            return false;
    return true;
}

std::string to_string_precise(double value) {
    std::ostringstream output;
    output << std::setprecision(15) << value;
    return output.str();
}

} // namespace

EpanetReader::EpanetReader(const std::string &filename, bool extended_period)
    : filename_(filename),
      flow_units_("LPS"),
      headloss_model_("H-W"),
      quality_mode_("NONE"),
      specific_gravity_(1.0),
      demand_multiplier_(1.0),
      us_customary_units_(false),
      extended_period_(extended_period) {
    parse_file();
    parse_options();
    parse_patterns();
    parse_curves();
    configure_settings();
}

void EpanetReader::parse_file() {
    document_ = EpanetDocument::read(filename_);
    for (const EpanetDocumentLine &line : document_.lines()) {
        if (line.section_header) {
            sections_[line.section];
            continue;
        }
        if (!line.data_record)
            continue;
        if (line.section.empty()) {
            add_warning("GLOBAL", "", line.line_number,
                        "Data before the first section was ignored.");
            continue;
        }
        sections_[line.section].push_back(
            Record{line.line_number, line.content, split_fields(line.content)});
    }

    if (sections_.find("END") == sections_.end())
        add_warning("END", "", 0, "The file has no [END] section.");
}

void EpanetReader::parse_options() {
    for (const Record &record : records("OPTIONS")) {
        if (record.fields.size() < 2) {
            add_warning("OPTIONS", "", record.line_number,
                        "Incomplete option was parsed but ignored.");
            continue;
        }
        const std::string key = upper(record.fields[0]);
        if (key == "UNITS")
            flow_units_ = upper(record.fields[1]);
        else if (key == "HEADLOSS")
            headloss_model_ = upper(record.fields[1]);
        else if (key == "PATTERN")
            default_pattern_ = record.fields[1];
        else if (key == "TRIALS")
            settings_["iter_max"] = record.fields[1];
        else if (key == "ACCURACY")
            settings_["e_p_max"] = record.fields[1];
        else if (key == "QUALITY") {
            quality_mode_ = upper(record.fields[1]);
            if (quality_mode_ == "AGE")
                add_warning("OPTIONS", "QUALITY", record.line_number,
                            "Initial water age is transferred, but EPANET extended-period "
                            "quality simulation is not executed.");
            else if (quality_mode_ == "CHEMICAL")
                add_warning("OPTIONS", "QUALITY", record.line_number,
                            "Initial chemical concentration is transferred; EPANET-specific "
                            "quality sources and reactions are not.");
            else if (quality_mode_ != "NONE")
                add_warning("OPTIONS", "QUALITY", record.line_number,
                            "Quality mode '" + record.fields[1] +
                            "' has no STACI equivalent and was not transferred.");
        }
        else if (starts_with(record.fields, {"SPECIFIC", "GRAVITY"}))
            specific_gravity_ = number(record, 2, "specific gravity", 1.0);
        else if (starts_with(record.fields, {"DEMAND", "MULTIPLIER"}))
            demand_multiplier_ = number(record, 2, "demand multiplier", 1.0);
        // Other OPTIONS records remain available in the lossless EPANET
        // document even when they do not affect the STACI hydraulic snapshot.
    }

    static const std::set<std::string> us_units = {"CFS", "GPM", "MGD", "IMGD", "AFD"};
    static const std::set<std::string> metric_units = {"LPS", "LPM", "MLD", "CMH", "CMD"};
    us_customary_units_ = us_units.count(flow_units_) != 0;
    if (!us_customary_units_ && metric_units.count(flow_units_) == 0) {
        add_warning("OPTIONS", "UNITS", 0,
                    "Unknown flow units '" + flow_units_ + "'; LPS was assumed.");
        flow_units_ = "LPS";
    }

    if (specific_gravity_ <= 0.0) {
        add_warning("OPTIONS", "SPECIFIC GRAVITY", 0,
                    "Specific gravity must be positive; 1.0 was used.");
        specific_gravity_ = 1.0;
    }
}

void EpanetReader::parse_patterns() {
    for (const Record &record : records("PATTERNS")) {
        if (record.fields.size() < 2) {
            add_warning("PATTERNS", record.fields.empty() ? "" : record.fields[0],
                        record.line_number, "Incomplete pattern row was ignored.");
            continue;
        }
        std::vector<double> &values = patterns_[record.fields[0]];
        for (std::size_t i = 1; i < record.fields.size(); ++i)
            values.push_back(number(record, i, "pattern multiplier", 1.0));
    }
}

void EpanetReader::parse_curves() {
    for (const Record &record : records("CURVES")) {
        if (record.fields.size() < 3) {
            add_warning("CURVES", record.fields.empty() ? "" : record.fields[0],
                        record.line_number, "Incomplete curve row was ignored.");
            continue;
        }
        curves_[record.fields[0]].push_back(
            std::make_pair(number(record, 1, "curve x value", 0.0),
                           number(record, 2, "curve y value", 0.0)));
    }
}

void EpanetReader::configure_settings() {
    settings_["debug_level"] = "1";
    settings_["out_file"] = "staci.out";
    if (settings_.find("iter_max") == settings_.end())
        settings_["iter_max"] = "100";
    if (settings_.find("e_p_max") == settings_.end())
        settings_["e_p_max"] = "0.001";
    else
        add_warning("OPTIONS", "ACCURACY", 0,
                    "EPANET's relative accuracy criterion is mapped approximately "
                    "to STACI's pressure tolerance.");
    settings_["e_mp_max"] = "1.0";
    settings_["relax"] = "1.0";
    settings_["relax_mul"] = "1.2";
    settings_["mp_init"] = "1.0";
    settings_["p_init"] = "100.0";

    std::string normalized = headloss_model_;
    normalized.erase(std::remove(normalized.begin(), normalized.end(), '-'), normalized.end());
    if (normalized == "HW")
        settings_["friction_model"] = "HW";
    else if (normalized == "DW")
        settings_["friction_model"] = "DW";
    else {
        settings_["friction_model"] = "DW";
        add_warning("OPTIONS", "HEADLOSS", 0,
                    "STACI does not implement the EPANET " + headloss_model_ +
                    " pipe formula; Darcy-Weisbach was selected.");
    }
}

const std::vector<EpanetReader::Record> &
EpanetReader::records(const std::string &section) const {
    static const std::vector<Record> empty;
    const auto found = sections_.find(section);
    return found == sections_.end() ? empty : found->second;
}

bool EpanetReader::has_records(const std::string &section) const {
    return !records(section).empty();
}

void EpanetReader::add_warning(const std::string &section,
                               const std::string &element,
                               std::size_t line_number,
                               const std::string &message) {
    warnings_.push_back(Warning{section, element, line_number, message});
}

double EpanetReader::number(const Record &record,
                            std::size_t field,
                            const std::string &label,
                            double default_value) {
    if (field >= record.fields.size()) {
        add_warning("PARSE", record.fields.empty() ? "" : record.fields[0],
                    record.line_number, "Missing " + label + "; default value used.");
        return default_value;
    }
    char *end = nullptr;
    const double value = std::strtod(record.fields[field].c_str(), &end);
    if (end == record.fields[field].c_str() || *end != '\0') {
        add_warning("PARSE", record.fields[0], record.line_number,
                    "Invalid " + label + " '" + record.fields[field] +
                    "'; default value used.");
        return default_value;
    }
    return value;
}

double EpanetReader::flow_to_m3_per_hour(double value) const {
    if (flow_units_ == "CFS") return value * 101.9406477312;
    if (flow_units_ == "GPM") return value * 0.22712470704;
    if (flow_units_ == "MGD") return value * 157.725491;
    if (flow_units_ == "IMGD") return value * 189.420142;
    if (flow_units_ == "AFD") return value * 51.3950766;
    if (flow_units_ == "LPS") return value * 3.6;
    if (flow_units_ == "LPM") return value * 0.06;
    if (flow_units_ == "MLD") return value * 41.6666666667;
    if (flow_units_ == "CMH") return value;
    if (flow_units_ == "CMD") return value / 24.0;
    return value * 3.6;
}

double EpanetReader::length_to_metres(double value) const {
    return us_customary_units_ ? value * 0.3048 : value;
}

double EpanetReader::pipe_diameter_to_metres(double value) const {
    return us_customary_units_ ? value * 0.0254 : value / 1000.0;
}

double EpanetReader::tank_diameter_to_metres(double value) const {
    return us_customary_units_ ? value * 0.3048 : value;
}

double EpanetReader::pump_power_to_watts(double value) const {
    return us_customary_units_ ? value * 745.699871582 : value * 1000.0;
}

double EpanetReader::first_pattern_multiplier(const std::string &pattern_id) const {
    if (pattern_id.empty())
        return 1.0;
    const auto found = patterns_.find(pattern_id);
    if (found == patterns_.end() || found->second.empty())
        return 1.0;
    return found->second.front();
}

void EpanetReader::load_system(std::vector<std::unique_ptr<Csomopont> > &nodes,
                               std::vector<std::unique_ptr<Agelem> > &edges) {
    const double density = 1000.0 * specific_gravity_;
    std::set<std::string> node_ids;
    std::set<std::string> edge_ids;
    std::map<std::string, double> extra_demands;
    std::map<std::string, double> initial_quality;

    const auto has_fields = [this](const Record &record, std::size_t count,
                                   const std::string &section) {
        if (record.fields.size() >= count)
            return true;
        add_warning(section, record.fields.empty() ? "" : record.fields[0],
                    record.line_number, "Incomplete row was ignored.");
        return false;
    };

    for (const Record &record : records("DEMANDS")) {
        if (!has_fields(record, 2, "DEMANDS"))
            continue;
        const std::string pattern = record.fields.size() > 2 ? record.fields[2] : default_pattern_;
        extra_demands[record.fields[0]] +=
            flow_to_m3_per_hour(number(record, 1, "demand", 0.0)) *
            first_pattern_multiplier(pattern) * demand_multiplier_;
    }

    for (const Record &record : records("QUALITY"))
        if (has_fields(record, 2, "QUALITY"))
            initial_quality[record.fields[0]] = number(record, 1, "initial quality", 0.0);

    for (const Record &record : records("JUNCTIONS")) {
        if (!has_fields(record, 2, "JUNCTIONS"))
            continue;
        const std::string &id = record.fields[0];
        if (!node_ids.insert(id).second) {
            add_warning("JUNCTIONS", id, record.line_number,
                        "Duplicate node ID; later definition ignored.");
            continue;
        }
        const double elevation = length_to_metres(number(record, 1, "elevation", 0.0));
        const double base_demand = record.fields.size() > 2
            ? number(record, 2, "base demand", 0.0) : 0.0;
        const std::string pattern = record.fields.size() > 3
            ? record.fields[3] : default_pattern_;
        const double demand = flow_to_m3_per_hour(base_demand) *
            first_pattern_multiplier(pattern) * demand_multiplier_ + extra_demands[id];
        const auto quality = initial_quality.find(id);
        const double initial_age = quality_mode_ == "AGE" && quality != initial_quality.end()
            ? quality->second : 0.0;
        auto node = std::make_unique<Csomopont>(id, elevation, demand, 0.0, 0.0,
                                                density, initial_age);
        if (quality_mode_ == "CHEMICAL" && quality != initial_quality.end())
            node->Set_dprop("concentration", quality->second);
        nodes.push_back(std::move(node));
    }

    double fixed_head_sum = 0.0;
    std::size_t fixed_head_count = 0;

    for (const Record &record : records("RESERVOIRS")) {
        if (!has_fields(record, 2, "RESERVOIRS"))
            continue;
        const std::string &id = record.fields[0];
        if (!node_ids.insert(id).second) {
            add_warning("RESERVOIRS", id, record.line_number,
                        "Duplicate node ID; reservoir ignored.");
            continue;
        }
        const std::string pattern = record.fields.size() > 2 ? record.fields[2] : "";
        const double head = length_to_metres(number(record, 1, "head", 0.0)) *
            first_pattern_multiplier(pattern);
        const auto quality = initial_quality.find(id);
        const double initial_age = quality_mode_ == "AGE" && quality != initial_quality.end()
            ? quality->second : 0.0;
        auto node = std::make_unique<Csomopont>(id, 0.0, 0.0, 0.0, head, density, initial_age);
        if (quality_mode_ == "CHEMICAL" && quality != initial_quality.end())
            node->Set_dprop("concentration", quality->second);
        nodes.push_back(std::move(node));
        const std::string boundary_id = "EPANET_RESERVOIR_" + id;
        edge_ids.insert(boundary_id);
        edges.push_back(std::make_unique<KonstNyomas>(
            boundary_id, 1.0, id, density, head * density * 9.81, 1.0, 0.0));
        fixed_head_sum += head;
        ++fixed_head_count;
    }

    for (const Record &record : records("TANKS")) {
        if (!has_fields(record, 6, "TANKS"))
            continue;
        const std::string &id = record.fields[0];
        if (!node_ids.insert(id).second) {
            add_warning("TANKS", id, record.line_number,
                        "Duplicate node ID; tank ignored.");
            continue;
        }
        const double elevation = length_to_metres(number(record, 1, "elevation", 0.0));
        const double initial_level = length_to_metres(number(record, 2, "initial level", 0.0));
        const double diameter = tank_diameter_to_metres(number(record, 5, "diameter", 1.0));
        const double area = 3.14159265358979323846 * diameter * diameter / 4.0;
        const auto quality = initial_quality.find(id);
        const double initial_age = quality_mode_ == "AGE" && quality != initial_quality.end()
            ? quality->second : 0.0;
        auto node = std::make_unique<Csomopont>(id, 0.0, 0.0, 0.0,
                                                elevation + initial_level, density, initial_age);
        if (quality_mode_ == "CHEMICAL" && quality != initial_quality.end())
            node->Set_dprop("concentration", quality->second);
        nodes.push_back(std::move(node));
        const std::string boundary_id = "EPANET_TANK_" + id;
        edge_ids.insert(boundary_id);
        edges.push_back(std::make_unique<Vegakna>(
            boundary_id, id, density, area, elevation, initial_level, 1.0, 0.0));
        if (!extended_period_)
            add_warning("TANKS", id, record.line_number,
                        "Imported as a fixed boundary at its initial level; min/max level, "
                        "volume curve, and dynamic storage are not simulated by steady-state STACI.");
        fixed_head_sum += elevation + initial_level;
        ++fixed_head_count;
    }

    if (fixed_head_count != 0)
        settings_["p_init"] = to_string_precise(fixed_head_sum / fixed_head_count);

    for (const Record &record : records("PIPES")) {
        if (!has_fields(record, 6, "PIPES"))
            continue;
        const std::string &id = record.fields[0];
        const std::string &from = record.fields[1];
        const std::string &to = record.fields[2];
        if (node_ids.count(from) == 0 || node_ids.count(to) == 0) {
            add_warning("PIPES", id, record.line_number,
                        "Unknown endpoint; pipe was not imported.");
            continue;
        }
        if (!edge_ids.insert(id).second) {
            add_warning("PIPES", id, record.line_number,
                        "Duplicate link ID; pipe was not imported.");
            continue;
        }
        const double length = length_to_metres(number(record, 3, "length", 1.0));
        const double diameter = pipe_diameter_to_metres(number(record, 4, "diameter", 0.1));
        double roughness = number(record, 5, "roughness", 100.0);
        std::string normalized = headloss_model_;
        normalized.erase(std::remove(normalized.begin(), normalized.end(), '-'), normalized.end());
        if (normalized == "DW" && us_customary_units_)
            roughness *= 0.3048; // EPANET US D-W roughness is in 0.001 ft; STACI uses mm.
        if (record.fields.size() > 6 && number(record, 6, "minor loss", 0.0) != 0.0)
            add_warning("PIPES", id, record.line_number,
                        "Minor-loss coefficient is not represented and was ignored.");
        if (record.fields.size() > 7 && upper(record.fields[7]) != "OPEN")
            add_warning("PIPES", id, record.line_number,
                        "Initial status '" + record.fields[7] +
                        "' is not represented; pipe was imported open.");
        edges.push_back(std::make_unique<Cso>(id, from, to, density, length, diameter,
                                              roughness, 0.0, 0.0, 1.0));
    }

    for (const Record &record : records("PUMPS")) {
        if (!has_fields(record, 5, "PUMPS"))
            continue;
        const std::string &id = record.fields[0];
        const std::string &from = record.fields[1];
        const std::string &to = record.fields[2];
        if (node_ids.count(from) == 0 || node_ids.count(to) == 0) {
            add_warning("PUMPS", id, record.line_number,
                        "Unknown endpoint; pump was not imported.");
            continue;
        }
        if (!edge_ids.insert(id).second) {
            add_warning("PUMPS", id, record.line_number,
                        "Duplicate link ID; pump was not imported.");
            continue;
        }
        const std::string type = upper(record.fields[3]);
        if (type == "POWER") {
            const double power = pump_power_to_watts(number(record, 4, "power", 0.0));
            edges.push_back(std::make_unique<EpanetPowerPump>(id, from, to, density, power, 1.0));
        } else if (type == "HEAD") {
            const auto curve = curves_.find(record.fields[4]);
            if (curve == curves_.end() || curve->second.size() < 3) {
                add_warning("PUMPS", id, record.line_number,
                            "Pump curve is missing or has fewer than three points; pump was not imported.");
                continue;
            }
            std::vector<double> flow;
            std::vector<double> head;
            for (const auto &point : curve->second) {
                flow.push_back(flow_to_m3_per_hour(point.first));
                head.push_back(length_to_metres(point.second));
            }
            edges.push_back(std::make_unique<Szivattyu>(id, from, to, density, 1.0, flow, head, 1.0));
        } else {
            add_warning("PUMPS", id, record.line_number,
                        "Unsupported pump parameter syntax; pump was not imported.");
        }
        if (record.fields.size() > 5)
            add_warning("PUMPS", id, record.line_number,
                        "Additional pump speed/pattern parameters are not represented.");
    }

    for (const Record &record : records("VALVES")) {
        const std::string id = record.fields.empty() ? "" : record.fields[0];
        add_warning("VALVES", id, record.line_number,
                    "EPANET valve types/settings have no equivalent STACI element; valve was not imported.");
    }

    if (!extended_period_)
        for (const Record &record : records("CONTROLS"))
            add_warning("CONTROLS", record.fields.size() > 1 ? record.fields[1] : "",
                        record.line_number,
                        "Dynamic link control was parsed but is not executed by steady-state STACI.");

    warn_for_unsupported_sections();
    print_report(nodes.size(), edges.size());
}

void EpanetReader::warn_for_unsupported_sections() {
    const std::vector<std::pair<std::string, std::string> > unsupported = {
        {"RULES", "Rule-based controls are not executed by STACI."},
        {"STATUS", "Initial link status/settings are not represented by the steady-state importer."},
        {"EMITTERS", "Pressure-dependent emitters are not represented."},
        {"SOURCES", "EPANET water-quality sources are not represented."},
        {"REACTIONS", "EPANET reaction options are not transferred to STACI's transport model."},
        {"MIXING", "Tank mixing models are not represented."},
    };
    for (const auto &entry : unsupported)
        if (has_records(entry.first) &&
            !(extended_period_ && (entry.first == "TIMES" || entry.first == "STATUS" ||
                                   entry.first == "COORDINATES")))
            add_warning(entry.first, "", records(entry.first).front().line_number,
                        entry.second);

    if (has_records("PATTERNS") && !extended_period_)
        add_warning("PATTERNS", "", records("PATTERNS").front().line_number,
                    "Only the first multiplier is applied to the imported steady-state snapshot.");

    if (has_records("DEMANDS"))
        for (const Record &record : records("DEMANDS"))
            if (record.fields.size() > 3) {
                add_warning("DEMANDS", "", record.line_number,
                            "Demand category labels are not stored; their demands were summed.");
                break;
            }

    if (has_records("QUALITY") && quality_mode_ != "AGE" && quality_mode_ != "CHEMICAL")
        add_warning("QUALITY", "", records("QUALITY").front().line_number,
                    "Initial quality values were parsed but cannot be represented for quality mode '" +
                    quality_mode_ + "'.");

    static const std::set<std::string> known = {
        "TITLE", "JUNCTIONS", "RESERVOIRS", "TANKS", "PIPES", "PUMPS",
        "VALVES", "TAGS", "DEMANDS", "STATUS", "PATTERNS", "CURVES",
        "CONTROLS", "RULES", "ENERGY", "EMITTERS", "QUALITY", "SOURCES",
        "REACTIONS", "MIXING", "TIMES", "REPORT", "OPTIONS", "COORDINATES",
        "VERTICES", "LABELS", "BACKDROP", "END"
    };
    for (const auto &section : sections_)
        if (known.count(section.first) == 0)
            add_warning(section.first, "",
                        section.second.empty() ? 0 : section.second.front().line_number,
                        "Unknown EPANET section was parsed but not imported.");
}

void EpanetReader::print_report(std::size_t node_count,
                                std::size_t edge_count) const {
    std::cerr << "\nEPANET import: " << filename_ << "\n"
              << "  units: " << flow_units_ << "\n"
              << "  imported STACI nodes: " << node_count << "\n"
              << "  imported STACI elements: " << edge_count << "\n"
              << "  compatibility warnings: " << warnings_.size() << "\n";
    for (const Warning &warning : warnings_) {
        std::cerr << "WARNING [EPANET][" << warning.section << "]";
        if (!warning.element.empty())
            std::cerr << "[" << warning.element << "]";
        if (warning.line_number != 0)
            std::cerr << " line " << warning.line_number;
        std::cerr << ": " << warning.message << "\n";
    }
    std::cerr << std::endl;
}

std::string EpanetReader::read_setting(const std::string &name) const {
    const auto found = settings_.find(name);
    return found == settings_.end() ? "nincs ilyen node!" : found->second;
}
