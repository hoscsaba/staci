#include "epanet_reader.h"

#include "Cso.h"
#include "Csomopont.h"
#include "EpanetPowerPump.h"
#include "KonstNyomas.h"
#include "JelleggorbesFojtas.h"
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

std::string join_fields(const std::vector<std::string> &fields,
                        std::size_t begin) {
    std::ostringstream output;
    for (std::size_t i = begin; i < fields.size(); ++i) {
        if (i != begin)
            output << ' ';
        output << fields[i];
    }
    return output.str();
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
            quality_chemical_name_.clear();
            quality_units_.clear();
            quality_trace_node_.clear();
            if (quality_mode_ == "AGE") {
                quality_units_ = "h";
                add_warning("OPTIONS", "QUALITY", record.line_number,
                            "Initial water age is transferred, but EPANET extended-period "
                            "quality simulation is not executed.");
            } else if (quality_mode_ == "CHEMICAL") {
                if (record.fields.size() > 2)
                    quality_chemical_name_ = record.fields[2];
                if (record.fields.size() > 3)
                    quality_units_ = record.fields[3];
                add_warning("OPTIONS", "QUALITY", record.line_number,
                            "Initial chemical concentration is transferred; EPANET-specific "
                            "quality sources and reactions are not.");
            } else if (quality_mode_ == "TRACE") {
                if (record.fields.size() > 2)
                    quality_trace_node_ = record.fields[2];
                quality_units_ = "%";
                add_warning("OPTIONS", "QUALITY", record.line_number,
                            "Trace initial values are retained as node metadata, but "
                            "EPANET trace simulation is not executed.");
            } else if (quality_mode_ != "NONE") {
                quality_units_ = record.fields.size() > 2
                    ? join_fields(record.fields, 2) : "";
                add_warning("OPTIONS", "QUALITY", record.line_number,
                            "Quality mode '" + record.fields[1] +
                            "' is retained as node metadata but is not simulated.");
            }
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
                    "to STACI's pressure and mass-flow residual tolerances.");
    // EPANET ACCURACY is dimensionless, so neither STACI dimensional
    // tolerance is an exact equivalent. Reusing its numeric value for the
    // RMS mass-flow residual keeps changing EPS demands from being accepted
    // as converged solely because they differ by less than the old 1 kg/s
    // threshold.
    settings_["e_mp_max"] = settings_["e_p_max"];
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
    std::map<std::string, std::string> status_overrides;
    std::map<std::string, double> setting_overrides;
    std::set<std::string> handled_status_ids;
    EpanetPumpMetadata global_pump_energy;
    std::map<std::string, EpanetPumpMetadata> pump_energy;
    std::map<std::string, std::vector<EpanetDemandComponent> > additional_demands;
    std::map<std::string, EpanetInitialQuality> initial_quality;

    for (const Record &record : records("STATUS")) {
        if (record.fields.size() < 2)
            continue;
        const std::string status = upper(record.fields[1]);
        if (status == "OPEN" || status == "CLOSED" || status == "ACTIVE")
            status_overrides[record.fields[0]] = status;
        else {
            const double setting = number(record, 1, "link setting", -1.0);
            if (setting >= 0.0)
                setting_overrides[record.fields[0]] = setting;
            else
                add_warning("STATUS", record.fields[0], record.line_number,
                            "Negative link setting was not applied.");
        }
    }

    for (const Record &record : records("ENERGY")) {
        if (record.fields.size() < 2)
            continue;
        const std::string scope = upper(record.fields[0]);
        EpanetPumpMetadata *energy = nullptr;
        std::size_t key_index = 1;
        if (scope == "GLOBAL") {
            energy = &global_pump_energy;
        } else if (scope == "PUMP" && record.fields.size() >= 4) {
            energy = &pump_energy[record.fields[1]];
            key_index = 2;
        } else if (scope == "DEMAND" && record.fields.size() >= 3 &&
                   upper(record.fields[1]) == "CHARGE") {
            global_pump_energy.demand_charge_specified = true;
            global_pump_energy.demand_charge = number(record, 2, "demand charge", 0.0);
            continue;
        }
        if (energy == nullptr || key_index + 1 >= record.fields.size())
            continue;
        const std::string key = upper(record.fields[key_index]);
        const std::size_t value_index = key_index + 1;
        const bool global = scope == "GLOBAL";
        if (key == "PRICE") {
            energy->energy_price_specified = true;
            energy->energy_price_global = global;
            energy->energy_price = number(record, value_index, "energy price", 0.0);
        } else if (key == "PATTERN") {
            energy->energy_pattern_id = record.fields[value_index];
            energy->energy_pattern_global = global;
        } else if (key == "EFFIC" || key == "EFFICIENCY") {
            energy->efficiency_curve_id = record.fields[value_index];
            energy->efficiency_curve_global = global;
        }
    }

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
        const auto values = patterns_.find(pattern);
        additional_demands[record.fields[0]].push_back(EpanetDemandComponent{
            flow_to_m3_per_hour(number(record, 1, "demand", 0.0)) / 3600.0,
            pattern,
            values == patterns_.end() ? std::vector<double>() : values->second,
            record.fields.size() > 3 ? join_fields(record.fields, 3) : "",
            false});
    }

    for (const Record &record : records("QUALITY"))
        if (has_fields(record, 2, "QUALITY")) {
            EpanetInitialQuality quality;
            quality.specified = true;
            quality.source_value = number(record, 1, "initial quality", 0.0);
            quality.mode = quality_mode_;
            quality.chemical_name = quality_chemical_name_;
            quality.units = quality_units_;
            quality.trace_node = quality_trace_node_;
            initial_quality[record.fields[0]] = quality;
        }

    const auto quality_for_node = [this, &initial_quality](const std::string &id) {
        const auto found = initial_quality.find(id);
        if (found != initial_quality.end())
            return found->second;
        EpanetInitialQuality quality;
        quality.mode = quality_mode_;
        quality.chemical_name = quality_chemical_name_;
        quality.units = quality_units_;
        quality.trace_node = quality_trace_node_;
        return quality;
    };

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
        const auto pattern_values = patterns_.find(pattern);
        std::vector<EpanetDemandComponent> components;
        components.push_back(EpanetDemandComponent{
            flow_to_m3_per_hour(base_demand) / 3600.0,
            pattern,
            pattern_values == patterns_.end() ? std::vector<double>() : pattern_values->second,
            "",
            true});
        const auto additional = additional_demands.find(id);
        if (additional != additional_demands.end())
            components.insert(components.end(), additional->second.begin(), additional->second.end());
        const bool demand_section_overrides_primary =
            additional != additional_demands.end() && !additional->second.empty();
        double demand = 0.0;
        for (const EpanetDemandComponent &component : components) {
            // In EPANET, one or more [DEMANDS] records for a junction replace
            // the primary demand from [JUNCTIONS]. Preserve that primary
            // record as metadata for lossless export, but do not solve it.
            if (demand_section_overrides_primary && component.primary)
                continue;
            const double pattern_multiplier = component.pattern_values.empty()
                ? 1.0 : component.pattern_values.front();
            demand += component.base_demand_m3s * 3600.0 * pattern_multiplier;
        }
        demand *= demand_multiplier_;
        const EpanetInitialQuality quality = quality_for_node(id);
        const double initial_age = quality_mode_ == "AGE" && quality.specified
            ? quality.source_value : 0.0;
        auto node = std::make_unique<Csomopont>(id, elevation, demand, 0.0, 0.0,
                                                density, initial_age);
        node->SetEpanetDemandComponents(components, demand_multiplier_);
        node->SetEpanetInitialQuality(quality);
        if (quality_mode_ == "CHEMICAL" && quality.specified)
            node->Set_dprop("concentration", quality.source_value);
        nodes.push_back(std::move(node));
    }

    double fixed_head_sum = 0.0;
    std::size_t fixed_head_count = 0;
    std::size_t reservoir_head_pattern_count = 0;

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
        const double base_head = length_to_metres(number(record, 1, "head", 0.0));
        const double head = base_head * first_pattern_multiplier(pattern);
        const EpanetInitialQuality quality = quality_for_node(id);
        const double initial_age = quality_mode_ == "AGE" && quality.specified
            ? quality.source_value : 0.0;
        auto node = std::make_unique<Csomopont>(id, 0.0, 0.0, 0.0, head, density, initial_age);
        node->SetEpanetInitialQuality(quality);
        if (quality_mode_ == "CHEMICAL" && quality.specified)
            node->Set_dprop("concentration", quality.source_value);
        nodes.push_back(std::move(node));
        const std::string boundary_id = "EPANET_RESERVOIR_" + id;
        edge_ids.insert(boundary_id);
        auto boundary = std::make_unique<KonstNyomas>(
            boundary_id, 1.0, id, density, head * density * 9.81, 1.0, 0.0);
        const auto pattern_values = patterns_.find(pattern);
        boundary->SetEpanetHeadPattern(EpanetHeadPattern{
            true,
            base_head,
            pattern,
            pattern_values == patterns_.end()
                ? std::vector<double>() : pattern_values->second});
        if (!pattern.empty())
            ++reservoir_head_pattern_count;
        edges.push_back(std::move(boundary));
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
        const EpanetInitialQuality quality = quality_for_node(id);
        const double initial_age = quality_mode_ == "AGE" && quality.specified
            ? quality.source_value : 0.0;
        auto node = std::make_unique<Csomopont>(id, 0.0, 0.0, 0.0,
                                                elevation + initial_level, density, initial_age);
        node->SetEpanetInitialQuality(quality);
        if (quality_mode_ == "CHEMICAL" && quality.specified)
            node->Set_dprop("concentration", quality.source_value);
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
        double minor_loss = record.fields.size() > 6
            ? number(record, 6, "minor loss", 0.0) : 0.0;
        if (minor_loss < 0.0) {
            add_warning("PIPES", id, record.line_number,
                        "Negative minor-loss coefficient was replaced with zero.");
            minor_loss = 0.0;
        }
        const std::string inline_status = record.fields.size() > 7
            ? upper(record.fields[7]) : "OPEN";
        const bool check_valve = inline_status == "CV";
        bool enabled = inline_status != "CLOSED";
        if (inline_status != "OPEN" && inline_status != "CLOSED" &&
            inline_status != "CV") {
            add_warning("PIPES", id, record.line_number,
                        "Unknown initial status '" + record.fields[7] +
                        "'; pipe was imported open.");
            enabled = true;
        }
        const auto status_override = status_overrides.find(id);
        if (status_override != status_overrides.end()) {
            enabled = status_override->second == "OPEN";
            handled_status_ids.insert(id);
        }
        auto pipe = std::make_unique<Cso>(id, from, to, density, length, diameter,
                                          roughness, 0.0, 0.0, 1.0);
        pipe->Set_dprop("minor_loss", minor_loss);
        pipe->SetCheckValve(check_valve);
        pipe->Set_enabled(enabled);
        edges.push_back(std::move(pipe));
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
        EpanetPumpMetadata metadata = global_pump_energy;
        metadata.base_speed = 1.0;
        const auto specific_energy = pump_energy.find(id);
        if (specific_energy != pump_energy.end()) {
            const EpanetPumpMetadata &specific = specific_energy->second;
            if (specific.energy_price_specified) {
                metadata.energy_price_specified = true;
                metadata.energy_price_global = false;
                metadata.energy_price = specific.energy_price;
            }
            if (!specific.energy_pattern_id.empty()) {
                metadata.energy_pattern_id = specific.energy_pattern_id;
                metadata.energy_pattern_global = false;
            }
            if (!specific.efficiency_curve_id.empty()) {
                metadata.efficiency_curve_id = specific.efficiency_curve_id;
                metadata.efficiency_curve_global = false;
            }
        }
        std::string definition_value;
        for (std::size_t field = 3; field + 1 < record.fields.size(); field += 2) {
            const std::string keyword = upper(record.fields[field]);
            if (keyword == "POWER" || keyword == "HEAD") {
                if (!metadata.definition.empty())
                    add_warning("PUMPS", id, record.line_number,
                                "Multiple POWER/HEAD definitions; the last one was used.");
                metadata.definition = keyword;
                definition_value = record.fields[field + 1];
            } else if (keyword == "SPEED") {
                metadata.base_speed = number(record, field + 1, "pump speed", 1.0);
                if (metadata.base_speed < 0.0) {
                    add_warning("PUMPS", id, record.line_number,
                                "Negative pump speed was replaced with zero.");
                    metadata.base_speed = 0.0;
                }
            } else if (keyword == "PATTERN") {
                metadata.speed_pattern_id = record.fields[field + 1];
                const auto pattern = patterns_.find(metadata.speed_pattern_id);
                if (pattern != patterns_.end())
                    metadata.speed_pattern_values = pattern->second;
                else
                    add_warning("PUMPS", id, record.line_number,
                                "Referenced pump speed pattern was not found; base speed was used.");
            } else {
                add_warning("PUMPS", id, record.line_number,
                            "Unknown pump keyword '" + record.fields[field] + "' was retained only in the INP document.");
            }
        }
        if ((record.fields.size() - 3) % 2 != 0)
            add_warning("PUMPS", id, record.line_number,
                        "Pump keyword without a value was ignored.");
        const std::string type = metadata.definition;
        std::unique_ptr<Agelem> pump;
        if (type == "POWER") {
            Record value_record = record;
            value_record.fields = {id, definition_value};
            const double power = pump_power_to_watts(number(value_record, 1, "power", 0.0));
            pump = std::make_unique<EpanetPowerPump>(id, from, to, density, power, 1.0);
        } else if (type == "HEAD") {
            metadata.head_curve_id = definition_value;
            const auto curve = curves_.find(metadata.head_curve_id);
            if (curve == curves_.end() || curve->second.empty()) {
                add_warning("PUMPS", id, record.line_number,
                            "Pump curve is missing or empty; pump was not imported.");
                continue;
            }
            std::vector<double> flow;
            std::vector<double> head;
            for (const auto &point : curve->second) {
                flow.push_back(flow_to_m3_per_hour(point.first));
                head.push_back(length_to_metres(point.second));
                metadata.head_curve_points.push_back(std::make_pair(
                    flow_to_m3_per_hour(point.first) / 3600.0,
                    length_to_metres(point.second)));
            }
            if (flow.size() == 1) {
                const double design_flow = flow.front();
                const double design_head = head.front();
                flow = {0.0, design_flow, 2.0 * design_flow};
                head = {4.0 * design_head / 3.0, design_head, 0.0};
            } else if (flow.size() == 2) {
                flow.push_back(2.0 * flow.back());
                head.push_back(std::max(0.0, 2.0 * head.back() - head.front()));
            }
            pump = std::make_unique<Szivattyu>(id, from, to, density, 1.0, flow, head, 1.0);
        } else {
            add_warning("PUMPS", id, record.line_number,
                        "Unsupported pump parameter syntax; pump was not imported.");
        }
        if (pump != nullptr) {
            if (!metadata.energy_pattern_id.empty()) {
                const auto pattern = patterns_.find(metadata.energy_pattern_id);
                if (pattern != patterns_.end())
                    metadata.energy_pattern_values = pattern->second;
            }
            if (!metadata.efficiency_curve_id.empty()) {
                const auto curve = curves_.find(metadata.efficiency_curve_id);
                if (curve != curves_.end())
                    for (const auto &point : curve->second)
                        metadata.efficiency_curve_points.push_back(std::make_pair(
                            flow_to_m3_per_hour(point.first) / 3600.0, point.second));
            }
            const auto setting_override = setting_overrides.find(id);
            if (setting_override != setting_overrides.end()) {
                metadata.initial_setting_specified = true;
                metadata.initial_setting = setting_override->second;
                pump->Set_enabled(setting_override->second > 0.0);
                handled_status_ids.insert(id);
            }
            auto *configuration = dynamic_cast<EpanetPumpConfigurable *>(pump.get());
            configuration->SetEpanetPumpMetadata(metadata);
            if (!metadata.speed_pattern_values.empty())
                configuration->SetOperatingSpeed(metadata.speed_pattern_values.front());
            else if (metadata.initial_setting_specified)
                configuration->SetOperatingSpeed(metadata.initial_setting);
            const auto status_override = status_overrides.find(id);
            if (status_override != status_overrides.end()) {
                pump->Set_enabled(status_override->second == "OPEN");
                handled_status_ids.insert(id);
            }
            edges.push_back(std::move(pump));
        }
    }

    for (const Record &record : records("VALVES")) {
        const std::string id = record.fields.empty() ? "" : record.fields[0];
        if (!has_fields(record, 6, "VALVES"))
            continue;
        const std::string &from = record.fields[1];
        const std::string &to = record.fields[2];
        const std::string valve_type = upper(record.fields[4]);
        if (valve_type != "TCV") {
            add_warning("VALVES", id, record.line_number,
                        "EPANET valve type '" + record.fields[4] +
                        "' has no equivalent STACI element; valve was not imported.");
            continue;
        }
        if (node_ids.count(from) == 0 || node_ids.count(to) == 0) {
            add_warning("VALVES", id, record.line_number,
                        "Unknown endpoint; TCV was not imported.");
            continue;
        }
        if (!edge_ids.insert(id).second) {
            add_warning("VALVES", id, record.line_number,
                        "Duplicate link ID; TCV was not imported.");
            continue;
        }

        const double diameter = pipe_diameter_to_metres(
            number(record, 3, "valve diameter", 0.1));
        if (diameter <= 0.0) {
            add_warning("VALVES", id, record.line_number,
                        "Non-positive TCV diameter; valve was not imported.");
            continue;
        }
        double setting = number(record, 5, "TCV loss coefficient", 0.0);
        double minor_loss = record.fields.size() > 6
            ? number(record, 6, "valve minor loss", 0.0) : 0.0;
        if (setting < 0.0 || minor_loss < 0.0) {
            add_warning("VALVES", id, record.line_number,
                        "Negative TCV loss coefficients were replaced with zero.");
            setting = std::max(0.0, setting);
            minor_loss = std::max(0.0, minor_loss);
        }

        const auto setting_override = setting_overrides.find(id);
        if (setting_override != setting_overrides.end()) {
            setting = setting_override->second;
            handled_status_ids.insert(id);
        }
        EpanetTcvStatus status = EpanetTcvStatus::Active;
        const auto status_override = status_overrides.find(id);
        if (status_override != status_overrides.end()) {
            if (status_override->second == "OPEN")
                status = EpanetTcvStatus::Open;
            else if (status_override->second == "CLOSED")
                status = EpanetTcvStatus::Closed;
            handled_status_ids.insert(id);
        }

        const double area = 3.14159265358979323846 * diameter * diameter / 4.0;
        auto valve = std::make_unique<JelleggorbesFojtas>(
            id, from, to, density, area,
            std::vector<double>{0.0, 100.0},
            std::vector<double>{0.0, 0.0}, 50.0, 1.0);
        valve->SetEpanetTcvMetadata(setting, minor_loss, status);
        edges.push_back(std::move(valve));
    }

    for (const auto &status : status_overrides)
        if (handled_status_ids.count(status.first) == 0)
            add_warning("STATUS", status.first, 0,
                        "The referenced imported link was not found or does not support status.");
    for (const auto &setting : setting_overrides)
        if (handled_status_ids.count(setting.first) == 0)
            add_warning("STATUS", setting.first, 0,
                        "Numeric settings are currently supported only for imported pumps and TCVs.");

    if (!extended_period_)
        for (const Record &record : records("CONTROLS"))
            add_warning("CONTROLS", record.fields.size() > 1 ? record.fields[1] : "",
                        record.line_number,
                        "Dynamic link control was parsed but is not executed by steady-state STACI.");
    if (!extended_period_ && has_records("RULES"))
        add_warning("RULES", "", records("RULES").front().line_number,
                    "Rule-based controls are retained but require EPANET EPS mode for execution.");

    std::size_t demand_component_count = 0;
    std::size_t initial_quality_count = 0;
    for (const std::unique_ptr<Csomopont> &node : nodes) {
        demand_component_count += node->GetEpanetDemandComponents().size();
        if (node->GetEpanetInitialQuality().specified)
            ++initial_quality_count;
    }
    warn_for_unsupported_sections();
    print_report(nodes.size(), edges.size(), demand_component_count,
                 initial_quality_count, reservoir_head_pattern_count);
}

void EpanetReader::warn_for_unsupported_sections() {
    const std::vector<std::pair<std::string, std::string> > unsupported = {
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
                    "All multipliers are retained; only the first is applied to the "
                    "imported steady-state snapshot.");

    if (has_records("QUALITY") && quality_mode_ != "AGE" && quality_mode_ != "CHEMICAL")
        add_warning("QUALITY", "", records("QUALITY").front().line_number,
                    "Initial quality values are retained as metadata for quality mode '" +
                    quality_mode_ + "', but are not simulated.");

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
                                std::size_t edge_count,
                                std::size_t demand_component_count,
                                std::size_t initial_quality_count,
                                std::size_t reservoir_head_pattern_count) const {
    std::cerr << "\nEPANET import: " << filename_ << "\n"
              << "  units: " << flow_units_ << "\n"
              << "  imported STACI nodes: " << node_count << "\n"
              << "  imported STACI elements: " << edge_count << "\n"
              << "  preserved demand components: " << demand_component_count << "\n"
              << "  preserved initial quality values: " << initial_quality_count << "\n"
              << "  preserved reservoir head patterns: "
              << reservoir_head_pattern_count << "\n"
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
