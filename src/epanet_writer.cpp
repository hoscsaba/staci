#include "epanet_writer.h"

#include "Agelem.h"
#include "Cso.h"
#include "Csomopont.h"
#include "EpanetPowerPump.h"
#include "EpanetPump.h"
#include "JelleggorbesFojtas.h"
#include "KonstNyomas.h"
#include "Szivattyu.h"
#include "epanet_document.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
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

std::vector<std::string> fields(const std::string &value) {
    std::istringstream input(value);
    std::vector<std::string> result;
    std::string field;
    while (input >> field)
        result.push_back(field);
    return result;
}

bool is_us_flow_unit(const std::string &unit) {
    const std::string normalized = upper(unit);
    return normalized == "CFS" || normalized == "GPM" || normalized == "MGD" ||
           normalized == "IMGD" || normalized == "AFD";
}

double flow_from_m3_per_hour(double value, const std::string &unit) {
    const std::string normalized = upper(unit);
    if (normalized == "CFS") return value / 101.9406477312;
    if (normalized == "GPM") return value / 0.22712470704;
    if (normalized == "MGD") return value / 157.725491;
    if (normalized == "IMGD") return value / 189.420142;
    if (normalized == "AFD") return value / 51.3950766;
    if (normalized == "LPS") return value / 3.6;
    if (normalized == "LPM") return value / 0.06;
    if (normalized == "MLD") return value / 41.6666666667;
    if (normalized == "CMH") return value;
    if (normalized == "CMD") return value * 24.0;
    return value / 3.6;
}

double length_from_metres(double value, bool us_units) {
    return us_units ? value / 0.3048 : value;
}

double pipe_diameter_from_metres(double value, bool us_units) {
    return us_units ? value / 0.0254 : value * 1000.0;
}

std::string safe_id(const std::string &source, std::set<std::string> &used) {
    std::string result = source;
    for (char &character : result) {
        const unsigned char c = static_cast<unsigned char>(character);
        if (std::isspace(c) || character == ';' || character == '[' || character == ']')
            character = '_';
    }
    if (result.empty())
        result = "STACI_ELEMENT";
    const std::string base = result;
    unsigned int suffix = 2;
    while (!used.insert(result).second)
        result = base + "_" + std::to_string(suffix++);
    if (result != source)
        std::cerr << "WARNING [EPANET][EXPORT][ID]: '" << source
                  << "' was exported as '" << result << "'.\n";
    return result;
}

std::string tank_node_id(Agelem *edge) {
    return edge->Get_Cspe_Nev() != "<nincs>" ? edge->Get_Cspe_Nev()
                                              : edge->Get_Cspv_Nev();
}

void replace_field(std::vector<std::string> &row, std::size_t index, double value) {
    while (row.size() <= index)
        row.push_back("0");
    std::ostringstream text;
    text << std::setprecision(15) << value;
    row[index] = text.str();
}

} // namespace

void EpanetWriter::write(const std::string &filename,
                         const std::vector<Csomopont *> &nodes,
                         const std::vector<Agelem *> &edges,
                         const std::string &friction_model) {
    std::ofstream output(filename.c_str(), std::ios::trunc);
    if (!output)
        throw std::runtime_error("Cannot create EPANET input file: " + filename);

    output << std::setprecision(15);
    std::set<std::string> used_node_ids;
    std::set<std::string> used_link_ids;
    std::map<std::string, std::string> node_ids;
    std::map<std::string, Csomopont *> nodes_by_name;
    for (Csomopont *node : nodes) {
        node_ids[node->Get_nev()] = safe_id(node->Get_nev(), used_node_ids);
        nodes_by_name[node->Get_nev()] = node;
    }

    std::set<std::string> reservoir_nodes;
    std::set<std::string> tank_nodes;
    for (Agelem *edge : edges) {
        if (edge->GetType() == "KonstNyomas")
            reservoir_nodes.insert(tank_node_id(edge));
        else if (edge->GetType() == "Vegakna")
            tank_nodes.insert(tank_node_id(edge));
    }

    output << "[TITLE]\nSTACI EPANET export\n\n[JUNCTIONS]\n"
           << ";ID\tElevation\tDemand\tPattern\n";
    for (Csomopont *node : nodes) {
        const std::string name = node->Get_nev();
        if (reservoir_nodes.count(name) || tank_nodes.count(name))
            continue;
        const std::vector<EpanetDemandComponent> &components =
            node->GetEpanetDemandComponents();
        const auto primary = std::find_if(
            components.begin(), components.end(),
            [](const EpanetDemandComponent &component) { return component.primary; });
        output << node_ids[name] << '\t' << node->Get_h() << '\t';
        if (primary == components.end()) {
            output << node->Get_dprop("demand") / 3.6;
        } else {
            output << primary->base_demand_m3s * 1000.0;
            if (!primary->pattern_id.empty())
                output << '\t' << primary->pattern_id;
        }
        output << '\n';
    }

    output << "\n[DEMANDS]\n;Junction\tDemand\tPattern\tCategory\n";
    for (Csomopont *node : nodes) {
        const std::string name = node->Get_nev();
        if (reservoir_nodes.count(name) || tank_nodes.count(name))
            continue;
        for (const EpanetDemandComponent &component : node->GetEpanetDemandComponents()) {
            if (component.primary)
                continue;
            output << node_ids[name] << '\t' << component.base_demand_m3s * 1000.0;
            if (!component.pattern_id.empty() || !component.category.empty())
                output << '\t' << component.pattern_id;
            if (!component.category.empty())
                output << '\t' << component.category;
            output << '\n';
        }
    }

    output << "\n[RESERVOIRS]\n;ID\tHead\n";
    for (Agelem *edge : edges) {
        if (edge->GetType() != "KonstNyomas")
            continue;
        const std::string node_name = tank_node_id(edge);
        if (node_ids.count(node_name) == 0) {
            std::cerr << "WARNING [EPANET][EXPORT][" << edge->Get_nev()
                      << "]: boundary node was not found; element skipped.\n";
            continue;
        }
        auto *boundary = dynamic_cast<KonstNyomas *>(edge);
        const EpanetHeadPattern *head_pattern = boundary == nullptr
            ? nullptr : &boundary->GetEpanetHeadPattern();
        output << node_ids[node_name] << '\t';
        if (head_pattern != nullptr && head_pattern->specified) {
            output << head_pattern->base_head_m;
            if (!head_pattern->pattern_id.empty())
                output << '\t' << head_pattern->pattern_id;
        } else {
            output << edge->Get_dprop("head") + nodes_by_name[node_name]->Get_h();
        }
        output << '\n';
    }

    output << "\n[TANKS]\n;ID\tElevation\tInitLevel\tMinLevel\tMaxLevel\tDiameter\tMinVolume\n";
    for (Agelem *edge : edges) {
        if (edge->GetType() != "Vegakna")
            continue;
        const std::string node_name = tank_node_id(edge);
        if (node_ids.count(node_name) == 0)
            continue;
        const double level = edge->Get_dprop("water_level");
        const double diameter = std::sqrt(4.0 * edge->Get_Aref() / 3.14159265358979323846);
        output << node_ids[node_name] << '\t' << edge->Get_dprop("bottom_level")
               << '\t' << level << "\t0\t" << std::max(level, 1.0)
               << '\t' << diameter << "\t0\n";
        std::cerr << "WARNING [EPANET][EXPORT][" << edge->Get_nev()
                  << "]: tank minimum and maximum levels were inferred.\n";
    }

    output << "\n[PIPES]\n;ID\tNode1\tNode2\tLength\tDiameter\tRoughness\tMinorLoss\tStatus\n";
    for (Agelem *edge : edges) {
        if (edge->GetType() != "Cso")
            continue;
        const std::string from = edge->Get_Cspe_Nev();
        const std::string to = edge->Get_Cspv_Nev();
        if (node_ids.count(from) == 0 || node_ids.count(to) == 0) {
            std::cerr << "WARNING [EPANET][EXPORT][" << edge->Get_nev()
                      << "]: pipe endpoint was not found; element skipped.\n";
            continue;
        }
        output << safe_id(edge->Get_nev(), used_link_ids) << '\t'
               << node_ids[from] << '\t' << node_ids[to] << '\t'
               << edge->Get_dprop("length") << '\t'
               << edge->Get_dprop("diameter") * 1000.0 << '\t'
               << edge->Get_dprop("friction_coeff") << '\t'
               << edge->Get_dprop("minor_loss") << '\t';
        auto *pipe = dynamic_cast<Cso *>(edge);
        if (pipe != nullptr && pipe->IsCheckValve())
            output << "CV\n";
        else
            output << (edge->Is_enabled() ? "Open\n" : "Closed\n");
    }

    struct ExportedPump {
        std::string id;
        Agelem *edge;
        EpanetPumpConfigurable *configuration;
        Szivattyu *curve_pump;
        std::string curve_id;
    };
    std::vector<ExportedPump> exported_pumps;
    output << "\n[PUMPS]\n;ID\tNode1\tNode2\tParameters\n";
    for (Agelem *edge : edges) {
        auto *power_pump = dynamic_cast<EpanetPowerPump *>(edge);
        auto *curve_pump = dynamic_cast<Szivattyu *>(edge);
        if (power_pump == nullptr && curve_pump == nullptr)
            continue;
        const std::string from = edge->Get_Cspe_Nev();
        const std::string to = edge->Get_Cspv_Nev();
        if (node_ids.count(from) == 0 || node_ids.count(to) == 0)
            continue;
        const std::string link_id = safe_id(edge->Get_nev(), used_link_ids);
        auto *configuration = dynamic_cast<EpanetPumpConfigurable *>(edge);
        const EpanetPumpMetadata *metadata = configuration == nullptr
            ? nullptr : &configuration->GetEpanetPumpMetadata();
        output << link_id << '\t' << node_ids[from] << '\t' << node_ids[to] << '\t';
        std::string curve_id;
        if (power_pump != nullptr) {
            output << "POWER\t" << power_pump->Get_dprop("power") / 1000.0;
        } else {
            curve_id = metadata != nullptr && !metadata->head_curve_id.empty()
                ? metadata->head_curve_id
                : "STACI_CURVE_" + std::to_string(exported_pumps.size() + 1);
            output << "HEAD\t" << curve_id;
        }
        if (metadata != nullptr && metadata->base_speed != 1.0)
            output << "\tSPEED\t" << metadata->base_speed;
        if (metadata != nullptr && !metadata->speed_pattern_id.empty())
            output << "\tPATTERN\t" << metadata->speed_pattern_id;
        output << '\n';
        exported_pumps.push_back(
            ExportedPump{link_id, edge, configuration, curve_pump, curve_id});
    }

    struct ExportedTcv {
        std::string id;
        JelleggorbesFojtas *valve;
    };
    std::vector<ExportedTcv> exported_tcvs;
    output << "\n[VALVES]\n;ID\tNode1\tNode2\tDiameter\tType\tSetting\tMinorLoss\n";
    for (Agelem *edge : edges) {
        auto *valve = dynamic_cast<JelleggorbesFojtas *>(edge);
        if (valve == nullptr || !valve->CanExportAsEpanetTcv())
            continue;
        const std::string from = edge->Get_Cspe_Nev();
        const std::string to = edge->Get_Cspv_Nev();
        if (node_ids.count(from) == 0 || node_ids.count(to) == 0) {
            std::cerr << "WARNING [EPANET][EXPORT][" << edge->Get_nev()
                      << "]: TCV endpoint was not found; element skipped.\n";
            continue;
        }
        const std::string link_id = safe_id(edge->Get_nev(), used_link_ids);
        const double diameter_mm =
            std::sqrt(4.0 * edge->Get_Aref() / 3.14159265358979323846) * 1000.0;
        output << link_id << '\t' << node_ids[from] << '\t' << node_ids[to]
               << '\t' << diameter_mm << "\tTCV\t"
               << valve->GetEpanetTcvSetting() << '\t'
               << valve->GetEpanetTcvMinorLoss() << '\n';
        exported_tcvs.push_back(ExportedTcv{link_id, valve});
    }

    output << "\n[CURVES]\n;ID\tFlow\tHead\n";
    std::set<std::string> written_curve_ids;
    for (const ExportedPump &pump : exported_pumps) {
        if (pump.curve_pump == nullptr)
            continue;
        const EpanetPumpMetadata *metadata = pump.configuration == nullptr
            ? nullptr : &pump.configuration->GetEpanetPumpMetadata();
        if (!written_curve_ids.insert(pump.curve_id).second)
            continue;
        std::vector<std::pair<double, double> > points;
        const bool retained_points =
            metadata != nullptr && !metadata->head_curve_points.empty();
        if (retained_points) {
            for (const auto &point : metadata->head_curve_points)
                points.push_back(std::make_pair(point.first * 1000.0, point.second));
        } else {
            const std::vector<double> flow = pump.curve_pump->GetCurveFlowM3PerHour();
            const std::vector<double> &head = pump.curve_pump->GetCurveHead();
            for (std::size_t i = 0; i < flow.size() && i < head.size(); ++i)
                points.push_back(std::make_pair(flow[i] / 3.6, head[i]));
        }
        if (!retained_points)
            std::sort(points.begin(), points.end());
        for (const auto &point : points)
            output << pump.curve_id << '\t' << point.first << '\t' << point.second << '\n';
    }
    for (const ExportedPump &pump : exported_pumps) {
        if (pump.configuration == nullptr)
            continue;
        const EpanetPumpMetadata &metadata = pump.configuration->GetEpanetPumpMetadata();
        if (metadata.efficiency_curve_id.empty() ||
            metadata.efficiency_curve_points.empty() ||
            !written_curve_ids.insert(metadata.efficiency_curve_id).second)
            continue;
        for (const auto &point : metadata.efficiency_curve_points)
            output << metadata.efficiency_curve_id << '\t'
                   << point.first * 1000.0 << '\t' << point.second << '\n';
    }

    std::map<std::string, std::vector<double> > demand_patterns;
    for (Csomopont *node : nodes) {
        for (const EpanetDemandComponent &component : node->GetEpanetDemandComponents()) {
            if (component.pattern_id.empty() || component.pattern_values.empty())
                continue;
            const auto inserted = demand_patterns.emplace(
                component.pattern_id, component.pattern_values);
            if (!inserted.second && inserted.first->second != component.pattern_values)
                std::cerr << "WARNING [EPANET][EXPORT][" << component.pattern_id
                          << "]: conflicting retained demand patterns; the first was used.\n";
        }
    }
    for (Agelem *edge : edges) {
        auto *boundary = dynamic_cast<KonstNyomas *>(edge);
        if (boundary == nullptr)
            continue;
        const EpanetHeadPattern &head_pattern = boundary->GetEpanetHeadPattern();
        if (head_pattern.pattern_id.empty() || head_pattern.pattern_values.empty())
            continue;
        const auto inserted = demand_patterns.emplace(
            head_pattern.pattern_id, head_pattern.pattern_values);
        if (!inserted.second && inserted.first->second != head_pattern.pattern_values)
            std::cerr << "WARNING [EPANET][EXPORT][" << head_pattern.pattern_id
                      << "]: conflicting retained patterns; the first was used.\n";
    }
    for (const ExportedPump &pump : exported_pumps) {
        if (pump.configuration == nullptr)
            continue;
        const EpanetPumpMetadata &metadata = pump.configuration->GetEpanetPumpMetadata();
        if (!metadata.speed_pattern_id.empty() && !metadata.speed_pattern_values.empty())
            demand_patterns.emplace(metadata.speed_pattern_id, metadata.speed_pattern_values);
        if (!metadata.energy_pattern_id.empty() && !metadata.energy_pattern_values.empty())
            demand_patterns.emplace(metadata.energy_pattern_id, metadata.energy_pattern_values);
    }
    output << "\n[PATTERNS]\n;ID\tMultipliers\n";
    for (const auto &pattern : demand_patterns) {
        output << pattern.first;
        for (double multiplier : pattern.second)
            output << '\t' << multiplier;
        output << '\n';
    }

    output << "\n[STATUS]\n;ID\tStatus/Setting\n";
    for (const ExportedPump &pump : exported_pumps) {
        const EpanetPumpMetadata *metadata = pump.configuration == nullptr
            ? nullptr : &pump.configuration->GetEpanetPumpMetadata();
        if (!pump.edge->Is_enabled())
            output << pump.id << "\tClosed\n";
        else if (metadata != nullptr && metadata->initial_setting_specified)
            output << pump.id << '\t' << metadata->initial_setting << '\n';
    }
    for (const ExportedTcv &tcv : exported_tcvs) {
        const EpanetTcvStatus status = tcv.valve->GetEpanetTcvStatus();
        if (status == EpanetTcvStatus::Closed)
            output << tcv.id << "\tClosed\n";
        else if (status == EpanetTcvStatus::Open)
            output << tcv.id << "\tOpen\n";
    }

    output << "\n[ENERGY]\n";
    bool global_price_written = false;
    bool global_pattern_written = false;
    bool global_efficiency_written = false;
    bool demand_charge_written = false;
    for (const ExportedPump &pump : exported_pumps) {
        if (pump.configuration == nullptr)
            continue;
        const EpanetPumpMetadata &metadata = pump.configuration->GetEpanetPumpMetadata();
        if (metadata.energy_price_specified) {
            if (metadata.energy_price_global && !global_price_written) {
                output << "GLOBAL\tPRICE\t" << metadata.energy_price << '\n';
                global_price_written = true;
            } else if (!metadata.energy_price_global) {
                output << "PUMP\t" << pump.id << "\tPRICE\t"
                       << metadata.energy_price << '\n';
            }
        }
        if (!metadata.energy_pattern_id.empty()) {
            if (metadata.energy_pattern_global && !global_pattern_written) {
                output << "GLOBAL\tPATTERN\t" << metadata.energy_pattern_id << '\n';
                global_pattern_written = true;
            } else if (!metadata.energy_pattern_global) {
                output << "PUMP\t" << pump.id << "\tPATTERN\t"
                       << metadata.energy_pattern_id << '\n';
            }
        }
        if (!metadata.efficiency_curve_id.empty()) {
            if (metadata.efficiency_curve_global && !global_efficiency_written) {
                output << "GLOBAL\tEFFIC\t" << metadata.efficiency_curve_id << '\n';
                global_efficiency_written = true;
            } else if (!metadata.efficiency_curve_global) {
                output << "PUMP\t" << pump.id << "\tEFFIC\t"
                       << metadata.efficiency_curve_id << '\n';
            }
        }
        if (metadata.demand_charge_specified && !demand_charge_written) {
            output << "DEMAND\tCHARGE\t" << metadata.demand_charge << '\n';
            demand_charge_written = true;
        }
    }

    output << "\n[QUALITY]\n;Node\tInitQual\n";
    for (Csomopont *node : nodes) {
        const EpanetInitialQuality &quality = node->GetEpanetInitialQuality();
        if (quality.specified) {
            output << node_ids[node->Get_nev()] << '\t' << quality.source_value << '\n';
        } else {
            const double concentration = node->Get_dprop("concentration");
            if (concentration != 0.0)
                output << node_ids[node->Get_nev()] << '\t' << concentration << '\n';
        }
    }

    for (Agelem *edge : edges) {
        const std::string_view type = edge->GetType();
        if (type == "JelleggorbesFojtas") {
            if (!static_cast<JelleggorbesFojtas *>(edge)->CanExportAsEpanetTcv())
                std::cerr << "WARNING [EPANET][EXPORT][" << edge->Get_nev()
                          << "]: STACI throttle valve was omitted from [VALVES], so the connection "
                          << edge->Get_Cspe_Nev() << " -> " << edge->Get_Cspv_Nev()
                          << " is absent from the exported INP. EPANET TCV stores one loss "
                          << "coefficient, but this STACI element stores a position-dependent loss "
                          << "curve (current position=" << edge->Get_dprop("position")
                          << ", current STACI loss parameter=" << edge->Get_dprop("veszt")
                          << "). Only a constant non-negative STACI loss curve can be exported "
                          << "losslessly as an EPANET TCV.\n";
        } else if (type != "Cso" && type != "Szivattyu" &&
                   type != "EpanetPowerPump" && type != "KonstNyomas" &&
                   type != "Vegakna") {
            std::cerr << "WARNING [EPANET][EXPORT][" << edge->Get_nev() << "]: STACI type '"
                      << type << "' has no lossless EPANET mapping and was skipped.\n";
        }
    }

    std::string headloss = upper(friction_model);
    if (headloss != "HW" && headloss != "DW")
        headloss = "DW";
    double demand_multiplier = 1.0;
    const EpanetInitialQuality *quality_metadata = nullptr;
    for (Csomopont *node : nodes) {
        if (!node->GetEpanetDemandComponents().empty())
            demand_multiplier = node->GetEpanetDemandMultiplier();
        if (quality_metadata == nullptr &&
            (node->GetEpanetInitialQuality().specified ||
             node->GetEpanetInitialQuality().mode != "NONE"))
            quality_metadata = &node->GetEpanetInitialQuality();
    }
    output << "\n[OPTIONS]\nUNITS\tLPS\nHEADLOSS\t"
           << (headloss == "HW" ? "H-W" : "D-W") << '\n';
    if (quality_metadata == nullptr) {
        output << "QUALITY\tCHEMICAL\tSTACI\tmg/L\n";
    } else if (quality_metadata->mode == "CHEMICAL") {
        output << "QUALITY\tCHEMICAL\t"
               << (quality_metadata->chemical_name.empty()
                       ? "STACI" : quality_metadata->chemical_name)
               << '\t'
               << (quality_metadata->units.empty() ? "mg/L" : quality_metadata->units)
               << '\n';
    } else if (quality_metadata->mode == "TRACE") {
        output << "QUALITY\tTRACE\t" << quality_metadata->trace_node << '\n';
    } else {
        output << "QUALITY\t" << quality_metadata->mode << '\n';
    }
    output << "DEMAND MULTIPLIER\t" << demand_multiplier
           << "\n\n[TIMES]\nDURATION\t0\n\n[END]\n";

    std::cerr << "EPANET export: wrote " << filename << " using metric LPS units.\n";
}

void EpanetWriter::copy(const std::string &input_filename,
                        const std::string &output_filename) {
    write_document(output_filename, EpanetDocument::read(input_filename));
}

void EpanetWriter::write_document(const std::string &output_filename,
                                  const EpanetDocument &document) {
    document.write(output_filename);
}

void EpanetWriter::write_modified_copy(const std::string &input_filename,
                                       const std::string &output_filename,
                                       const std::string &element_id,
                                       const std::string &property,
                                       double value_si) {
    std::ifstream input(input_filename.c_str());
    if (!input)
        throw std::runtime_error("Cannot open EPANET input file: " + input_filename);

    std::vector<std::string> lines;
    std::string line;
    std::string section;
    std::string flow_units = "LPS";
    std::string default_pattern;
    double demand_multiplier = 1.0;
    std::map<std::string, double> first_pattern_values;
    std::vector<std::vector<std::string> > demand_rows;
    while (std::getline(input, line)) {
        lines.push_back(line);
        const std::size_t comment = line.find(';');
        const std::string data = trim(line.substr(0, comment));
        if (!data.empty() && data.front() == '[' && data.back() == ']') {
            section = upper(trim(data.substr(1, data.size() - 2)));
            continue;
        }
        const std::vector<std::string> row = fields(data);
        if (section == "OPTIONS" && row.size() > 1) {
            if (upper(row[0]) == "UNITS")
                flow_units = upper(row[1]);
            else if (upper(row[0]) == "PATTERN")
                default_pattern = row[1];
            else if (row.size() > 2 && upper(row[0]) == "DEMAND" &&
                     upper(row[1]) == "MULTIPLIER")
                demand_multiplier = std::atof(row[2].c_str());
        } else if (section == "PATTERNS" && row.size() > 1 &&
                   first_pattern_values.count(row[0]) == 0) {
            first_pattern_values[row[0]] = std::atof(row[1].c_str());
        } else if (section == "DEMANDS" && row.size() > 1 && row[0] == element_id) {
            demand_rows.push_back(row);
        }
    }

    const bool us_units = is_us_flow_unit(flow_units);
    std::string target_section;
    std::size_t target_field = 0;
    double output_value = value_si;
    std::string target_id = element_id;
    if (property == "diameter") {
        target_section = "PIPES";
        target_field = 4;
        output_value = pipe_diameter_from_metres(value_si, us_units);
    } else if (property == "minor_loss") {
        target_section = "PIPES";
        target_field = 6;
        if (value_si < 0.0)
            throw std::runtime_error("Pipe minor-loss coefficient cannot be negative.");
    } else if (property == "tcv_setting") {
        target_section = "VALVES";
        target_field = 5;
        if (value_si < 0.0)
            throw std::runtime_error("TCV setting cannot be negative.");
    } else if (property == "tcv_minor_loss") {
        target_section = "VALVES";
        target_field = 6;
        if (value_si < 0.0)
            throw std::runtime_error("TCV minor-loss coefficient cannot be negative.");
    } else if (property == "demand") {
        target_section = "JUNCTIONS";
        target_field = 2;
        output_value = flow_from_m3_per_hour(value_si, flow_units);
    } else if (property == "bottom_level") {
        target_section = "TANKS";
        target_field = 1;
        output_value = length_from_metres(value_si, us_units);
        const std::string prefix = "EPANET_TANK_";
        if (target_id.compare(0, prefix.size(), prefix) == 0)
            target_id.erase(0, prefix.size());
    } else if (property == "water_level") {
        target_section = "TANKS";
        target_field = 2;
        output_value = length_from_metres(value_si, us_units);
        const std::string prefix = "EPANET_TANK_";
        if (target_id.compare(0, prefix.size(), prefix) == 0)
            target_id.erase(0, prefix.size());
    } else {
        throw std::runtime_error("Property '" + property +
                                 "' cannot be written back to an EPANET .inp file.");
    }

    bool replaced = false;
    section.clear();
    for (std::string &source_line : lines) {
        const bool had_carriage_return = !source_line.empty() && source_line.back() == '\r';
        const std::size_t comment_at = source_line.find(';');
        const std::string comment = comment_at == std::string::npos ? "" : source_line.substr(comment_at);
        const std::string data = trim(source_line.substr(0, comment_at));
        if (!data.empty() && data.front() == '[' && data.back() == ']') {
            section = upper(trim(data.substr(1, data.size() - 2)));
            continue;
        }
        std::vector<std::string> row = fields(data);
        if (section != target_section || row.empty() || row[0] != target_id)
            continue;
        if (property == "demand") {
            const std::string pattern_id = row.size() > 3 ? row[3] : default_pattern;
            const auto pattern = first_pattern_values.find(pattern_id);
            const double base_pattern = pattern_id.empty() || pattern == first_pattern_values.end()
                ? 1.0 : pattern->second;
            double additional_demand = 0.0;
            for (const std::vector<std::string> &demand_row : demand_rows) {
                const std::string additional_pattern_id = demand_row.size() > 2
                    ? demand_row[2] : default_pattern;
                const auto additional_pattern = first_pattern_values.find(additional_pattern_id);
                const double factor = additional_pattern_id.empty() ||
                    additional_pattern == first_pattern_values.end()
                    ? 1.0 : additional_pattern->second;
                additional_demand += std::atof(demand_row[1].c_str()) * factor;
            }
            if (std::fabs(demand_multiplier * base_pattern) < 1.0e-12)
                throw std::runtime_error("The requested demand cannot be represented because its "
                                         "time-zero pattern multiplier is zero.");
            output_value = (output_value / demand_multiplier - additional_demand) / base_pattern;
        }
        replace_field(row, target_field, output_value);
        std::ostringstream replacement;
        replacement << ' ';
        for (std::size_t i = 0; i < row.size(); ++i)
            replacement << (i == 0 ? "" : "\t") << row[i];
        if (!comment.empty())
            replacement << '\t' << comment;
        if (had_carriage_return && (comment.empty() || comment.back() != '\r'))
            replacement << '\r';
        source_line = replacement.str();
        replaced = true;
        break;
    }

    if (!replaced)
        throw std::runtime_error("Element '" + target_id + "' was not found in EPANET section [" +
                                 target_section + "].");

    std::ofstream output(output_filename.c_str(), std::ios::trunc);
    if (!output)
        throw std::runtime_error("Cannot create EPANET input file: " + output_filename);
    for (const std::string &output_line : lines)
        output << output_line << '\n';

    if (property == "demand")
        std::cerr << "WARNING [EPANET][MODIFY][" << element_id
                  << "]: the base demand was adjusted so the imported time-zero demand matches "
                     "the requested value; patterns and demand categories remain active.\n";
    std::cerr << "EPANET modification: wrote " << output_filename
              << " and preserved all other input sections.\n";
}
