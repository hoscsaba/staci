#include "epanet_writer.h"

#include "Agelem.h"
#include "Csomopont.h"
#include "EpanetPowerPump.h"
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
           << ";ID\tElevation\tDemand\n";
    for (Csomopont *node : nodes) {
        const std::string name = node->Get_nev();
        if (reservoir_nodes.count(name) || tank_nodes.count(name))
            continue;
        output << node_ids[name] << '\t' << node->Get_h() << '\t'
               << node->Get_dprop("demand") / 3.6 << '\n';
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
        output << node_ids[node_name] << '\t'
               << edge->Get_dprop("head") + nodes_by_name[node_name]->Get_h() << '\n';
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
               << edge->Get_dprop("friction_coeff") << "\t0\tOpen\n";
    }

    std::vector<std::pair<std::string, Szivattyu *> > curve_pumps;
    output << "\n[PUMPS]\n;ID\tNode1\tNode2\tParameters\n";
    for (Agelem *edge : edges) {
        if (edge->GetType() != "Szivattyu" && edge->GetType() != "EpanetPowerPump")
            continue;
        const std::string from = edge->Get_Cspe_Nev();
        const std::string to = edge->Get_Cspv_Nev();
        if (node_ids.count(from) == 0 || node_ids.count(to) == 0)
            continue;
        const std::string link_id = safe_id(edge->Get_nev(), used_link_ids);
        output << link_id << '\t' << node_ids[from] << '\t' << node_ids[to] << '\t';
        if (edge->GetType() == "EpanetPowerPump") {
            output << "POWER\t" << edge->Get_dprop("power") / 1000.0 << '\n';
        } else {
            const std::string curve_id = "STACI_CURVE_" + std::to_string(curve_pumps.size() + 1);
            output << "HEAD\t" << curve_id << '\n';
            curve_pumps.push_back(std::make_pair(curve_id, dynamic_cast<Szivattyu *>(edge)));
        }
    }

    output << "\n[CURVES]\n;ID\tFlow\tHead\n";
    for (const auto &entry : curve_pumps) {
        const std::vector<double> flow = entry.second->GetCurveFlowM3PerHour();
        const std::vector<double> head = entry.second->GetCurveHead();
        std::vector<std::pair<double, double> > points;
        for (std::size_t i = 0; i < flow.size() && i < head.size(); ++i)
            points.push_back(std::make_pair(flow[i], head[i]));
        std::sort(points.begin(), points.end());
        for (const auto &point : points)
            output << entry.first << '\t' << point.first / 3.6 << '\t' << point.second << '\n';
    }

    output << "\n[QUALITY]\n;Node\tInitQual\n";
    for (Csomopont *node : nodes) {
        const double concentration = node->Get_dprop("concentration");
        if (concentration != 0.0)
            output << node_ids[node->Get_nev()] << '\t' << concentration << '\n';
    }

    for (Agelem *edge : edges) {
        const std::string type = edge->GetType();
        if (type != "Cso" && type != "Szivattyu" && type != "EpanetPowerPump" &&
            type != "KonstNyomas" && type != "Vegakna")
            std::cerr << "WARNING [EPANET][EXPORT][" << edge->Get_nev() << "]: STACI type '"
                      << type << "' has no lossless EPANET mapping and was skipped.\n";
    }

    std::string headloss = upper(friction_model);
    if (headloss != "HW" && headloss != "DW")
        headloss = "DW";
    output << "\n[OPTIONS]\nUNITS\tLPS\nHEADLOSS\t"
           << (headloss == "HW" ? "H-W" : "D-W")
           << "\nQUALITY\tCHEMICAL\tSTACI\tmg/L\n\n[TIMES]\nDURATION\t0\n\n[END]\n";

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
