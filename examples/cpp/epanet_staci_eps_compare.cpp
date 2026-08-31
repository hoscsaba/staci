#include "Staci.h"
#include "epanet_extended_simulation.h"

#include "epanet2_2.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef STACI_SOURCE_DIR
#define STACI_SOURCE_DIR "."
#endif

namespace {

struct NodeRow {
    double head_m = 0.0;
    double pressure_m = 0.0;
    double demand_m3s = 0.0;
};

struct LinkRow {
    double flow_m3s = 0.0;
    double velocity_ms = 0.0;
};

struct ErrorSummary {
    double head_m = 0.0;
    double pressure_m = 0.0;
    double demand_m3s = 0.0;
    double flow_m3s = 0.0;
    double velocity_ms = 0.0;
    std::size_t node_values = 0;
    std::size_t link_values = 0;
};

std::vector<std::string> split_csv(const std::string &line) {
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for (std::size_t index = 0; index < line.size(); ++index) {
        const char value = line[index];
        if (value == '"') {
            if (quoted && index + 1 < line.size() && line[index + 1] == '"') {
                field.push_back('"');
                ++index;
            } else {
                quoted = !quoted;
            }
        } else if (value == ',' && !quoted) {
            fields.push_back(field);
            field.clear();
        } else {
            field.push_back(value);
        }
    }
    fields.push_back(field);
    return fields;
}

std::map<std::string, std::size_t> header_map(const std::string &line) {
    std::map<std::string, std::size_t> result;
    const std::vector<std::string> fields = split_csv(line);
    for (std::size_t index = 0; index < fields.size(); ++index)
        result[fields[index]] = index;
    return result;
}

std::string key(long time_s, const std::string &id) {
    return std::to_string(time_s) + "\x1f" + id;
}

std::map<std::string, NodeRow> read_staci_nodes(const std::string &filename) {
    std::ifstream input(filename);
    std::string line;
    if (!input || !std::getline(input, line))
        throw std::runtime_error("Cannot read STACI node CSV: " + filename);
    const auto columns = header_map(line);
    std::map<std::string, NodeRow> rows;
    while (std::getline(input, line)) {
        const auto values = split_csv(line);
        if (values.size() < columns.size())
            continue;
        const long time_s = std::stol(values.at(columns.at("time_seconds")));
        const std::string id = values.at(columns.at("node_id"));
        rows[key(time_s, id)] = NodeRow{
            std::stod(values.at(columns.at("total_head_m"))),
            std::stod(values.at(columns.at("pressure_head_m"))),
            std::stod(values.at(columns.at("demand_m3s")))};
    }
    return rows;
}

std::map<std::string, LinkRow> read_staci_links(const std::string &filename) {
    std::ifstream input(filename);
    std::string line;
    if (!input || !std::getline(input, line))
        throw std::runtime_error("Cannot read STACI link CSV: " + filename);
    const auto columns = header_map(line);
    std::map<std::string, LinkRow> rows;
    while (std::getline(input, line)) {
        const auto values = split_csv(line);
        if (values.size() < columns.size())
            continue;
        const long time_s = std::stol(values.at(columns.at("time_seconds")));
        const std::string id = values.at(columns.at("link_id"));
        rows[key(time_s, id)] = LinkRow{
            std::stod(values.at(columns.at("flow_m3s"))),
            std::stod(values.at(columns.at("velocity_mps")))};
    }
    return rows;
}

void epanet_check(int code, const char *operation) {
    if (code == 0)
        return;
    std::array<char, 256> message{};
    EN_geterror(code, message.data(), static_cast<int>(message.size()));
    throw std::runtime_error(std::string(operation) + ": " + message.data());
}

double flow_factor(int units) {
    static const std::array<double, 10> factors = {
        0.028316846592,                 // CFS
        0.003785411784 / 60.0,          // GPM
        1.0e6 * 0.003785411784 / 86400, // MGD
        1.0e6 * 0.00454609 / 86400,     // IMGD
        1233.48183754752 / 86400,        // AFD
        1.0e-3,                          // LPS
        1.0e-3 / 60.0,                   // LPM
        1000.0 / 86400,                  // MLD
        1.0 / 3600.0,                    // CMH
        1.0 / 86400.0                    // CMD
    };
    if (units < 0 || units >= static_cast<int>(factors.size()))
        throw std::runtime_error("Unknown EPANET flow unit code.");
    return factors[static_cast<std::size_t>(units)];
}

ErrorSummary compare_with_epanet(
    const std::string &input_file,
    const std::string &report_file,
    const std::map<std::string, NodeRow> &staci_nodes,
    const std::map<std::string, LinkRow> &staci_links) {
    EN_Project project = nullptr;
    epanet_check(EN_createproject(&project), "EN_createproject");
    try {
        epanet_check(EN_open(project, input_file.c_str(), report_file.c_str(), ""),
                     "EN_open");
        int units = 0;
        int node_count = 0;
        int link_count = 0;
        epanet_check(EN_getflowunits(project, &units), "EN_getflowunits");
        epanet_check(EN_getcount(project, EN_NODECOUNT, &node_count), "EN_getcount nodes");
        epanet_check(EN_getcount(project, EN_LINKCOUNT, &link_count), "EN_getcount links");
        const double q_to_si = flow_factor(units);
        const double length_to_si = units <= EN_AFD ? 0.3048 : 1.0;

        epanet_check(EN_openH(project), "EN_openH");
        epanet_check(EN_initH(project, EN_NOSAVE), "EN_initH");
        ErrorSummary errors;
        long step = 0;
        do {
            long time_s = 0;
            epanet_check(EN_runH(project, &time_s), "EN_runH");
            for (int index = 1; index <= node_count; ++index) {
                std::array<char, EN_MAXID + 1> id{};
                epanet_check(EN_getnodeid(project, index, id.data()), "EN_getnodeid");
                const auto found = staci_nodes.find(key(time_s, id.data()));
                if (found == staci_nodes.end())
                    continue;
                double head = 0.0, pressure = 0.0, demand = 0.0;
                int node_type = 0;
                epanet_check(EN_getnodetype(project, index, &node_type),
                             "EN_getnodetype");
                epanet_check(EN_getnodevalue(project, index, EN_HEAD, &head), "EN_HEAD");
                errors.head_m = std::max(errors.head_m,
                    std::abs(found->second.head_m - head * length_to_si));
                if (node_type == EN_JUNCTION) {
                    epanet_check(EN_getnodevalue(project, index, EN_PRESSURE,
                                                 &pressure), "EN_PRESSURE");
                    epanet_check(EN_getnodevalue(project, index, EN_DEMAND,
                                                 &demand), "EN_DEMAND");
                    errors.pressure_m = std::max(errors.pressure_m,
                        std::abs(found->second.pressure_m -
                                 pressure * length_to_si));
                    errors.demand_m3s = std::max(errors.demand_m3s,
                        std::abs(found->second.demand_m3s - demand * q_to_si));
                }
                ++errors.node_values;
            }
            for (int index = 1; index <= link_count; ++index) {
                std::array<char, EN_MAXID + 1> id{};
                epanet_check(EN_getlinkid(project, index, id.data()), "EN_getlinkid");
                const auto found = staci_links.find(key(time_s, id.data()));
                if (found == staci_links.end())
                    continue;
                double flow = 0.0, velocity = 0.0;
                epanet_check(EN_getlinkvalue(project, index, EN_FLOW, &flow), "EN_FLOW");
                epanet_check(EN_getlinkvalue(project, index, EN_VELOCITY, &velocity), "EN_VELOCITY");
                errors.flow_m3s = std::max(errors.flow_m3s,
                    std::abs(found->second.flow_m3s - flow * q_to_si));
                errors.velocity_ms = std::max(errors.velocity_ms,
                    std::abs(found->second.velocity_ms - velocity * length_to_si));
                ++errors.link_values;
            }
            epanet_check(EN_nextH(project, &step), "EN_nextH");
        } while (step > 0);
        epanet_check(EN_closeH(project), "EN_closeH");
        epanet_check(EN_close(project), "EN_close");
        EN_deleteproject(project);
        return errors;
    } catch (...) {
        if (project != nullptr) {
            EN_closeH(project);
            EN_close(project);
            EN_deleteproject(project);
        }
        throw;
    }
}

} // namespace

int main(int argc, char *argv[]) {
    try {
        const std::string input = argc > 1
            ? argv[1]
            : STACI_SOURCE_DIR "/tests/epanet_eps_smoke.inp";
        const std::filesystem::path output_directory = argc > 2
            ? std::filesystem::path(argv[2])
            : std::filesystem::current_path() / "epanet-staci-eps-results";
        std::filesystem::create_directories(output_directory);
        const std::string prefix = (output_directory / "staci").string();

        Staci staci(input, true);
        staci.Set_debug_level(0);
        staci.build_system();
        EpanetExtendedSimulation simulation(input, prefix);
        simulation.run(staci);

        const auto staci_nodes = read_staci_nodes(prefix + "-nodes.csv");
        const auto staci_links = read_staci_links(prefix + "-links.csv");
        const ErrorSummary errors = compare_with_epanet(
            input, (output_directory / "epanet.rpt").string(),
            staci_nodes, staci_links);
        if (errors.node_values == 0 || errors.link_values == 0)
            throw std::runtime_error(
                "No matching EPANET/STACI EPS result rows were found.");

        std::cout << std::setprecision(9)
                  << "\nEPANET 2.2 versus STACI EPS maximum absolute differences (SI)\n"
                  << "Compared node states: " << errors.node_values
                  << ", link states: " << errors.link_values << '\n'
                  << "  total head [m]       : " << errors.head_m << '\n'
                  << "  pressure head [m]    : " << errors.pressure_m << '\n'
                  << "  demand [m3/s]        : " << errors.demand_m3s << '\n'
                  << "  flow rate [m3/s]     : " << errors.flow_m3s << '\n'
                  << "  velocity [m/s]       : " << errors.velocity_ms << '\n'
                  << "STACI EPS files: " << output_directory << '\n';

        const bool passed = errors.head_m <= 0.003 &&
                            errors.pressure_m <= 0.003 &&
                            errors.demand_m3s <= 1.0e-6 &&
                            errors.flow_m3s <= 3.0e-4 &&
                            errors.velocity_ms <= 0.01;
        std::cout << "Benchmark status: " << (passed ? "PASS" : "FAIL") << '\n';
        return passed ? 0 : 2;
    } catch (const std::exception &error) {
        std::cerr << "EPS comparison example failed: " << error.what() << '\n';
        return 1;
    }
}
