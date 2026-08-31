#include "StaciSession.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

#ifndef STACI_SOURCE_DIR
#define STACI_SOURCE_DIR "."
#endif

int main(int argc, char *argv[]) {
    try {
        const std::string input = argc > 1
            ? argv[1]
            : STACI_SOURCE_DIR "/tests/epanet_eps_smoke.inp";
        StaciSession network(input);
        const StaciSessionInfo &info = network.info();
        std::map<std::string, std::size_t> type_counts;
        for (const std::string &type : info.link_types)
            ++type_counts[type];

        std::cout << "Imported EPANET INP: " << input << '\n'
                  << "STACI nodes: " << info.node_ids.size()
                  << ", STACI elements: " << info.link_ids.size() << '\n';
        for (const auto &entry : type_counts)
            std::cout << "  " << entry.first << ": " << entry.second << '\n';

        const StaciHydraulicStatus status = network.solve_hydraulics();
        if (!status.converged)
            throw std::runtime_error("The imported hydraulic snapshot did not converge.");
        const StaciNodeResults nodes = network.node_results();
        const StaciLinkResults links = network.link_results();

        std::cout << std::fixed << std::setprecision(6)
                  << "\nImported snapshot node results (SI):\n"
                  << "id, total_head_m, pressure_head_m, demand_m3s\n";
        for (std::size_t index = 0; index < nodes.ids.size(); ++index)
            std::cout << nodes.ids[index] << ", " << nodes.total_head_m[index]
                      << ", " << nodes.pressure_head_m[index] << ", "
                      << nodes.demand_m3s[index] << '\n';

        std::cout << "\nPhysical EPANET links represented by STACI:\n"
                  << "id, type, from, to, flow_rate_m3s\n";
        for (std::size_t index = 0; index < links.ids.size(); ++index) {
            if (links.to_nodes[index].empty())
                continue; // synthetic reservoir/tank boundary element
            std::cout << links.ids[index] << ", " << links.types[index]
                      << ", " << links.from_nodes[index] << ", "
                      << links.to_nodes[index] << ", "
                      << links.flow_rate_m3s[index] << '\n';
        }
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "EPANET import example failed: " << error.what() << '\n';
        return 1;
    }
}
