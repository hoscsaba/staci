#include "StaciSession.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#ifndef STACI_SOURCE_DIR
#define STACI_SOURCE_DIR "."
#endif

int main(int argc, char *argv[]) {
    try {
        const std::string input = argc > 1
            ? argv[1]
            : STACI_SOURCE_DIR "/tests/anytown_1med.spr";
        StaciSession network(input);
        const StaciSessionInfo &info = network.info();
        std::cout << "Loaded SPR network: " << input << '\n'
                  << "Nodes: " << info.node_ids.size()
                  << ", links: " << info.link_ids.size() << '\n';

        const StaciHydraulicStatus status = network.solve_hydraulics();
        if (!status.converged)
            throw std::runtime_error("The hydraulic solution did not converge.");

        const StaciNodeResults nodes = network.node_results();
        const StaciLinkResults links = network.link_results();
        const auto pressure_range = std::minmax_element(
            nodes.pressure_head_m.begin(), nodes.pressure_head_m.end());
        std::cout << std::fixed << std::setprecision(6)
                  << "Pressure-head range [m]: " << *pressure_range.first
                  << " ... " << *pressure_range.second << '\n';

        std::cout << "\nFirst five node results (SI):\n"
                  << "id, elevation_m, pressure_head_m, demand_m3s\n";
        for (std::size_t index = 0;
             index < std::min<std::size_t>(5, nodes.ids.size()); ++index)
            std::cout << nodes.ids[index] << ", " << nodes.elevation_m[index]
                      << ", " << nodes.pressure_head_m[index] << ", "
                      << nodes.demand_m3s[index] << '\n';

        std::cout << "\nFirst five link results (SI):\n"
                  << "id, type, flow_rate_m3s, velocity_ms, headloss_m\n";
        for (std::size_t index = 0;
             index < std::min<std::size_t>(5, links.ids.size()); ++index)
            std::cout << links.ids[index] << ", " << links.types[index]
                      << ", " << links.flow_rate_m3s[index] << ", "
                      << links.velocity_ms[index] << ", "
                      << links.headloss_m[index] << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "STACI example failed: " << error.what() << '\n';
        return 1;
    }
}
