#include "StaciSession.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char **argv) {
    if (argc != 2) return 2;
    StaciSession model(argv[1]);
    const auto bulk = model.get_link_property({"P1"}, "bulk_reaction_per_s");
    const auto wall = model.get_link_property({"P1"}, "wall_reaction_m_s");
    const auto source = model.get_node_property({"R1"}, "source_concentration_kgm3");
    if (std::abs(bulk.front() + 0.5 / 86400.0) > 1e-14 ||
        std::abs(wall.front() + 0.001 / 86400.0) > 1e-14 ||
        std::abs(source.front() - 0.001) > 1e-12)
        throw std::runtime_error("INP chemical settings were not imported in SI units");
    if (!model.solve_hydraulics().converged)
        throw std::runtime_error("Hydraulics did not converge");
    const auto first = model.solve_steady_quality(true, true);
    const auto index = std::find(first.node_ids.begin(), first.node_ids.end(), "J1")
        - first.node_ids.begin();
    if (index >= first.node_ids.size() ||
        !(first.node_concentration_kgm3[index] > 0.0 &&
          first.node_concentration_kgm3[index] < 0.001) ||
        !(first.node_water_age_s[index] > 0.0))
        throw std::runtime_error("Imported source did not reach J1");

    model.set_link_property({"P1"}, "bulk_reaction_per_s", {-1.0 / 86400.0});
    model.set_node_property({"R1"}, "source_concentration_kgm3", {0.0008});
    if (!model.solve_hydraulics().converged)
        throw std::runtime_error("Hydraulics failed after changing quality data");
    const auto second = model.solve_steady_quality(false, true);
    if (!(second.node_concentration_kgm3[index] <
          first.node_concentration_kgm3[index]))
        throw std::runtime_error("In-memory chemical setters had no effect");
    std::cout << "INP chemistry and in-memory steady quality: passed\n";
}
