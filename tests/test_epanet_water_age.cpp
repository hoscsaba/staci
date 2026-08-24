#include "epanet_water_age.h"

#include <cmath>
#include <iostream>
#include <vector>

namespace {
bool close_to(double actual, double expected, double tolerance,
              const char *label) {
    if (std::abs(actual - expected) <= tolerance)
        return true;
    std::cerr << label << ": expected " << expected << ", got " << actual << '\n';
    return false;
}
}

int main() {
    bool ok = true;

    EpanetWaterAgeModel forward(
        {{true}, {false}}, {{0, 1, 10.0}}, 5.0);
    forward.advance(20.0, {1.0});
    ok = close_to(forward.node_age_s()[0], 0.0, 1.0e-12,
                  "forward reservoir age") && ok;
    ok = close_to(forward.node_age_s()[1], 5.0, 1.0e-9,
                  "forward pipe travel time") && ok;

    EpanetWaterAgeModel reverse(
        {{false}, {true}}, {{0, 1, 10.0}}, 5.0);
    reverse.advance(20.0, {-1.0});
    ok = close_to(reverse.node_age_s()[1], 0.0, 1.0e-12,
                  "reverse reservoir age") && ok;
    ok = close_to(reverse.node_age_s()[0], 5.0, 1.0e-9,
                  "reverse pipe travel time") && ok;

    EpanetWaterAgeModel stagnant(
        {{true}, {false}}, {{0, 1, 10.0}}, 5.0);
    stagnant.advance(100.0, {0.0});
    ok = close_to(stagnant.node_age_s()[1], 100.0, 1.0e-9,
                  "stagnant node age") && ok;
    ok = close_to(stagnant.link_average_age_s({0.0})[0], 100.0, 1.0e-9,
                  "stagnant pipe age") && ok;

    EpanetChemicalModel chlorine(
        {{0.001, true}, {0.0, false}},
        {{0, 1, 10.0, -0.001, 0.0}}, 1.0);
    chlorine.advance(20.0, {1.0}, {0.0, 0.0}, {{}, {}});
    ok = close_to(chlorine.node_concentration_kgm3()[0], 0.001, 1.0e-12,
                  "fixed chlorine boundary") && ok;
    ok = close_to(chlorine.node_concentration_kgm3()[1],
                  0.001 * std::exp(-0.010), 3.0e-6,
                  "advective chlorine decay") && ok;

    EpanetChemicalModel reverse_chlorine(
        {{0.0, false}, {0.002, true}},
        {{0, 1, 5.0, 0.0, 0.0}}, 1.0);
    reverse_chlorine.advance(10.0, {-1.0}, {0.0, 0.0}, {{}, {}});
    ok = close_to(reverse_chlorine.node_concentration_kgm3()[0], 0.002, 1.0e-9,
                  "reverse chlorine transport") && ok;

    EpanetChemicalModel reservoir_source(
        {{0.0, true}, {0.0, false}},
        {{0, 1, 0.0, 0.0, 0.0}}, 1.0);
    reservoir_source.advance(
        1.0, {1.0}, {0.0, 0.0},
        {{EpanetChemicalSourceType::Concentration, 0.003}, {}});
    ok = close_to(reservoir_source.node_concentration_kgm3()[0], 0.003, 1.0e-12,
                  "reservoir concentration source") && ok;
    ok = close_to(reservoir_source.node_concentration_kgm3()[1], 0.003, 1.0e-12,
                  "zero-volume chemical link") && ok;

    return ok ? 0 : 1;
}
