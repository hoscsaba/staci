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

    // A 100 m3 tank at 1 mg/L receives 10 m3 at 0 mg/L and
    // releases 20 m3: mixed concentration is 100/110 mg/L.
    EpanetChemicalModel tank(
        {{0.0, true}, {0.001, false, 100.0, 0.0}, {0.0, false}},
        {{0, 1, 0.0, 0.0, 0.0}, {1, 2, 0.0, 0.0, 0.0}}, 10.0);
    tank.advance(10.0, {1.0, 2.0}, {0.0, 0.0, 0.0}, {{}, {}, {}});
    ok = close_to(tank.node_concentration_kgm3()[1], 0.1 / 110.0, 1e-12,
                  "stored tank dilution") && ok;
    ok = close_to(tank.node_concentration_kgm3()[2], 0.1 / 110.0, 1e-12,
                  "tank outlet mass") && ok;
    EpanetChemicalModel batch(
        {{0.001, false, 100.0, -0.01}}, {}, 10.0);
    batch.advance(10.0, {}, {0.0}, {{}});
    ok = close_to(batch.node_concentration_kgm3()[0], 0.0009, 1e-12,
                  "isolated tank reaction") && ok;
    EpanetChemicalModel terminal_source(
        {{0.001, true}, {0.0, false}}, {{0, 1, 0.0, 0.0, 0.0}}, 1.0);
    terminal_source.advance(1.0, {1.0}, {0.0, -1.0},
        {{}, {EpanetChemicalSourceType::Mass, 0.001}});
    ok = close_to(terminal_source.node_concentration_kgm3()[1], 0.002, 1e-12,
                  "mass source includes consumed flow") && ok;
    EpanetChemicalModel inert_stagnant(
        {{0.001, false}, {0.001, false}}, {{0, 1, 10.0, 0.0, 0.002}}, 1.0);
    inert_stagnant.advance(10.0, {0.0}, {0.0, 0.0}, {{}, {}});
    ok = close_to(inert_stagnant.node_concentration_kgm3()[0], 0.001, 1e-12,
                  "inert stagnant node retains initial value") && ok;
    return ok ? 0 : 1;
}
