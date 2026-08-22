#include "Cso.h"
#include "Csomopont.h"
#include "EpanetPowerPump.h"
#include "EpanetPump.h"
#include "KonstNyomas.h"
#include "Szivattyu.h"
#include "epanet_writer.h"

#include <fstream>
#include <cmath>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

namespace {

bool contains(const std::string &text, const std::string &expected) {
    if (text.find(expected) != std::string::npos)
        return true;
    std::cerr << "Missing exported text: " << expected << '\n';
    return false;
}

} // namespace

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: test_epanet_node_metadata output.inp\n";
        return 2;
    }

    Csomopont node("J1", 10.0, 0.0, 0.0, 0.0, 1000.0, 0.0);
    const std::vector<EpanetDemandComponent> demands = {
        {0.001, "P-BASE", {0.5, 1.0}, "", true},
        {0.002, "P-DOM", {1.5, 0.75}, "Domestic users", false},
        {0.003, "P-IND", {2.0, 0.25}, "Industrial process", false},
    };
    node.SetEpanetDemandComponents(demands, 1.2);

    EpanetInitialQuality quality;
    quality.specified = true;
    quality.source_value = 0.0;
    quality.mode = "CHEMICAL";
    quality.chemical_name = "Chlorine";
    quality.units = "mg/L";
    node.SetEpanetInitialQuality(quality);

    Csomopont reservoir("R1", 0.0, 0.0, 0.0, 30.0, 1000.0, 0.0);
    KonstNyomas reservoir_boundary(
        "EPANET_RESERVOIR_R1", 1.0, "R1", 1000.0,
        30.0 * 1000.0 * 9.81, 1.0, 0.0);
    reservoir_boundary.SetEpanetHeadPattern(
        EpanetHeadPattern{true, 30.0, "P-HEAD", {1.0, 1.1}});

    Cso check_valve_pipe("P-CV", "R1", "J1", 1000.0, 100.0, 0.15,
                         120.0, 0.0, 0.0, 1.0);
    check_valve_pipe.Set_friction_model("HW");
    check_valve_pipe.Set_dprop("minor_loss", 2.5);
    check_valve_pipe.SetCheckValve(true);
    Cso closed_pipe("P-CLOSED", "R1", "J1", 1000.0, 50.0, 0.10,
                    120.0, 0.0, 0.0, 1.0);
    closed_pipe.Set_friction_model("HW");
    closed_pipe.Set_dprop("minor_loss", 1.25);
    closed_pipe.Set_enabled(false);

    Cso loss_only("LOSS", "R1", "J1", 1000.0, 0.0, 0.10,
                  -0.02, 0.0, 0.0, 1.0);
    loss_only.Set_dprop("minor_loss", 2.0);
    const double area = 3.14159265358979323846 * 0.10 * 0.10 / 4.0;
    const double velocity = 1.0 / 1000.0 / area;
    const double expected_minor_headloss = 2.0 * velocity * velocity / (2.0 * 9.81);
    const double actual_minor_headloss =
        loss_only.ComputeHeadloss() / (1000.0 * 9.81);
    if (std::fabs(actual_minor_headloss - expected_minor_headloss) > 1.0e-12) {
        std::cerr << "Pipe minor-loss equation does not match K*v^2/(2g)\n";
        return 1;
    }
    closed_pipe.Ini(0, 0.0);
    if (closed_pipe.Get_mp() != 0.0 ||
        closed_pipe.f(std::vector<double>{20.0, 10.0, 0.0, 0.0}) != 0.0 ||
        closed_pipe.df(std::vector<double>{20.0, 10.0, 0.0, 0.0}) !=
            std::vector<double>({0.0, 0.0, 1.0, 0.0})) {
        std::cerr << "Closed pipe does not enforce zero flow\n";
        return 1;
    }
    Cso directional_pipe("DIRECTIONAL", "R1", "J1", 1000.0, 0.0, 0.10,
                         -0.02, 0.0, 0.0, -1.0);
    directional_pipe.SetCheckValve(true);
    if (directional_pipe.f(std::vector<double>{10.0, 20.0, 0.0, 0.0}) != -1.0) {
        std::cerr << "CV pipe did not close against reverse flow and head\n";
        return 1;
    }
    directional_pipe.Set_mp(0.0);
    if (std::fabs(directional_pipe.f(
            std::vector<double>{20.0, 10.0, 0.0, 0.0}) + 10.0) > 1.0e-12) {
        std::cerr << "CV pipe did not open for forward driving head\n";
        return 1;
    }

    EpanetPowerPump power_pump("P-POWER", "R1", "J1", 1000.0, 10000.0, 1.0);
    EpanetPumpMetadata power_metadata;
    power_metadata.definition = "POWER";
    power_metadata.base_speed = 0.5;
    power_metadata.initial_setting_specified = true;
    power_metadata.initial_setting = 0.7;
    power_metadata.speed_pattern_id = "P-SPEED";
    power_metadata.speed_pattern_values = {0.5, 0.8};
    power_metadata.energy_price_specified = true;
    power_metadata.energy_price = 0.22;
    power_metadata.energy_pattern_id = "P-ENERGY";
    power_metadata.energy_pattern_values = {1.0, 1.5};
    power_metadata.efficiency_curve_id = "EFF-POWER";
    power_metadata.efficiency_curve_points = {
        {0.0, 50.0}, {0.01, 80.0}, {0.02, 60.0}};
    power_pump.SetEpanetPumpMetadata(power_metadata);
    if (std::fabs(power_pump.Get_dprop("effective_power") - 1250.0) > 1.0e-12) {
        std::cerr << "POWER pump does not follow the cubic speed affinity law\n";
        return 1;
    }

    const std::vector<double> pump_flow = {0.0, 36.0, 72.0};
    const std::vector<double> pump_head = {40.0, 30.0, 0.0};
    Szivattyu head_pump("P-HEAD-PUMP", "R1", "J1", 1000.0, 1.0,
                        pump_flow, pump_head, 1.0);
    const double nominal_head = head_pump.Get_PumpHeadAt(0.01);
    EpanetPumpMetadata head_metadata;
    head_metadata.definition = "HEAD";
    head_metadata.head_curve_id = "ORIGINAL-HEAD";
    head_metadata.head_curve_points = {
        {0.0, 40.0}, {0.01, 30.0}, {0.02, 0.0}};
    head_metadata.base_speed = 0.5;
    head_pump.SetEpanetPumpMetadata(head_metadata);
    if (std::fabs(head_pump.Get_PumpHeadAt(0.005) - 0.25 * nominal_head) > 1.0e-10) {
        std::cerr << "HEAD pump does not follow H_s(Q)=s^2 H(Q/s)\n";
        return 1;
    }

    const auto &stored_demands = node.GetEpanetDemandComponents();
    if (stored_demands.size() != 3 || stored_demands[0].pattern_id != "P-BASE" ||
        stored_demands[1].category != "Domestic users" ||
        stored_demands[2].pattern_values != std::vector<double>({2.0, 0.25}) ||
        node.GetEpanetDemandMultiplier() != 1.2 ||
        !node.GetEpanetInitialQuality().specified ||
        node.GetEpanetInitialQuality().source_value != 0.0) {
        std::cerr << "Csomopont did not retain structured EPANET metadata\n";
        return 1;
    }

    const EpanetHeadPattern &stored_head_pattern =
        reservoir_boundary.GetEpanetHeadPattern();
    if (!stored_head_pattern.specified || stored_head_pattern.base_head_m != 30.0 ||
        stored_head_pattern.pattern_id != "P-HEAD" ||
        stored_head_pattern.pattern_values != std::vector<double>({1.0, 1.1})) {
        std::cerr << "KonstNyomas did not retain the EPANET reservoir head pattern\n";
        return 1;
    }

    std::vector<Csomopont *> nodes = {&node, &reservoir};
    std::vector<Agelem *> edges = {
        &reservoir_boundary, &check_valve_pipe, &closed_pipe,
        &power_pump, &head_pump};
    EpanetWriter::write(argv[1], nodes, edges, "HW");

    std::ifstream input(argv[1], std::ios::binary);
    const std::string exported((std::istreambuf_iterator<char>(input)),
                               std::istreambuf_iterator<char>());
    const bool valid =
        contains(exported, "J1\t10\t1\tP-BASE") &&
        contains(exported, "J1\t2\tP-DOM\tDomestic users") &&
        contains(exported, "J1\t3\tP-IND\tIndustrial process") &&
        contains(exported, "R1\t30\tP-HEAD") &&
        contains(exported, "P-CV\tR1\tJ1\t100\t150\t120\t2.5\tCV") &&
        contains(exported, "P-CLOSED\tR1\tJ1\t50\t100\t120\t1.25\tClosed") &&
        contains(exported, "P-POWER\tR1\tJ1\tPOWER\t10\tSPEED\t0.5\tPATTERN\tP-SPEED") &&
        contains(exported, "P-HEAD-PUMP\tR1\tJ1\tHEAD\tORIGINAL-HEAD\tSPEED\t0.5") &&
        contains(exported, "ORIGINAL-HEAD\t10\t30") &&
        contains(exported, "EFF-POWER\t10\t80") &&
        contains(exported, "P-SPEED\t0.5\t0.8") &&
        contains(exported, "P-ENERGY\t1\t1.5") &&
        contains(exported, "P-POWER\t0.7") &&
        contains(exported, "PUMP\tP-POWER\tPRICE\t0.22") &&
        contains(exported, "PUMP\tP-POWER\tPATTERN\tP-ENERGY") &&
        contains(exported, "PUMP\tP-POWER\tEFFIC\tEFF-POWER") &&
        contains(exported, "P-BASE\t0.5\t1") &&
        contains(exported, "P-DOM\t1.5\t0.75") &&
        contains(exported, "P-IND\t2\t0.25") &&
        contains(exported, "P-HEAD\t1\t1.1") &&
        contains(exported, "[QUALITY]\n;Node\tInitQual\nJ1\t0") &&
        contains(exported, "QUALITY\tCHEMICAL\tChlorine\tmg/L") &&
        contains(exported, "DEMAND MULTIPLIER\t1.2");
    return valid ? 0 : 1;
}
