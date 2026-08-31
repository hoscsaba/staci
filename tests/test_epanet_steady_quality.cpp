#include "epanet_steady_quality.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void require_close(double actual, double expected, double relative_tolerance,
                   const char *label) {
    const double scale = std::max(1.0e-12, std::abs(expected));
    if (std::abs(actual - expected) > relative_tolerance * scale)
        throw std::runtime_error(std::string(label) + " mismatch: " +
                                 std::to_string(actual) + " vs " +
                                 std::to_string(expected));
}

EpanetSteadyQualityModel make_model(double first_volume) {
    std::vector<EpanetSteadyQualityNode> nodes(3);
    nodes[0].id = "R";
    nodes[0].fixed_age = true;
    nodes[0].fixed_chemical = true;
    nodes[0].fixed_concentration_kgm3 = 1.0e-3;
    nodes[1].id = "J1";
    nodes[2].id = "J2";

    const double reaction = -1.0e-4;
    std::vector<EpanetSteadyQualityLink> links{
        {"P1", 0, 1, 0.01, first_volume, reaction},
        {"P2", 1, 2, 0.01, 10.0, reaction},
        // Negative signed flow verifies automatic direction reversal: R -> J2.
        {"P3", 2, 0, -0.01, 5.0, reaction}};
    return EpanetSteadyQualityModel(std::move(nodes), std::move(links));
}

} // namespace

int main() {
    try {
        EpanetSteadyQualitySensitivityInput sensitivity;
        sensitivity.parameter_element = "P1";
        sensitivity.parameter_property = "volume";
        sensitivity.link_volume_derivative_m3 = {1.0, 0.0, 0.0};
        const EpanetSteadyQualityResult result =
            make_model(10.0).solve(true, true, &sensitivity);

        const double beta1000 = std::exp(-0.1);
        const double beta500 = std::exp(-0.05);
        require_close(result.node_age_s[0], 0.0, 1.0e-12, "reservoir age");
        require_close(result.node_age_s[1], 1000.0, 1.0e-12, "J1 age");
        require_close(result.node_age_s[2], 1250.0, 1.0e-12, "mixed J2 age");
        require_close(result.node_concentration_kgm3[1], 1.0e-3 * beta1000,
                      1.0e-12, "J1 concentration");
        require_close(result.node_concentration_kgm3[2],
                      0.5e-3 * (beta1000 * beta1000 + beta500),
                      1.0e-12, "mixed J2 concentration");

        const double delta = 1.0e-4;
        const EpanetSteadyQualityResult plus =
            make_model(10.0 + delta).solve(true, true);
        const EpanetSteadyQualityResult minus =
            make_model(10.0 - delta).solve(true, true);
        for (std::size_t index = 0; index < result.node_age_s.size(); ++index) {
            const double finite_difference =
                (plus.node_age_s[index] - minus.node_age_s[index]) /
                (2.0 * delta);
            require_close(result.node_age_sensitivity[index], finite_difference,
                          2.0e-7, "water-age sensitivity");
        }
        for (std::size_t index = 0;
             index < result.node_concentration_kgm3.size(); ++index) {
            const double finite_difference =
                (plus.node_concentration_kgm3[index] -
                 minus.node_concentration_kgm3[index]) / (2.0 * delta);
            require_close(result.node_concentration_sensitivity[index],
                          finite_difference, 2.0e-7,
                          "chemical sensitivity");
        }
        std::cout << "steady water-quality base solution and analytic sensitivity: PASS\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "steady water-quality test: FAIL: " << error.what() << '\n';
        return 1;
    }
}
