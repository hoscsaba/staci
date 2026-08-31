#ifndef STACI_EPANET_STEADY_QUALITY_H
#define STACI_EPANET_STEADY_QUALITY_H

#include <cstddef>
#include <string>
#include <vector>

struct EpanetSteadyQualityNode {
    std::string id;
    bool fixed_age = false;
    bool fixed_chemical = false;
    double fixed_concentration_kgm3 = 0.0;
    double external_inflow_m3s = 0.0;
    double external_concentration_kgm3 = 0.0;
    double mass_source_kgs = 0.0;
};

struct EpanetSteadyQualityLink {
    std::string id;
    std::size_t from_node = 0;
    std::size_t to_node = 0;
    double flow_m3s = 0.0;
    double volume_m3 = 0.0;
    double reaction_coefficient_per_s = 0.0;
};

struct EpanetSteadyQualitySensitivityInput {
    std::string parameter_element;
    std::string parameter_property;
    std::vector<double> link_flow_derivative_m3s;
    std::vector<double> link_volume_derivative_m3;
    std::vector<double> link_reaction_derivative_per_s;
    std::vector<double> fixed_concentration_derivative_kgm3;
    std::vector<double> external_inflow_derivative_m3s;
    std::vector<double> external_concentration_derivative_kgm3;
    std::vector<double> mass_source_derivative_kgs;
};

struct EpanetSteadyQualityResult {
    std::vector<double> node_age_s;
    std::vector<double> link_average_age_s;
    std::vector<double> link_travel_time_s;
    std::vector<double> node_concentration_kgm3;
    std::vector<double> link_average_concentration_kgm3;
    std::vector<double> link_transfer_factor;

    std::vector<double> node_age_sensitivity;
    std::vector<double> link_average_age_sensitivity;
    std::vector<double> node_concentration_sensitivity;
    std::vector<double> link_average_concentration_sensitivity;
};

// Solves the asymptotic water-age and first-order chemical transport problem
// for a fixed hydraulic state. Pipes are directed by the signed flow supplied
// in EpanetSteadyQualityLink. Node mixing is instantaneous and complete.
class EpanetSteadyQualityModel {
public:
    EpanetSteadyQualityModel(std::vector<EpanetSteadyQualityNode> nodes,
                             std::vector<EpanetSteadyQualityLink> links);

    EpanetSteadyQualityResult solve(
        bool solve_age, bool solve_chemical,
        const EpanetSteadyQualitySensitivityInput *sensitivity = nullptr) const;

private:
    std::vector<EpanetSteadyQualityNode> nodes_;
    std::vector<EpanetSteadyQualityLink> links_;
};

#endif
