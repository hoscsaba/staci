#ifndef STACI_SESSION_H
#define STACI_SESSION_H

#include "Staci.h"

#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

struct StaciSessionInfo {
    std::vector<std::string> node_ids;
    std::vector<std::string> link_ids;
    std::vector<std::string> link_types;
    std::vector<std::string> link_from_nodes;
    std::vector<std::string> link_to_nodes;
};

struct StaciHydraulicStatus {
    bool converged = false;
};

struct StaciNodeResults {
    std::vector<std::string> ids;
    std::vector<double> elevation_m;
    std::vector<double> pressure_head_m;
    std::vector<double> total_head_m;
    std::vector<double> demand_m3s;
};

struct StaciLinkResults {
    std::vector<std::string> ids;
    std::vector<std::string> types;
    std::vector<std::string> from_nodes;
    std::vector<std::string> to_nodes;
    std::vector<double> flow_rate_m3s;
    std::vector<double> velocity_ms;
    std::vector<double> headloss_m;
    std::vector<double> status;
};

struct StaciSteadyQualityResults {
    std::vector<std::string> node_ids;
    std::vector<std::string> link_ids;
    std::vector<double> node_water_age_s;
    std::vector<double> link_average_water_age_s;
    std::vector<double> link_travel_time_s;
    std::vector<double> node_concentration_kgm3;
    std::vector<double> link_average_concentration_kgm3;
};

struct StaciSensitivityResults {
    std::string parameter_element_id;
    std::string parameter_property;
    std::vector<std::string> link_ids;
    std::vector<std::string> node_ids;
    std::vector<double> link_flow_rate_derivative;
    std::vector<double> node_pressure_head_derivative;
};

// Stable, MATLAB-independent facade around the legacy STACI object model.
// All public property names and values use SI units.
class StaciSession {
public:
    explicit StaciSession(const std::string &network_file);

    const StaciSessionInfo &info() const noexcept { return info_; }
    void reset_hydraulic_state();
    StaciHydraulicStatus solve_hydraulics();

    std::vector<double> get_node_property(
        const std::vector<std::string> &ids,
        const std::string &property) const;
    void set_node_property(const std::vector<std::string> &ids,
                           const std::string &property,
                           const std::vector<double> &values);
    std::vector<double> get_link_property(
        const std::vector<std::string> &ids,
        const std::string &property) const;
    void set_link_property(const std::vector<std::string> &ids,
                           const std::string &property,
                           const std::vector<double> &values);

    StaciNodeResults node_results() const;
    StaciLinkResults link_results() const;
    StaciSteadyQualityResults solve_steady_quality(bool solve_age,
                                                   bool solve_chemical) const;
    StaciSensitivityResults hydraulic_sensitivity(
        const std::string &element_id,
        const std::string &property);

private:
    Csomopont &node(const std::string &id) const;
    Agelem &link(const std::string &id) const;
    static double value_for(const std::vector<double> &values,
                            std::size_t index, std::size_t expected);
    void require_solved() const;

    std::unique_ptr<Staci> system_;
    StaciSessionInfo info_;
    std::unordered_map<std::string, Csomopont *> nodes_;
    std::unordered_map<std::string, Agelem *> links_;
    bool solved_ = false;
};

#endif
