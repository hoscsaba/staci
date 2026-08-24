#ifndef STACI_EPANET_WATER_AGE_H
#define STACI_EPANET_WATER_AGE_H

#include <cstddef>
#include <deque>
#include <vector>

struct EpanetWaterAgeNode {
    bool fixed_zero_age = false;
};

struct EpanetWaterAgeLink {
    std::size_t from_node = 0;
    std::size_t to_node = 0;
    double volume_m3 = 0.0;
};

enum class EpanetChemicalSourceType {
    None,
    Concentration,
    Mass,
    Setpoint,
    FlowPaced
};

struct EpanetChemicalNode {
    double initial_concentration_kgm3 = 0.0;
    bool fixed_external_concentration = false;
};

struct EpanetChemicalLink {
    std::size_t from_node = 0;
    std::size_t to_node = 0;
    double volume_m3 = 0.0;
    // First-order total reaction coefficient in SI [1/s]. Negative values
    // consume the chemical and positive values produce it.
    double reaction_coefficient_per_s = 0.0;
    double initial_concentration_kgm3 = 0.0;
};

struct EpanetChemicalSource {
    EpanetChemicalSourceType type = EpanetChemicalSourceType::None;
    // kg/m3 for CONCEN, SETPOINT, and FLOWPACED; kg/s for MASS.
    double strength_si = 0.0;
};

// Chemical transport uses the same segment displacement and instantaneous
// node-mixing model as the existing water-age calculation. Concentrations are
// stored exclusively in SI kg/m3.
class EpanetChemicalModel {
public:
    EpanetChemicalModel(std::vector<EpanetChemicalNode> nodes,
                        std::vector<EpanetChemicalLink> links,
                        double quality_timestep_s);
    void advance(double duration_s,
                 const std::vector<double> &link_flows_m3s,
                 const std::vector<double> &external_inflows_m3s,
                 const std::vector<EpanetChemicalSource> &sources);
    const std::vector<double> &node_concentration_kgm3() const {
        return node_concentration_kgm3_;
    }
    std::vector<double> link_average_concentration_kgm3(
        const std::vector<double> &link_flows_m3s) const;

private:
    struct Segment { double volume_m3; double concentration_kgm3; };
    std::vector<EpanetChemicalNode> nodes_;
    std::vector<EpanetChemicalLink> links_;
    std::vector<std::deque<Segment> > segments_;
    std::vector<double> node_concentration_kgm3_;
    double quality_timestep_s_;

    void step(double dt_s, const std::vector<double> &link_flows_m3s,
              const std::vector<double> &external_inflows_m3s,
              const std::vector<EpanetChemicalSource> &sources);
    void mix_nodes(const std::vector<double> &link_flows_m3s,
                   const std::vector<double> &external_inflows_m3s,
                   const std::vector<EpanetChemicalSource> &sources);
    double outlet_concentration(std::size_t link_index, double flow_m3s) const;
    void displace(std::size_t link_index, double flow_m3s, double dt_s,
                  double inlet_concentration_kgm3);
};

class EpanetWaterAgeModel {
public:
    EpanetWaterAgeModel(std::vector<EpanetWaterAgeNode> nodes,
                        std::vector<EpanetWaterAgeLink> links,
                        double quality_timestep_s);
    void advance(double duration_s, const std::vector<double> &link_flows_m3s);
    const std::vector<double> &node_age_s() const { return node_age_s_; }
    std::vector<double> link_average_age_s(
        const std::vector<double> &link_flows_m3s) const;

private:
    struct Segment { double volume_m3; double age_s; };
    std::vector<EpanetWaterAgeNode> nodes_;
    std::vector<EpanetWaterAgeLink> links_;
    std::vector<std::deque<Segment> > segments_;
    std::vector<double> node_age_s_;
    double quality_timestep_s_;

    void step(double dt_s, const std::vector<double> &link_flows_m3s);
    void mix_nodes(double elapsed_s, const std::vector<double> &link_flows_m3s);
    double outlet_age(std::size_t link_index, double flow_m3s) const;
    void displace(std::size_t link_index, double flow_m3s, double dt_s,
                  double inlet_age_s);
};

#endif
