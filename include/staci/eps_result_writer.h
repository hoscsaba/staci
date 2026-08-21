#ifndef STACI_EPS_RESULT_WRITER_H
#define STACI_EPS_RESULT_WRITER_H

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

struct EpsNodeInfo {
    std::string id;
    std::string type;
    double elevation_m;
    double x_m;
    double y_m;
};

struct EpsLinkInfo {
    std::string id;
    std::string type;
    std::uint32_t from_node;
    std::uint32_t to_node;
    double length_m;
    double diameter_m;
};

struct EpsTankInfo {
    std::string id;
    std::uint32_t node_index;
    double min_level_m;
    double max_level_m;
    double diameter_m;
};

struct EpsOutputMetadata {
    std::string input_file;
    std::int64_t duration_s;
    std::int64_t hydraulic_timestep_s;
    std::int64_t effective_timestep_s;
    std::int64_t pattern_timestep_s;
    std::int64_t report_timestep_s;
};

struct EpsResultFrame {
    std::int64_t time_s;
    std::vector<double> node_head_m;
    std::vector<double> node_pressure_head_m;
    std::vector<double> node_demand_m3s;
    std::vector<double> link_flow_rate_m3s;
    std::vector<double> link_velocity_ms;
    std::vector<double> link_headloss_m;
    std::vector<std::uint8_t> link_status;
    std::vector<double> tank_level_m;
    std::vector<double> tank_volume_m3;
    std::vector<double> tank_inflow_m3s;
    bool converged;
    std::uint32_t iterations;
};

class EpsResultWriter {
public:
    EpsResultWriter(const std::string &output_prefix,
                    const std::vector<EpsNodeInfo> &nodes,
                    const std::vector<EpsLinkInfo> &links,
                    const std::vector<EpsTankInfo> &tanks,
                    const EpsOutputMetadata &metadata,
                    std::size_t chunk_frames = 16);
    ~EpsResultWriter();

    EpsResultWriter(const EpsResultWriter &) = delete;
    EpsResultWriter &operator=(const EpsResultWriter &) = delete;

    void append(const EpsResultFrame &frame);
    void finish(std::size_t hydraulic_states,
                std::size_t failed_states,
                const std::vector<std::string> &warnings);
    bool hdf5_enabled() const;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

#endif
