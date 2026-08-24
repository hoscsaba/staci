#include "eps_result_writer.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>

#ifdef STACI_HAS_HDF5
#include <hdf5.h>
#endif

namespace {

std::string json_escape(const std::string &value) {
    std::ostringstream out;
    for (unsigned char character : value) {
        switch (character) {
        case '"': out << "\\\""; break;
        case '\\': out << "\\\\"; break;
        case '\b': out << "\\b"; break;
        case '\f': out << "\\f"; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default:
            if (character < 0x20)
                out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                    << static_cast<unsigned int>(character) << std::dec;
            else
                out << character;
        }
    }
    return out.str();
}

struct Range {
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();

    void include(const std::vector<double> &values) {
        for (double value : values) {
            if (!std::isfinite(value))
                continue;
            minimum = std::min(minimum, value);
            maximum = std::max(maximum, value);
        }
    }

    bool valid() const { return minimum <= maximum; }
};

std::string utc_timestamp() {
    const std::time_t now = std::time(nullptr);
    std::tm utc{};
#ifdef _WIN32
    gmtime_s(&utc, &now);
#else
    gmtime_r(&now, &utc);
#endif
    std::ostringstream output;
    output << std::put_time(&utc, "%Y-%m-%dT%H:%M:%SZ");
    return output.str();
}

#ifdef STACI_HAS_HDF5

void h5_check(herr_t status, const std::string &operation) {
    if (status < 0)
        throw std::runtime_error("HDF5 error while " + operation + ".");
}

hid_t h5_check_id(hid_t id, const std::string &operation) {
    if (id < 0)
        throw std::runtime_error("HDF5 error while " + operation + ".");
    return id;
}

void write_string_attribute(hid_t object, const char *name, const std::string &value) {
    const hid_t type = h5_check_id(H5Tcopy(H5T_C_S1), "creating string attribute type");
    h5_check(H5Tset_size(type, H5T_VARIABLE), "configuring string attribute type");
    h5_check(H5Tset_cset(type, H5T_CSET_UTF8), "configuring UTF-8 attribute");
    const hid_t space = h5_check_id(H5Screate(H5S_SCALAR), "creating attribute dataspace");
    hid_t attribute = -1;
    if (H5Aexists(object, name) > 0)
        attribute = h5_check_id(H5Aopen(object, name, H5P_DEFAULT), "opening attribute");
    else
        attribute = h5_check_id(H5Acreate2(object, name, type, space,
                                           H5P_DEFAULT, H5P_DEFAULT),
                                "creating attribute");
    const char *text = value.c_str();
    h5_check(H5Awrite(attribute, type, &text), "writing attribute");
    H5Aclose(attribute);
    H5Sclose(space);
    H5Tclose(type);
}

template <typename T>
void write_scalar_attribute(hid_t object, const char *name, hid_t type, T value) {
    const hid_t space = h5_check_id(H5Screate(H5S_SCALAR), "creating attribute dataspace");
    hid_t attribute = -1;
    if (H5Aexists(object, name) > 0)
        attribute = h5_check_id(H5Aopen(object, name, H5P_DEFAULT), "opening attribute");
    else
        attribute = h5_check_id(H5Acreate2(object, name, type, space,
                                           H5P_DEFAULT, H5P_DEFAULT),
                                "creating attribute");
    h5_check(H5Awrite(attribute, type, &value), "writing attribute");
    H5Aclose(attribute);
    H5Sclose(space);
}

hid_t create_group(hid_t parent, const char *name) {
    return h5_check_id(H5Gcreate2(parent, name, H5P_DEFAULT, H5P_DEFAULT,
                                  H5P_DEFAULT),
                       std::string("creating group ") + name);
}

hid_t create_dynamic_1d(hid_t parent, const char *name, hid_t type,
                        hsize_t chunk_frames, const char *unit,
                        const char *description) {
    const hsize_t dimensions[1] = {0};
    const hsize_t maximum[1] = {H5S_UNLIMITED};
    const hsize_t chunks[1] = {chunk_frames};
    const hid_t space = h5_check_id(H5Screate_simple(1, dimensions, maximum),
                                     "creating dynamic dataspace");
    const hid_t properties = h5_check_id(H5Pcreate(H5P_DATASET_CREATE),
                                          "creating dataset properties");
    h5_check(H5Pset_chunk(properties, 1, chunks), "setting one-dimensional chunks");
    const hid_t dataset = h5_check_id(H5Dcreate2(parent, name, type, space,
                                                  H5P_DEFAULT, properties,
                                                  H5P_DEFAULT),
                                       std::string("creating dataset ") + name);
    write_string_attribute(dataset, "unit", unit);
    write_string_attribute(dataset, "description", description);
    H5Pclose(properties);
    H5Sclose(space);
    return dataset;
}

hid_t create_dynamic_2d(hid_t parent, const char *name, hid_t type,
                        hsize_t columns, hsize_t chunk_frames,
                        const char *unit, const char *description) {
    if (columns == 0)
        return -1;
    const hsize_t dimensions[2] = {0, columns};
    const hsize_t maximum[2] = {H5S_UNLIMITED, columns};
    const hsize_t chunks[2] = {chunk_frames, columns};
    const hid_t space = h5_check_id(H5Screate_simple(2, dimensions, maximum),
                                     "creating dynamic matrix dataspace");
    const hid_t properties = h5_check_id(H5Pcreate(H5P_DATASET_CREATE),
                                          "creating matrix dataset properties");
    h5_check(H5Pset_chunk(properties, 2, chunks), "setting matrix chunks");
    if (H5Zfilter_avail(H5Z_FILTER_DEFLATE) > 0) {
        h5_check(H5Pset_shuffle(properties), "enabling shuffle filter");
        h5_check(H5Pset_deflate(properties, 1), "enabling gzip compression");
    }
    const hid_t dataset = h5_check_id(H5Dcreate2(parent, name, type, space,
                                                  H5P_DEFAULT, properties,
                                                  H5P_DEFAULT),
                                       std::string("creating dataset ") + name);
    write_string_attribute(dataset, "unit", unit);
    write_string_attribute(dataset, "description", description);
    write_string_attribute(dataset, "array_order", "time,element");
    H5Pclose(properties);
    H5Sclose(space);
    return dataset;
}

template <typename T>
void append_1d(hid_t dataset, hid_t memory_type, hsize_t old_rows,
               const std::vector<T> &values) {
    if (dataset < 0 || values.empty())
        return;
    const hsize_t new_dimensions[1] = {old_rows + values.size()};
    h5_check(H5Dset_extent(dataset, new_dimensions), "extending vector dataset");
    const hid_t file_space = h5_check_id(H5Dget_space(dataset), "opening vector dataspace");
    const hsize_t start[1] = {old_rows};
    const hsize_t count[1] = {values.size()};
    h5_check(H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, nullptr,
                                 count, nullptr),
             "selecting vector output block");
    const hid_t memory_space = h5_check_id(H5Screate_simple(1, count, nullptr),
                                            "creating vector memory space");
    h5_check(H5Dwrite(dataset, memory_type, memory_space, file_space,
                      H5P_DEFAULT, values.data()),
             "writing vector output block");
    H5Sclose(memory_space);
    H5Sclose(file_space);
}

template <typename T>
void append_2d(hid_t dataset, hid_t memory_type, hsize_t old_rows,
               hsize_t columns, const std::vector<T> &values) {
    if (dataset < 0 || columns == 0 || values.empty())
        return;
    const hsize_t rows = values.size() / columns;
    const hsize_t new_dimensions[2] = {old_rows + rows, columns};
    h5_check(H5Dset_extent(dataset, new_dimensions), "extending matrix dataset");
    const hid_t file_space = h5_check_id(H5Dget_space(dataset), "opening matrix dataspace");
    const hsize_t start[2] = {old_rows, 0};
    const hsize_t count[2] = {rows, columns};
    h5_check(H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, nullptr,
                                 count, nullptr),
             "selecting matrix output block");
    const hid_t memory_space = h5_check_id(H5Screate_simple(2, count, nullptr),
                                            "creating matrix memory space");
    h5_check(H5Dwrite(dataset, memory_type, memory_space, file_space,
                      H5P_DEFAULT, values.data()),
             "writing matrix output block");
    H5Sclose(memory_space);
    H5Sclose(file_space);
}

template <typename T>
void write_vector(hid_t parent, const char *name, hid_t file_type,
                  hid_t memory_type,
                  const std::vector<T> &values, const char *unit,
                  const char *description) {
    const hsize_t dimensions[1] = {values.size()};
    const hid_t space = h5_check_id(H5Screate_simple(1, dimensions, nullptr),
                                     "creating static vector dataspace");
    const hid_t dataset = h5_check_id(H5Dcreate2(parent, name, file_type, space,
                                                  H5P_DEFAULT, H5P_DEFAULT,
                                                  H5P_DEFAULT),
                                       std::string("creating dataset ") + name);
    if (!values.empty())
        h5_check(H5Dwrite(dataset, memory_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          values.data()),
                 "writing static vector");
    write_string_attribute(dataset, "unit", unit);
    write_string_attribute(dataset, "description", description);
    H5Dclose(dataset);
    H5Sclose(space);
}

void write_strings(hid_t parent, const char *name,
                   const std::vector<std::string> &values,
                   const char *description) {
    const hsize_t dimensions[1] = {values.size()};
    const hid_t space = h5_check_id(H5Screate_simple(1, dimensions, nullptr),
                                     "creating string dataspace");
    const hid_t type = h5_check_id(H5Tcopy(H5T_C_S1), "creating string type");
    h5_check(H5Tset_size(type, H5T_VARIABLE), "configuring variable strings");
    h5_check(H5Tset_cset(type, H5T_CSET_UTF8), "configuring UTF-8 strings");
    const hid_t dataset = h5_check_id(H5Dcreate2(parent, name, type, space,
                                                  H5P_DEFAULT, H5P_DEFAULT,
                                                  H5P_DEFAULT),
                                       std::string("creating string dataset ") + name);
    if (!values.empty()) {
        std::vector<const char *> pointers;
        pointers.reserve(values.size());
        for (const std::string &value : values)
            pointers.push_back(value.c_str());
        h5_check(H5Dwrite(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          pointers.data()),
                 "writing strings");
    }
    write_string_attribute(dataset, "description", description);
    H5Dclose(dataset);
    H5Tclose(type);
    H5Sclose(space);
}

#endif

} // namespace

class EpsResultWriter::Impl {
public:
    Impl(std::string output_prefix,
         const std::vector<EpsNodeInfo> &nodes,
         const std::vector<EpsLinkInfo> &links,
         const std::vector<EpsTankInfo> &tanks,
         EpsOutputMetadata metadata,
         std::size_t chunk_frames)
        : prefix_(std::move(output_prefix)), nodes_(nodes), links_(links), tanks_(tanks),
          metadata_(std::move(metadata)), chunk_frames_(std::max<std::size_t>(1, chunk_frames)) {
#ifdef STACI_HAS_HDF5
        open_hdf5();
#endif
    }

    ~Impl() {
#ifdef STACI_HAS_HDF5
        close_hdf5();
#endif
    }

    void append(const EpsResultFrame &frame) {
        validate(frame);
        ++frame_count_;
        if (!frame.converged)
            ++failed_frame_count_;
        if (frame.converged) {
            ranges_["head"].include(frame.node_head_m);
            ranges_["pressure_head"].include(frame.node_pressure_head_m);
            ranges_["demand"].include(frame.node_demand_m3s);
            ranges_["node_water_age"].include(frame.node_water_age_s);
            ranges_["node_chlorine"].include(frame.node_chlorine_kgm3);
            ranges_["flow_rate"].include(frame.link_flow_rate_m3s);
            ranges_["velocity"].include(frame.link_velocity_ms);
            ranges_["headloss"].include(frame.link_headloss_m);
            ranges_["link_water_age"].include(frame.link_water_age_s);
            ranges_["link_chlorine"].include(frame.link_chlorine_kgm3);
            ranges_["tank_level"].include(frame.tank_level_m);
            ranges_["tank_volume"].include(frame.tank_volume_m3);
            ranges_["tank_inflow"].include(frame.tank_inflow_m3s);
        }
#ifdef STACI_HAS_HDF5
        pending_.push_back(frame);
        if (pending_.size() >= chunk_frames_)
            flush();
#else
        (void)frame;
#endif
    }

    void finish(std::size_t hydraulic_states, std::size_t failed_states,
                const std::vector<std::string> &warnings) {
#ifdef STACI_HAS_HDF5
        flush();
        write_strings(simulation_group_, "warnings", warnings,
                      "Compatibility and simulation warnings");
        write_string_attribute(file_, "simulation_status", "complete");
        write_scalar_attribute(file_, "frames_written", H5T_NATIVE_ULLONG,
                               static_cast<unsigned long long>(written_frames_));
        h5_check(H5Fflush(file_, H5F_SCOPE_GLOBAL), "flushing completed EPS file");
#endif
        write_json(hydraulic_states, failed_states, warnings);
    }

    bool hdf5_enabled() const {
#ifdef STACI_HAS_HDF5
        return true;
#else
        return false;
#endif
    }

private:
    std::string prefix_;
    std::vector<EpsNodeInfo> nodes_;
    std::vector<EpsLinkInfo> links_;
    std::vector<EpsTankInfo> tanks_;
    EpsOutputMetadata metadata_;
    std::size_t chunk_frames_;
    std::map<std::string, Range> ranges_;
    std::size_t frame_count_ = 0;
    std::size_t failed_frame_count_ = 0;
    const std::string created_utc_ = utc_timestamp();

#ifdef STACI_HAS_HDF5
    hid_t file_ = -1;
    hid_t nodes_group_ = -1;
    hid_t links_group_ = -1;
    hid_t tanks_group_ = -1;
    hid_t simulation_group_ = -1;
    hid_t time_dataset_ = -1;
    hid_t node_head_dataset_ = -1;
    hid_t node_pressure_dataset_ = -1;
    hid_t node_demand_dataset_ = -1;
    hid_t node_water_age_dataset_ = -1;
    hid_t node_chlorine_dataset_ = -1;
    hid_t link_flow_dataset_ = -1;
    hid_t link_velocity_dataset_ = -1;
    hid_t link_headloss_dataset_ = -1;
    hid_t link_status_dataset_ = -1;
    hid_t link_water_age_dataset_ = -1;
    hid_t link_chlorine_dataset_ = -1;
    hid_t tank_level_dataset_ = -1;
    hid_t tank_volume_dataset_ = -1;
    hid_t tank_inflow_dataset_ = -1;
    hid_t converged_dataset_ = -1;
    hid_t iterations_dataset_ = -1;
    hsize_t written_frames_ = 0;
    std::vector<EpsResultFrame> pending_;

    void open_hdf5() {
        file_ = h5_check_id(H5Fcreate((prefix_ + ".h5").c_str(), H5F_ACC_TRUNC,
                                      H5P_DEFAULT, H5P_DEFAULT),
                            "creating EPS HDF5 file");
        write_string_attribute(file_, "format_name", "STACI EPS OUTPUT");
        write_string_attribute(file_, "format_version", "1.0");
        write_string_attribute(file_, "generator", "STACI");
        write_string_attribute(file_, "generator_version", "2.0");
        write_string_attribute(file_, "created_utc", created_utc_);
        write_string_attribute(file_, "input_file", metadata_.input_file);
        write_string_attribute(file_, "array_order", "time,element");
        write_string_attribute(file_, "unit_system", "SI");
        write_string_attribute(file_, "simulation_status", "in_progress");
        write_scalar_attribute(file_, "frames_written", H5T_NATIVE_ULLONG, 0ULL);
        write_scalar_attribute(file_, "node_count", H5T_NATIVE_ULLONG,
                               static_cast<unsigned long long>(nodes_.size()));
        write_scalar_attribute(file_, "link_count", H5T_NATIVE_ULLONG,
                               static_cast<unsigned long long>(links_.size()));
        write_scalar_attribute(file_, "tank_count", H5T_NATIVE_ULLONG,
                               static_cast<unsigned long long>(tanks_.size()));

        nodes_group_ = create_group(file_, "nodes");
        links_group_ = create_group(file_, "links");
        tanks_group_ = create_group(file_, "tanks");
        simulation_group_ = create_group(file_, "simulation");
        write_static_data();

        const hsize_t chunk = static_cast<hsize_t>(chunk_frames_);
        time_dataset_ = create_dynamic_1d(file_, "time", H5T_STD_I64LE, chunk,
                                          "s", "Elapsed simulation time");
        node_head_dataset_ = create_dynamic_2d(nodes_group_, "head", H5T_IEEE_F64LE,
                                               nodes_.size(), chunk, "m",
                                               "Total hydraulic head");
        node_pressure_dataset_ = create_dynamic_2d(nodes_group_, "pressure_head",
                                                   H5T_IEEE_F64LE, nodes_.size(), chunk,
                                                   "m", "Pressure head");
        node_demand_dataset_ = create_dynamic_2d(nodes_group_, "demand", H5T_IEEE_F64LE,
                                                 nodes_.size(), chunk, "m3/s",
                                                 "Node demand; positive means withdrawal");
        node_water_age_dataset_ = create_dynamic_2d(
            nodes_group_, "water_age", H5T_IEEE_F64LE, nodes_.size(), chunk,
            "s", "Water residence time; NaN unless EPANET QUALITY AGE is enabled");
        node_chlorine_dataset_ = create_dynamic_2d(
            nodes_group_, "chlorine", H5T_IEEE_F64LE, nodes_.size(), chunk,
            "kg/m3", "Chemical concentration in SI; NaN unless EPANET QUALITY CHEMICAL is enabled");
        link_flow_dataset_ = create_dynamic_2d(links_group_, "flow_rate", H5T_IEEE_F64LE,
                                               links_.size(), chunk, "m3/s",
                                               "Signed flow from from_node to to_node");
        link_velocity_dataset_ = create_dynamic_2d(links_group_, "velocity", H5T_IEEE_F64LE,
                                                   links_.size(), chunk, "m/s",
                                                   "Signed mean velocity from from_node to to_node");
        link_headloss_dataset_ = create_dynamic_2d(links_group_, "headloss", H5T_IEEE_F64LE,
                                                   links_.size(), chunk, "m",
                                                   "Signed head(from_node)-head(to_node)");
        link_status_dataset_ = create_dynamic_2d(links_group_, "status", H5T_STD_U8LE,
                                                 links_.size(), chunk, "1",
                                                 "0=CLOSED, 1=OPEN, 2=ACTIVE, 255=UNKNOWN");
        link_water_age_dataset_ = create_dynamic_2d(
            links_group_, "water_age", H5T_IEEE_F64LE, links_.size(), chunk,
            "s", "Volume-averaged link water residence time");
        link_chlorine_dataset_ = create_dynamic_2d(
            links_group_, "chlorine", H5T_IEEE_F64LE, links_.size(), chunk,
            "kg/m3", "Volume-averaged chemical concentration in SI");
        tank_level_dataset_ = create_dynamic_2d(tanks_group_, "level", H5T_IEEE_F64LE,
                                                tanks_.size(), chunk, "m",
                                                "Water level above tank bottom");
        tank_volume_dataset_ = create_dynamic_2d(tanks_group_, "volume", H5T_IEEE_F64LE,
                                                 tanks_.size(), chunk, "m3",
                                                 "Stored water volume");
        tank_inflow_dataset_ = create_dynamic_2d(tanks_group_, "inflow", H5T_IEEE_F64LE,
                                                 tanks_.size(), chunk, "m3/s",
                                                 "Signed inflow into tank");
        converged_dataset_ = create_dynamic_1d(simulation_group_, "converged",
                                               H5T_STD_U8LE, chunk, "1",
                                               "1 when the hydraulic state converged");
        iterations_dataset_ = create_dynamic_1d(simulation_group_, "iterations",
                                                H5T_STD_U32LE, chunk, "1",
                                                "Hydraulic nonlinear iteration count; 0 when unavailable");
    }

    void write_static_data() {
        std::vector<std::string> node_ids, node_types;
        std::vector<double> elevations, x, y;
        for (const EpsNodeInfo &node : nodes_) {
            node_ids.push_back(node.id);
            node_types.push_back(node.type);
            elevations.push_back(node.elevation_m);
            x.push_back(node.x_m);
            y.push_back(node.y_m);
        }
        write_strings(nodes_group_, "id", node_ids, "EPANET node identifier");
        write_strings(nodes_group_, "type", node_types, "JUNCTION, RESERVOIR, or TANK");
        write_vector(nodes_group_, "elevation", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, elevations, "m",
                     "Node elevation in SI units");
        write_vector(nodes_group_, "x", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, x, "m",
                     "EPANET x coordinate converted with the network length unit");
        write_vector(nodes_group_, "y", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, y, "m",
                     "EPANET y coordinate converted with the network length unit");

        std::vector<std::string> link_ids, link_types;
        std::vector<std::uint32_t> from, to;
        std::vector<double> lengths, diameters;
        for (const EpsLinkInfo &link : links_) {
            link_ids.push_back(link.id);
            link_types.push_back(link.type);
            from.push_back(link.from_node);
            to.push_back(link.to_node);
            lengths.push_back(link.length_m);
            diameters.push_back(link.diameter_m);
        }
        write_strings(links_group_, "id", link_ids, "EPANET link identifier");
        write_strings(links_group_, "type", link_types, "STACI link type");
        write_vector(links_group_, "from_node", H5T_STD_U32LE, H5T_NATIVE_UINT32, from, "1",
                     "Index into nodes/id");
        write_vector(links_group_, "to_node", H5T_STD_U32LE, H5T_NATIVE_UINT32, to, "1",
                     "Index into nodes/id");
        write_vector(links_group_, "length", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, lengths, "m",
                     "Link length; NaN when not applicable");
        write_vector(links_group_, "diameter", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, diameters, "m",
                     "Link hydraulic diameter; NaN when not applicable");

        std::vector<std::string> tank_ids;
        std::vector<std::uint32_t> tank_nodes;
        std::vector<double> minimum, maximum, diameter;
        for (const EpsTankInfo &tank : tanks_) {
            tank_ids.push_back(tank.id);
            tank_nodes.push_back(tank.node_index);
            minimum.push_back(tank.min_level_m);
            maximum.push_back(tank.max_level_m);
            diameter.push_back(tank.diameter_m);
        }
        write_strings(tanks_group_, "id", tank_ids, "EPANET tank identifier");
        write_vector(tanks_group_, "node_index", H5T_STD_U32LE, H5T_NATIVE_UINT32, tank_nodes, "1",
                     "Index into nodes/id");
        write_vector(tanks_group_, "min_level", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, minimum, "m",
                     "Minimum water level above tank bottom");
        write_vector(tanks_group_, "max_level", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, maximum, "m",
                     "Maximum water level above tank bottom");
        write_vector(tanks_group_, "diameter", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, diameter, "m",
                     "Equivalent cylindrical tank diameter");

        write_scalar_attribute(simulation_group_, "duration_seconds", H5T_NATIVE_LLONG,
                               static_cast<long long>(metadata_.duration_s));
        write_scalar_attribute(simulation_group_, "hydraulic_timestep_seconds", H5T_NATIVE_LLONG,
                               static_cast<long long>(metadata_.hydraulic_timestep_s));
        write_scalar_attribute(simulation_group_, "effective_timestep_seconds", H5T_NATIVE_LLONG,
                               static_cast<long long>(metadata_.effective_timestep_s));
        write_scalar_attribute(simulation_group_, "pattern_timestep_seconds", H5T_NATIVE_LLONG,
                               static_cast<long long>(metadata_.pattern_timestep_s));
        write_scalar_attribute(simulation_group_, "report_timestep_seconds", H5T_NATIVE_LLONG,
                               static_cast<long long>(metadata_.report_timestep_s));
        write_scalar_attribute(simulation_group_, "quality_timestep_seconds", H5T_NATIVE_LLONG,
                               static_cast<long long>(metadata_.quality_timestep_s));
        write_scalar_attribute(simulation_group_, "water_age_enabled", H5T_NATIVE_UCHAR,
                               static_cast<unsigned char>(metadata_.water_age_enabled ? 1 : 0));
        write_scalar_attribute(simulation_group_, "chemical_enabled", H5T_NATIVE_UCHAR,
                               static_cast<unsigned char>(metadata_.chemical_enabled ? 1 : 0));
        write_string_attribute(simulation_group_, "chemical_name", metadata_.chemical_name);
    }

    void flush() {
        if (pending_.empty())
            return;
        std::vector<std::int64_t> time;
        std::vector<double> node_head, node_pressure, node_demand, node_water_age, node_chlorine;
        std::vector<double> link_flow, link_velocity, link_headloss, link_water_age, link_chlorine;
        std::vector<std::uint8_t> link_status, converged;
        std::vector<double> tank_level, tank_volume, tank_inflow;
        std::vector<std::uint32_t> iterations;
        for (const EpsResultFrame &frame : pending_) {
            time.push_back(frame.time_s);
            node_head.insert(node_head.end(), frame.node_head_m.begin(), frame.node_head_m.end());
            node_pressure.insert(node_pressure.end(), frame.node_pressure_head_m.begin(),
                                 frame.node_pressure_head_m.end());
            node_demand.insert(node_demand.end(), frame.node_demand_m3s.begin(),
                               frame.node_demand_m3s.end());
            node_water_age.insert(node_water_age.end(), frame.node_water_age_s.begin(),
                                  frame.node_water_age_s.end());
            node_chlorine.insert(node_chlorine.end(), frame.node_chlorine_kgm3.begin(),
                                 frame.node_chlorine_kgm3.end());
            link_flow.insert(link_flow.end(), frame.link_flow_rate_m3s.begin(),
                             frame.link_flow_rate_m3s.end());
            link_velocity.insert(link_velocity.end(), frame.link_velocity_ms.begin(),
                                 frame.link_velocity_ms.end());
            link_headloss.insert(link_headloss.end(), frame.link_headloss_m.begin(),
                                 frame.link_headloss_m.end());
            link_status.insert(link_status.end(), frame.link_status.begin(), frame.link_status.end());
            link_water_age.insert(link_water_age.end(), frame.link_water_age_s.begin(),
                                  frame.link_water_age_s.end());
            link_chlorine.insert(link_chlorine.end(), frame.link_chlorine_kgm3.begin(),
                                 frame.link_chlorine_kgm3.end());
            tank_level.insert(tank_level.end(), frame.tank_level_m.begin(), frame.tank_level_m.end());
            tank_volume.insert(tank_volume.end(), frame.tank_volume_m3.begin(), frame.tank_volume_m3.end());
            tank_inflow.insert(tank_inflow.end(), frame.tank_inflow_m3s.begin(), frame.tank_inflow_m3s.end());
            converged.push_back(frame.converged ? 1 : 0);
            iterations.push_back(frame.iterations);
        }
        append_1d(time_dataset_, H5T_NATIVE_LLONG, written_frames_, time);
        append_2d(node_head_dataset_, H5T_NATIVE_DOUBLE, written_frames_, nodes_.size(), node_head);
        append_2d(node_pressure_dataset_, H5T_NATIVE_DOUBLE, written_frames_, nodes_.size(), node_pressure);
        append_2d(node_demand_dataset_, H5T_NATIVE_DOUBLE, written_frames_, nodes_.size(), node_demand);
        append_2d(node_water_age_dataset_, H5T_NATIVE_DOUBLE, written_frames_, nodes_.size(), node_water_age);
        append_2d(node_chlorine_dataset_, H5T_NATIVE_DOUBLE, written_frames_, nodes_.size(), node_chlorine);
        append_2d(link_flow_dataset_, H5T_NATIVE_DOUBLE, written_frames_, links_.size(), link_flow);
        append_2d(link_velocity_dataset_, H5T_NATIVE_DOUBLE, written_frames_, links_.size(), link_velocity);
        append_2d(link_headloss_dataset_, H5T_NATIVE_DOUBLE, written_frames_, links_.size(), link_headloss);
        append_2d(link_status_dataset_, H5T_NATIVE_UCHAR, written_frames_, links_.size(), link_status);
        append_2d(link_water_age_dataset_, H5T_NATIVE_DOUBLE, written_frames_, links_.size(), link_water_age);
        append_2d(link_chlorine_dataset_, H5T_NATIVE_DOUBLE, written_frames_, links_.size(), link_chlorine);
        append_2d(tank_level_dataset_, H5T_NATIVE_DOUBLE, written_frames_, tanks_.size(), tank_level);
        append_2d(tank_volume_dataset_, H5T_NATIVE_DOUBLE, written_frames_, tanks_.size(), tank_volume);
        append_2d(tank_inflow_dataset_, H5T_NATIVE_DOUBLE, written_frames_, tanks_.size(), tank_inflow);
        append_1d(converged_dataset_, H5T_NATIVE_UCHAR, written_frames_, converged);
        append_1d(iterations_dataset_, H5T_NATIVE_UINT, written_frames_, iterations);
        written_frames_ += pending_.size();
        write_scalar_attribute(file_, "frames_written", H5T_NATIVE_ULLONG,
                               static_cast<unsigned long long>(written_frames_));
        h5_check(H5Fflush(file_, H5F_SCOPE_LOCAL), "flushing EPS output chunk");
        pending_.clear();
    }

    void close_hdf5() {
        const hid_t datasets[] = {
            time_dataset_, node_head_dataset_, node_pressure_dataset_, node_demand_dataset_,
            node_water_age_dataset_, node_chlorine_dataset_, link_flow_dataset_, link_velocity_dataset_,
            link_headloss_dataset_, link_status_dataset_, link_water_age_dataset_, link_chlorine_dataset_,
            tank_level_dataset_, tank_volume_dataset_, tank_inflow_dataset_, converged_dataset_,
            iterations_dataset_};
        for (hid_t dataset : datasets)
            if (dataset >= 0)
                H5Dclose(dataset);
        const hid_t groups[] = {nodes_group_, links_group_, tanks_group_, simulation_group_};
        for (hid_t group : groups)
            if (group >= 0)
                H5Gclose(group);
        if (file_ >= 0)
            H5Fclose(file_);
        file_ = -1;
    }
#endif

    void validate(const EpsResultFrame &frame) const {
        const bool node_sizes = frame.node_head_m.size() == nodes_.size() &&
                                frame.node_pressure_head_m.size() == nodes_.size() &&
                                frame.node_demand_m3s.size() == nodes_.size() &&
                                frame.node_water_age_s.size() == nodes_.size() &&
                                frame.node_chlorine_kgm3.size() == nodes_.size();
        const bool link_sizes = frame.link_flow_rate_m3s.size() == links_.size() &&
                                frame.link_velocity_ms.size() == links_.size() &&
                                frame.link_headloss_m.size() == links_.size() &&
                                frame.link_status.size() == links_.size() &&
                                frame.link_water_age_s.size() == links_.size() &&
                                frame.link_chlorine_kgm3.size() == links_.size();
        const bool tank_sizes = frame.tank_level_m.size() == tanks_.size() &&
                                frame.tank_volume_m3.size() == tanks_.size() &&
                                frame.tank_inflow_m3s.size() == tanks_.size();
        if (!node_sizes || !link_sizes || !tank_sizes)
            throw std::runtime_error("EPS result frame dimensions do not match the network index.");
    }

    void write_json(std::size_t hydraulic_states, std::size_t failed_states,
                    const std::vector<std::string> &warnings) const {
        std::ofstream output(prefix_ + ".meta.json", std::ios::trunc);
        if (!output)
            throw std::runtime_error("Cannot create EPS metadata JSON: " + prefix_ + ".meta.json");
        output << std::setprecision(15)
               << "{\n"
               << "  \"format\": {\"name\": \"STACI EPS OUTPUT\", \"version\": \"1.0\", "
               << "\"data_file\": ";
        if (hdf5_enabled())
            output << "\"" << json_escape(std::filesystem::path(prefix_ + ".h5").filename().string()) << "\"";
        else
            output << "null";
        output << ", \"hdf5_available\": " << (hdf5_enabled() ? "true" : "false") << "},\n"
               << "  \"generator\": {\"name\": \"STACI\", \"version\": \"2.0\"},\n"
               << "  \"created_utc\": \"" << created_utc_ << "\",\n"
               << "  \"source\": {\"input_file\": \"" << json_escape(metadata_.input_file) << "\"},\n"
               << "  \"units\": \"SI\",\n"
               << "  \"simulation\": {\n"
               << "    \"status\": \"complete\",\n"
               << "    \"duration_seconds\": " << metadata_.duration_s << ",\n"
               << "    \"hydraulic_timestep_seconds\": " << metadata_.hydraulic_timestep_s << ",\n"
               << "    \"effective_timestep_seconds\": " << metadata_.effective_timestep_s << ",\n"
               << "    \"pattern_timestep_seconds\": " << metadata_.pattern_timestep_s << ",\n"
               << "    \"report_timestep_seconds\": " << metadata_.report_timestep_s << ",\n"
               << "    \"quality_timestep_seconds\": " << metadata_.quality_timestep_s << ",\n"
               << "    \"water_age_enabled\": " << (metadata_.water_age_enabled ? "true" : "false") << ",\n"
               << "    \"chemical_enabled\": " << (metadata_.chemical_enabled ? "true" : "false") << ",\n"
               << "    \"chemical_name\": \"" << json_escape(metadata_.chemical_name) << "\",\n"
               << "    \"frames\": " << frame_count_ << ",\n"
               << "    \"hydraulic_states\": " << hydraulic_states << ",\n"
               << "    \"failed_frames\": " << failed_frame_count_ << ",\n"
               << "    \"failed_hydraulic_states\": " << failed_states << "\n"
               << "  },\n"
               << "  \"network\": {\"nodes\": " << nodes_.size()
               << ", \"links\": " << links_.size() << ", \"tanks\": " << tanks_.size() << "},\n"
               << "  \"ranges\": {\n";
        const std::vector<std::pair<std::string, std::string> > range_order = {
            {"head", "m"}, {"pressure_head", "m"}, {"demand", "m3/s"},
            {"node_water_age", "s"}, {"node_chlorine", "kg/m3"}, {"flow_rate", "m3/s"},
            {"velocity", "m/s"}, {"headloss", "m"}, {"link_water_age", "s"},
            {"link_chlorine", "kg/m3"},
            {"tank_level", "m"}, {"tank_volume", "m3"}, {"tank_inflow", "m3/s"}};
        for (std::size_t index = 0; index < range_order.size(); ++index) {
            const auto found = ranges_.find(range_order[index].first);
            output << "    \"" << range_order[index].first << "\": {\"unit\": \""
                   << range_order[index].second << "\", \"min\": ";
            if (found != ranges_.end() && found->second.valid())
                output << found->second.minimum << ", \"max\": " << found->second.maximum;
            else
                output << "null, \"max\": null";
            output << '}' << (index + 1 == range_order.size() ? "\n" : ",\n");
        }
        output << "  },\n"
               << "  \"status_codes\": {\"0\": \"closed\", \"1\": \"open\", "
               << "\"2\": \"active\", \"255\": \"unknown\"},\n"
               << "  \"warnings\": [";
        for (std::size_t index = 0; index < warnings.size(); ++index) {
            if (index != 0)
                output << ',';
            output << "\n    \"" << json_escape(warnings[index]) << "\"";
        }
        if (!warnings.empty())
            output << '\n';
        output << "  ]\n}\n";
    }
};

EpsResultWriter::EpsResultWriter(const std::string &output_prefix,
                                 const std::vector<EpsNodeInfo> &nodes,
                                 const std::vector<EpsLinkInfo> &links,
                                 const std::vector<EpsTankInfo> &tanks,
                                 const EpsOutputMetadata &metadata,
                                 std::size_t chunk_frames)
    : impl_(new Impl(output_prefix, nodes, links, tanks, metadata, chunk_frames)) {}

EpsResultWriter::~EpsResultWriter() = default;

void EpsResultWriter::append(const EpsResultFrame &frame) { impl_->append(frame); }

void EpsResultWriter::finish(std::size_t hydraulic_states,
                             std::size_t failed_states,
                             const std::vector<std::string> &warnings) {
    impl_->finish(hydraulic_states, failed_states, warnings);
}

bool EpsResultWriter::hdf5_enabled() const { return impl_->hdf5_enabled(); }
