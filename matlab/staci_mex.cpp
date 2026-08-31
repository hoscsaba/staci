#include "mex.h"

#include "StaciSession.h"

#include <algorithm>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

std::unordered_map<std::uint64_t, std::unique_ptr<StaciSession> > sessions;
std::uint64_t next_handle = 1;

void cleanup() {
    sessions.clear();
    while (mexIsLocked())
        mexUnlock();
}

std::string text_value(const mxArray *value, const char *name) {
    if (!mxIsChar(value))
        throw std::invalid_argument(std::string(name) + " must be text.");
    char *text = mxArrayToUTF8String(value);
    if (text == nullptr)
        throw std::runtime_error(std::string("Cannot decode ") + name + '.');
    std::string result(text);
    mxFree(text);
    return result;
}

std::vector<std::string> text_vector(const mxArray *value, const char *name) {
    if (mxIsChar(value))
        return {text_value(value, name)};
    if (!mxIsCell(value))
        throw std::invalid_argument(
            std::string(name) + " must be text or a cell array of text.");
    std::vector<std::string> result;
    result.reserve(mxGetNumberOfElements(value));
    for (std::size_t index = 0; index < mxGetNumberOfElements(value); ++index) {
        const mxArray *item = mxGetCell(value, index);
        if (item == nullptr)
            throw std::invalid_argument(std::string(name) + " contains an empty cell.");
        result.push_back(text_value(item, name));
    }
    return result;
}

std::vector<double> numeric_vector(const mxArray *value, const char *name) {
    if (!mxIsDouble(value) || mxIsComplex(value))
        throw std::invalid_argument(std::string(name) + " must be a real double array.");
    const double *data = mxGetPr(value);
    return std::vector<double>(data, data + mxGetNumberOfElements(value));
}

std::uint64_t handle_value(const mxArray *value) {
    if (!mxIsUint64(value) || mxGetNumberOfElements(value) != 1)
        throw std::invalid_argument("The STACI handle must be a scalar uint64.");
    return *static_cast<const std::uint64_t *>(mxGetData(value));
}

StaciSession &session(const mxArray *value) {
    const auto found = sessions.find(handle_value(value));
    if (found == sessions.end())
        throw std::invalid_argument("Unknown or already released STACI handle.");
    return *found->second;
}

mxArray *double_vector(const std::vector<double> &values) {
    mxArray *result = mxCreateDoubleMatrix(values.size(), 1, mxREAL);
    std::copy(values.begin(), values.end(), mxGetPr(result));
    return result;
}

mxArray *string_cells(const std::vector<std::string> &values) {
    mxArray *result = mxCreateCellMatrix(values.size(), 1);
    for (std::size_t index = 0; index < values.size(); ++index)
        mxSetCell(result, index, mxCreateString(values[index].c_str()));
    return result;
}

mxArray *info_struct(const StaciSessionInfo &info) {
    const char *fields[] = {"api_version", "node_ids", "link_ids", "link_types",
                            "link_from_nodes", "link_to_nodes"};
    mxArray *result = mxCreateStructMatrix(1, 1, 6, fields);
    mxSetField(result, 0, fields[0], mxCreateString("1.0"));
    mxSetField(result, 0, fields[1], string_cells(info.node_ids));
    mxSetField(result, 0, fields[2], string_cells(info.link_ids));
    mxSetField(result, 0, fields[3], string_cells(info.link_types));
    mxSetField(result, 0, fields[4], string_cells(info.link_from_nodes));
    mxSetField(result, 0, fields[5], string_cells(info.link_to_nodes));
    return result;
}

mxArray *node_results_struct(const StaciNodeResults &values) {
    const char *fields[] = {"id", "elevation_m", "pressure_head_m",
                            "total_head_m", "demand_m3s"};
    mxArray *result = mxCreateStructMatrix(1, 1, 5, fields);
    mxSetField(result, 0, fields[0], string_cells(values.ids));
    mxSetField(result, 0, fields[1], double_vector(values.elevation_m));
    mxSetField(result, 0, fields[2], double_vector(values.pressure_head_m));
    mxSetField(result, 0, fields[3], double_vector(values.total_head_m));
    mxSetField(result, 0, fields[4], double_vector(values.demand_m3s));
    return result;
}

mxArray *link_results_struct(const StaciLinkResults &values) {
    const char *fields[] = {"id", "type", "from_node", "to_node",
                            "flow_rate_m3s", "velocity_ms", "headloss_m", "status"};
    mxArray *result = mxCreateStructMatrix(1, 1, 8, fields);
    mxSetField(result, 0, fields[0], string_cells(values.ids));
    mxSetField(result, 0, fields[1], string_cells(values.types));
    mxSetField(result, 0, fields[2], string_cells(values.from_nodes));
    mxSetField(result, 0, fields[3], string_cells(values.to_nodes));
    mxSetField(result, 0, fields[4], double_vector(values.flow_rate_m3s));
    mxSetField(result, 0, fields[5], double_vector(values.velocity_ms));
    mxSetField(result, 0, fields[6], double_vector(values.headloss_m));
    mxSetField(result, 0, fields[7], double_vector(values.status));
    return result;
}

mxArray *quality_results_struct(const StaciSteadyQualityResults &values) {
    const char *fields[] = {"node_ids", "link_ids", "node_water_age_s",
                            "link_average_water_age_s", "link_travel_time_s",
                            "node_concentration_kgm3",
                            "link_average_concentration_kgm3"};
    mxArray *result = mxCreateStructMatrix(1, 1, 7, fields);
    mxSetField(result, 0, fields[0], string_cells(values.node_ids));
    mxSetField(result, 0, fields[1], string_cells(values.link_ids));
    mxSetField(result, 0, fields[2], double_vector(values.node_water_age_s));
    mxSetField(result, 0, fields[3], double_vector(values.link_average_water_age_s));
    mxSetField(result, 0, fields[4], double_vector(values.link_travel_time_s));
    mxSetField(result, 0, fields[5], double_vector(values.node_concentration_kgm3));
    mxSetField(result, 0, fields[6], double_vector(values.link_average_concentration_kgm3));
    return result;
}

mxArray *sensitivity_results_struct(const StaciSensitivityResults &values) {
    const char *fields[] = {"parameter_element_id", "parameter_property",
                            "link_ids", "node_ids",
                            "link_flow_rate_derivative",
                            "node_pressure_head_derivative"};
    mxArray *result = mxCreateStructMatrix(1, 1, 6, fields);
    mxSetField(result, 0, fields[0],
               mxCreateString(values.parameter_element_id.c_str()));
    mxSetField(result, 0, fields[1],
               mxCreateString(values.parameter_property.c_str()));
    mxSetField(result, 0, fields[2], string_cells(values.link_ids));
    mxSetField(result, 0, fields[3], string_cells(values.node_ids));
    mxSetField(result, 0, fields[4],
               double_vector(values.link_flow_rate_derivative));
    mxSetField(result, 0, fields[5],
               double_vector(values.node_pressure_head_derivative));
    return result;
}

void require_arguments(int actual, int expected, const char *usage) {
    if (actual != expected)
        throw std::invalid_argument(usage);
}

} // namespace

extern "C" void staciMexFunction(int nlhs, mxArray *plhs[], int nrhs,
                                 const mxArray *prhs[]) {
    try {
        if (nrhs < 1)
            throw std::invalid_argument("Usage: staci_mex(command, ...)");
        const std::string command = text_value(prhs[0], "command");
        if (command == "load") {
            require_arguments(nrhs, 2, "[handle,info] = staci_mex('load', networkFile)");
            if (nlhs < 1)
                throw std::invalid_argument("The load command requires an output handle.");
            if (sessions.empty())
                mexAtExit(cleanup);
            const std::uint64_t handle = next_handle++;
            auto value = std::make_unique<StaciSession>(
                text_value(prhs[1], "networkFile"));
            StaciSession *pointer = value.get();
            sessions.emplace(handle, std::move(value));
            mexLock();
            plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
            *static_cast<std::uint64_t *>(mxGetData(plhs[0])) = handle;
            if (nlhs > 1)
                plhs[1] = info_struct(pointer->info());
        } else if (command == "release") {
            require_arguments(nrhs, 2, "staci_mex('release', handle)");
            const std::uint64_t handle = handle_value(prhs[1]);
            if (sessions.erase(handle) != 0)
                mexUnlock();
        } else if (command == "info") {
            require_arguments(nrhs, 2, "info = staci_mex('info', handle)");
            plhs[0] = info_struct(session(prhs[1]).info());
        } else if (command == "reset") {
            require_arguments(nrhs, 2, "staci_mex('reset', handle)");
            session(prhs[1]).reset_hydraulic_state();
        } else if (command == "solve_hydraulics") {
            require_arguments(nrhs, 2, "status = staci_mex('solve_hydraulics', handle)");
            const StaciHydraulicStatus status = session(prhs[1]).solve_hydraulics();
            const char *fields[] = {"converged"};
            plhs[0] = mxCreateStructMatrix(1, 1, 1, fields);
            mxSetField(plhs[0], 0, fields[0], mxCreateLogicalScalar(status.converged));
        } else if (command == "get_node_property" ||
                   command == "get_link_property") {
            require_arguments(nrhs, 4,
                "values = staci_mex('get_node_property'|'get_link_property', handle, ids, property)");
            const std::vector<std::string> ids = text_vector(prhs[2], "ids");
            const std::string property = text_value(prhs[3], "property");
            const std::vector<double> values = command == "get_node_property"
                ? session(prhs[1]).get_node_property(ids, property)
                : session(prhs[1]).get_link_property(ids, property);
            plhs[0] = double_vector(values);
        } else if (command == "set_node_property" ||
                   command == "set_link_property") {
            require_arguments(nrhs, 5,
                "staci_mex('set_node_property'|'set_link_property', handle, ids, property, values)");
            const std::vector<std::string> ids = text_vector(prhs[2], "ids");
            const std::string property = text_value(prhs[3], "property");
            const std::vector<double> values = numeric_vector(prhs[4], "values");
            if (command == "set_node_property")
                session(prhs[1]).set_node_property(ids, property, values);
            else
                session(prhs[1]).set_link_property(ids, property, values);
        } else if (command == "node_results") {
            require_arguments(nrhs, 2, "results = staci_mex('node_results', handle)");
            plhs[0] = node_results_struct(session(prhs[1]).node_results());
        } else if (command == "link_results") {
            require_arguments(nrhs, 2, "results = staci_mex('link_results', handle)");
            plhs[0] = link_results_struct(session(prhs[1]).link_results());
        } else if (command == "solve_steady_quality") {
            require_arguments(nrhs, 3,
                "results = staci_mex('solve_steady_quality', handle, mode)");
            const std::string mode = text_value(prhs[2], "mode");
            const bool age = mode == "age" || mode == "both";
            const bool chemical = mode == "chemical" || mode == "chlorine" ||
                                  mode == "both";
            if (!age && !chemical)
                throw std::invalid_argument(
                    "Quality mode must be 'age', 'chemical', 'chlorine' or 'both'.");
            plhs[0] = quality_results_struct(
                session(prhs[1]).solve_steady_quality(age, chemical));
        } else if (command == "hydraulic_sensitivity") {
            require_arguments(nrhs, 4,
                "result = staci_mex('hydraulic_sensitivity', handle, elementId, property)");
            plhs[0] = sensitivity_results_struct(
                session(prhs[1]).hydraulic_sensitivity(
                    text_value(prhs[2], "elementId"),
                    text_value(prhs[3], "property")));
        } else if (command == "version") {
            require_arguments(nrhs, 1, "version = staci_mex('version')");
            plhs[0] = mxCreateString("STACI MATLAB API 1.0");
        } else {
            throw std::invalid_argument("Unknown staci_mex command: " + command);
        }
    } catch (const std::exception &error) {
        mexErrMsgIdAndTxt("STACI:MatlabAPI", "%s", error.what());
    }
}
