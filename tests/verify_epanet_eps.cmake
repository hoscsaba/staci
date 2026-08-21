if(NOT DEFINED RESULT_PREFIX)
    message(FATAL_ERROR "RESULT_PREFIX was not supplied")
endif()

foreach(kind nodes links tanks summary)
    set(path "${RESULT_PREFIX}-${kind}.csv")
    if(NOT EXISTS "${path}")
        message(FATAL_ERROR "Missing EPS output: ${path}")
    endif()
endforeach()

file(READ "${RESULT_PREFIX}-nodes.csv" nodes)
if(NOT nodes MATCHES "time_seconds,node_id,elevation_m,pressure_head_m,total_head_m,demand_m3s,converged")
    message(FATAL_ERROR "The node CSV schema is not SI-only")
endif()
if(NOT nodes MATCHES "0,J1,0,[^\n]*,0\\.001,1")
    message(FATAL_ERROR "The time-zero patterned demand was not written correctly")
endif()
if(NOT nodes MATCHES "3600,J1,0,[^\n]*,0\\.002,1")
    message(FATAL_ERROR "The second pattern multiplier was not applied")
endif()

file(READ "${RESULT_PREFIX}-tanks.csv" tanks)
if(NOT tanks MATCHES "0,T1,10,")
    message(FATAL_ERROR "The initial tank level is missing")
endif()

file(READ "${RESULT_PREFIX}-links.csv" links)
if(NOT links MATCHES "0,Pump1,PUMP,[^\n]*,0,1")
    message(FATAL_ERROR "The tank-level pump control was not applied")
endif()
if(NOT tanks MATCHES "7200,T1,")
    message(FATAL_ERROR "The final tank state is missing")
endif()

file(READ "${RESULT_PREFIX}-summary.csv" summary)
if(NOT summary MATCHES "hydraulic_states,3")
    message(FATAL_ERROR "Unexpected EPS state count")
endif()
if(NOT summary MATCHES "failed_states,0")
    message(FATAL_ERROR "The EPS smoke test did not converge")
endif()

if(NOT EXISTS "${RESULT_PREFIX}.meta.json")
    message(FATAL_ERROR "Missing EPS metadata JSON")
endif()
file(READ "${RESULT_PREFIX}.meta.json" metadata)
if(NOT metadata MATCHES "\"units\": \"SI\"")
    message(FATAL_ERROR "EPS metadata does not declare SI units")
endif()
if(NOT metadata MATCHES "\"frames\": 3")
    message(FATAL_ERROR "Unexpected number of reported EPS frames")
endif()
if(NOT metadata MATCHES "\"demand\": \{\"unit\": \"m3/s\"")
    message(FATAL_ERROR "EPS demand metadata is not in SI units")
endif()

if(EXPECT_HDF5)
    if(NOT EXISTS "${RESULT_PREFIX}.h5")
        message(FATAL_ERROR "Missing chunked EPS HDF5 output")
    endif()
    file(READ "${RESULT_PREFIX}.h5" signature OFFSET 0 LIMIT 8 HEX)
    string(TOLOWER "${signature}" signature)
    if(NOT signature STREQUAL "894844460d0a1a0a")
        message(FATAL_ERROR "EPS result does not have a valid HDF5 signature")
    endif()
    if(H5DUMP_EXECUTABLE AND EXISTS "${H5DUMP_EXECUTABLE}")
        execute_process(
            COMMAND "${H5DUMP_EXECUTABLE}" -n "${RESULT_PREFIX}.h5"
            RESULT_VARIABLE h5dump_status
            OUTPUT_VARIABLE h5_names
            ERROR_VARIABLE h5dump_error
        )
        if(NOT h5dump_status EQUAL 0)
            message(FATAL_ERROR "h5dump failed: ${h5dump_error}")
        endif()
        foreach(dataset
                "/time"
                "/nodes/head"
                "/nodes/pressure_head"
                "/nodes/demand"
                "/links/flow_rate"
                "/links/velocity"
                "/links/headloss"
                "/links/status"
                "/tanks/level"
                "/tanks/volume"
                "/tanks/inflow"
                "/simulation/converged")
            if(NOT h5_names MATCHES "dataset    ${dataset}")
                message(FATAL_ERROR "Missing HDF5 dataset: ${dataset}")
            endif()
        endforeach()
    endif()
endif()
