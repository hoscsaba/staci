if(NOT EXISTS "${NODES_CSV}")
    message(FATAL_ERROR "Missing pump EPS node results: ${NODES_CSV}")
endif()

file(READ "${NODES_CSV}" node_results)
foreach(expected
        "0,J1,0,539."
        "3600,J1,0,290."
        "7200,J1,0,93."
        "0,J2,0,87.5"
        "3600,J2,0,69.9"
        "7200,J2,0,49.5")
    string(FIND "${node_results}" "${expected}" position)
    if(position EQUAL -1)
        message(FATAL_ERROR
            "Pump speed-pattern result does not contain expected value: ${expected}")
    endif()
endforeach()
