if(NOT DEFINED RESULT_PREFIX)
    message(FATAL_ERROR "RESULT_PREFIX was not supplied")
endif()

foreach(kind nodes links tanks summary)
    if(NOT EXISTS "${RESULT_PREFIX}-${kind}.csv")
        message(FATAL_ERROR "Missing volume-curve EPS output: ${RESULT_PREFIX}-${kind}.csv")
    endif()
endforeach()

file(READ "${RESULT_PREFIX}-summary.csv" summary)
if(NOT summary MATCHES "failed_states,0")
    message(FATAL_ERROR "The volume-curve EPS test did not converge")
endif()

file(READ "${RESULT_PREFIX}-tanks.csv" tanks)
# The curve explicitly maps initial level 10 m to 50 m3. A cylindrical
# approximation using the nominal 1 m diameter would instead report 3.927 m3.
if(NOT tanks MATCHES "0,T1,10,50,")
    message(FATAL_ERROR "Initial tank volume was not evaluated from the EPANET volume curve")
endif()
