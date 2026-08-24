if(NOT DEFINED RESULT_PREFIX)
    message(FATAL_ERROR "RESULT_PREFIX was not supplied")
endif()

foreach(kind nodes links tanks summary)
    if(NOT EXISTS "${RESULT_PREFIX}-${kind}.csv")
        message(FATAL_ERROR "Missing rule-based EPS output: ${RESULT_PREFIX}-${kind}.csv")
    endif()
endforeach()

file(READ "${RESULT_PREFIX}-summary.csv" summary)
if(NOT summary MATCHES "rules,15")
    message(FATAL_ERROR "Not all EPANET rules were parsed")
endif()
if(NOT summary MATCHES "rule_timestep_seconds,600")
    message(FATAL_ERROR "RULE TIMESTEP was not parsed")
endif()
if(NOT summary MATCHES "hydraulic_states,10")
    message(FATAL_ERROR "Rule events did not shorten the hydraulic timestep")
endif()
if(NOT summary MATCHES "failed_states,0")
    message(FATAL_ERROR "The rule-based EPS test did not converge")
endif()

file(READ "${RESULT_PREFIX}-links.csv" links)
foreach(link TimePipe ClockPipe PressurePipe FlowPipe EqualPipe TimeEqualPipe
             TankRulePipe FillTimePipe)
    if(NOT links MATCHES "1800,${link},PIPE,[^\n]*,0,1")
        message(FATAL_ERROR "Expected rule did not close ${link}")
    endif()
endforeach()

# T1 is filling throughout this fixture. EPANET treats DRAINTIME as an
# inapplicable/false premise unless a tank is actually draining, so this link
# must remain open. The official EPANET reference comparison checks the same
# status directly.
if(NOT links MATCHES "1800,DrainTimePipe,PIPE,[^\n]*,1,1")
    message(FATAL_ERROR "DRAINTIME rule fired while the tank was filling")
endif()

if(NOT links MATCHES "1800,ClockEqualPipe,PIPE,[^\n]*,1,1" OR
   NOT links MATCHES "3600,ClockEqualPipe,PIPE,[^\n]*,0,1")
    message(FATAL_ERROR "CLOCKTIME equality was not detected inside a rule interval")
endif()
if(NOT links MATCHES "1800,LogicPipe,PIPE,[^\n]*,1,1" OR
   NOT links MATCHES "3600,LogicPipe,PIPE,[^\n]*,0,1")
    message(FATAL_ERROR "EPANET OR-before-AND rule precedence is incorrect")
endif()
if(NOT links MATCHES "1800,ConflictPipe,PIPE,[^\n]*,1,1")
    message(FATAL_ERROR "Higher-priority rule did not win")
endif()
if(NOT links MATCHES "1800,EqualPipe,PIPE,[^\n]*,0,1")
    message(FATAL_ERROR "First equal-priority rule did not win")
endif()
if(NOT links MATCHES "1800,PSet,PUMP,[^\n]*,0\\.071[0-9]*,nan,[^\n]*,1,1")
    message(FATAL_ERROR "Numeric pump rule setting was not applied")
endif()
if(NOT links MATCHES "1800,ClockPipe,PIPE,[^\n]*,0,1" OR
   NOT links MATCHES "5400,ClockPipe,PIPE,[^\n]*,1,1")
    message(FATAL_ERROR "THEN/ELSE clock-window actions are incorrect")
endif()
