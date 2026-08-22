cmake_minimum_required(VERSION 3.16)

project(staci VERSION 1.0 LANGUAGES C CXX)

# Collect all .cpp files.
file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
)

# Staci.cpp is the main implementation file. These utility programs have their
# own entry points and are not part of the staci executable.
set(MAIN_FILE "${PROJECT_SOURCE_DIR}/src/Staci.cpp")
set(EXCLUDE_FILES
    "${PROJECT_SOURCE_DIR}/src/staci_calibrate.cpp"
    "${PROJECT_SOURCE_DIR}/src/staci_split.cpp"
    "${PROJECT_SOURCE_DIR}/src/data_io_2.cpp"
)
list(REMOVE_ITEM SRC_FILES ${MAIN_FILE} ${EXCLUDE_FILES})

add_executable(${PROJECT_NAME} ${MAIN_FILE} ${SRC_FILES})

option(STACI_BUILD_OPTIMIZERS "Build the pagmo2-based calibrate and split tools" ON)
set(STACI_EXECUTABLE_TARGETS ${PROJECT_NAME})

if(STACI_BUILD_OPTIMIZERS)
    find_package(Pagmo CONFIG REQUIRED)
    find_package(Eigen3 CONFIG REQUIRED)
    find_package(igraph CONFIG REQUIRED)

    set(STACI_OPTIMIZER_COMMON_SOURCES ${SRC_FILES} ${MAIN_FILE})
    list(REMOVE_ITEM STACI_OPTIMIZER_COMMON_SOURCES
        "${PROJECT_SOURCE_DIR}/src/staci_main.cpp"
    )

    add_executable(staci_calibrate
        "${PROJECT_SOURCE_DIR}/src/staci_calibrate.cpp"
        ${STACI_OPTIMIZER_COMMON_SOURCES}
    )
    target_link_libraries(staci_calibrate PRIVATE Pagmo::pagmo)

    add_executable(staci_split
        "${PROJECT_SOURCE_DIR}/src/staci_split.cpp"
        ${STACI_OPTIMIZER_COMMON_SOURCES}
    )
    target_link_libraries(staci_split PRIVATE Pagmo::pagmo Eigen3::Eigen igraph::igraph)
    list(APPEND STACI_EXECUTABLE_TARGETS staci_calibrate staci_split)
endif()

foreach(staci_target IN LISTS STACI_EXECUTABLE_TARGETS)
    target_compile_features(${staci_target} PRIVATE cxx_std_17)
    set_target_properties(${staci_target} PROPERTIES CXX_EXTENSIONS OFF)
    target_include_directories(${staci_target} PRIVATE
        "${PROJECT_SOURCE_DIR}/include/staci"
    )
endforeach()

# M_PI is available by default with GCC and Clang, but MSVC requires this
# definition before including <cmath> or <math.h>.
if(MSVC)
    foreach(staci_target IN LISTS STACI_EXECUTABLE_TARGETS)
        target_compile_definitions(${staci_target} PRIVATE _USE_MATH_DEFINES)
    endforeach()
endif()

# HDF5 is optional at build time so the legacy executable remains portable on
# minimal systems. When present, EPS runs write the chunked STACI EPS OUTPUT v1
# file in addition to the SI CSV and JSON outputs.
option(STACI_ENABLE_HDF5 "Enable chunked HDF5 output for EPANET EPS runs" ON)
set(STACI_HAS_HDF5 OFF)
if(STACI_ENABLE_HDF5)
    find_package(HDF5 COMPONENTS C QUIET)
    if(HDF5_FOUND)
        set(STACI_HAS_HDF5 ON)
        foreach(staci_target IN LISTS STACI_EXECUTABLE_TARGETS)
            target_compile_definitions(${staci_target} PRIVATE STACI_HAS_HDF5=1)
        endforeach()
        if(TARGET HDF5::HDF5)
            foreach(staci_target IN LISTS STACI_EXECUTABLE_TARGETS)
                target_link_libraries(${staci_target} PRIVATE HDF5::HDF5)
            endforeach()
        else()
            foreach(staci_target IN LISTS STACI_EXECUTABLE_TARGETS)
                target_include_directories(${staci_target} PRIVATE ${HDF5_INCLUDE_DIRS})
                target_link_libraries(${staci_target} PRIVATE ${HDF5_LIBRARIES})
            endforeach()
        endif()
        find_program(STACI_H5DUMP_EXECUTABLE NAMES h5dump
            HINTS "${HDF5_ROOT}/bin")
        message(STATUS "HDF5 EPS output: enabled (${HDF5_VERSION})")
    else()
        message(WARNING
            "HDF5 was not found. STACI will build, but EPS runs will only write SI CSV and JSON files. "
            "Install the HDF5 development package or set HDF5_ROOT to enable .h5 output.")
    endif()
endif()

# Modern SuiteSparse installations provide an imported UMFPACK target that
# carries all required include paths and transitive dependencies. This is the
# preferred path on Windows (for example with vcpkg), Linux, and macOS.
find_package(UMFPACK CONFIG QUIET)

if(TARGET SuiteSparse::UMFPACK)
    foreach(staci_target IN LISTS STACI_EXECUTABLE_TARGETS)
        target_link_libraries(${staci_target} PRIVATE SuiteSparse::UMFPACK)
    endforeach()
elseif(TARGET SuiteSparse::UMFPACK_static)
    foreach(staci_target IN LISTS STACI_EXECUTABLE_TARGETS)
        target_link_libraries(${staci_target} PRIVATE SuiteSparse::UMFPACK_static)
    endforeach()
else()
    # Some Linux distributions ship SuiteSparse without CMake package files.
    # Fall back to normal CMake searches; callers can use CMAKE_PREFIX_PATH or
    # the SUITESPARSE_ROOT environment variable for a non-standard install.
    find_path(SUITESPARSE_INCLUDE_DIR
        NAMES umfpack.h
        HINTS "$ENV{SUITESPARSE_ROOT}"
        PATH_SUFFIXES include/suitesparse include SuiteSparse
    )

    find_library(SUITESPARSE_UMFPACK_LIBRARY
        NAMES umfpack
        HINTS "$ENV{SUITESPARSE_ROOT}"
        PATH_SUFFIXES lib lib64
    )
    find_library(SUITESPARSE_CHOLMOD_LIBRARY
        NAMES cholmod
        HINTS "$ENV{SUITESPARSE_ROOT}"
        PATH_SUFFIXES lib lib64
    )
    find_library(SUITESPARSE_AMD_LIBRARY
        NAMES amd
        HINTS "$ENV{SUITESPARSE_ROOT}"
        PATH_SUFFIXES lib lib64
    )
    find_library(SUITESPARSE_COLAMD_LIBRARY
        NAMES colamd
        HINTS "$ENV{SUITESPARSE_ROOT}"
        PATH_SUFFIXES lib lib64
    )
    find_library(SUITESPARSE_CONFIG_LIBRARY
        NAMES suitesparseconfig
        HINTS "$ENV{SUITESPARSE_ROOT}"
        PATH_SUFFIXES lib lib64
    )

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(SuiteSparse REQUIRED_VARS
        SUITESPARSE_INCLUDE_DIR
        SUITESPARSE_UMFPACK_LIBRARY
        SUITESPARSE_CHOLMOD_LIBRARY
        SUITESPARSE_AMD_LIBRARY
        SUITESPARSE_COLAMD_LIBRARY
        SUITESPARSE_CONFIG_LIBRARY
    )

    foreach(staci_target IN LISTS STACI_EXECUTABLE_TARGETS)
        target_include_directories(${staci_target} PRIVATE "${SUITESPARSE_INCLUDE_DIR}")
        target_link_libraries(${staci_target} PRIVATE
            "${SUITESPARSE_UMFPACK_LIBRARY}"
            "${SUITESPARSE_CHOLMOD_LIBRARY}"
            "${SUITESPARSE_AMD_LIBRARY}"
            "${SUITESPARSE_COLAMD_LIBRARY}"
            "${SUITESPARSE_CONFIG_LIBRARY}"
        )
    endforeach()
endif()

include(GNUInstallDirs)
install(TARGETS ${STACI_EXECUTABLE_TARGETS}
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

include(CTest)
if(BUILD_TESTING)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        add_test(NAME channel_reference_suite
            COMMAND "${Python3_EXECUTABLE}"
                    "${PROJECT_SOURCE_DIR}/tests/run_channel_tests.py"
                    --binary $<TARGET_FILE:${PROJECT_NAME}>
                    --results "${CMAKE_CURRENT_BINARY_DIR}/channel-reference-tests"
        )
        set_tests_properties(channel_reference_suite PROPERTIES
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            TIMEOUT 180
            LABELS "channel;hydraulics;reference"
        )
        add_test(NAME channel_network_merge_split
            COMMAND "${Python3_EXECUTABLE}"
                    "${PROJECT_SOURCE_DIR}/tests/run_channel_network_test.py"
                    --binary $<TARGET_FILE:${PROJECT_NAME}>
                    --network "${PROJECT_SOURCE_DIR}/tests/channel_network_merge_split.spr"
                    --results "${CMAKE_CURRENT_BINARY_DIR}/channel-network-test"
        )
        set_tests_properties(channel_network_merge_split PROPERTIES
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            TIMEOUT 120
            LABELS "channel;hydraulics;network;reference"
        )
    endif()

    if(STACI_BUILD_OPTIMIZERS)
        set(STACI_CALIBRATE_TEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/optimizer-tests/calibrate")
        set(STACI_SPLIT_TEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/optimizer-tests/split")
        file(MAKE_DIRECTORY "${STACI_CALIBRATE_TEST_DIR}" "${STACI_SPLIT_TEST_DIR}")
        configure_file(
            "${PROJECT_SOURCE_DIR}/tests/staci_calibrate_settings.xml.in"
            "${STACI_CALIBRATE_TEST_DIR}/staci_calibrate_settings.xml"
            @ONLY
        )
        configure_file(
            "${PROJECT_SOURCE_DIR}/tests/calibrate_targets.csv"
            "${STACI_CALIBRATE_TEST_DIR}/calibrate_targets.csv"
            COPYONLY
        )
        configure_file(
            "${PROJECT_SOURCE_DIR}/tests/anytown_1med.spr"
            "${STACI_CALIBRATE_TEST_DIR}/calibrate_network_0.spr"
            COPYONLY
        )
        configure_file(
            "${PROJECT_SOURCE_DIR}/tests/staci_split_settings.xml.in"
            "${STACI_SPLIT_TEST_DIR}/staci_split_settings.xml"
            @ONLY
        )
        configure_file(
            "${PROJECT_SOURCE_DIR}/tests/LOV-LOVOTV-2-input_mod.spr"
            "${STACI_SPLIT_TEST_DIR}/split-network.spr"
            COPYONLY
        )

        add_test(NAME pagmo_calibrate_small_network
            COMMAND $<TARGET_FILE:staci_calibrate> --seed 12345
        )
        set_tests_properties(pagmo_calibrate_small_network PROPERTIES
            WORKING_DIRECTORY "${STACI_CALIBRATE_TEST_DIR}"
            PASS_REGULAR_EXPRESSION "Optimization finished"
            TIMEOUT 120
        )

        add_test(NAME pagmo_split_large_network
            COMMAND $<TARGET_FILE:staci_split> --seed 12345
        )
        set_tests_properties(pagmo_split_large_network PROPERTIES
            WORKING_DIRECTORY "${STACI_SPLIT_TEST_DIR}"
            PASS_REGULAR_EXPRESSION "# of pipes : 2479"
            TIMEOUT 300
        )
    endif()

    add_executable(staci_epanet_document_test
        "${PROJECT_SOURCE_DIR}/tests/test_epanet_document.cpp"
        "${PROJECT_SOURCE_DIR}/src/epanet_document.cpp"
    )
    target_compile_features(staci_epanet_document_test PRIVATE cxx_std_17)
    set_target_properties(staci_epanet_document_test PROPERTIES CXX_EXTENSIONS OFF)
    target_include_directories(staci_epanet_document_test PRIVATE
        "${PROJECT_SOURCE_DIR}/include/staci"
    )

    add_test(NAME epanet_document_associations
        COMMAND $<TARGET_FILE:staci_epanet_document_test>
                "${PROJECT_SOURCE_DIR}/tests/epanet_metadata_roundtrip.inp"
    )

    add_executable(staci_epanet_node_metadata_test
        "${PROJECT_SOURCE_DIR}/tests/test_epanet_node_metadata.cpp"
        "${PROJECT_SOURCE_DIR}/src/Csomopont.cpp"
        "${PROJECT_SOURCE_DIR}/src/Agelem.cpp"
        "${PROJECT_SOURCE_DIR}/src/Cso.cpp"
        "${PROJECT_SOURCE_DIR}/src/EpanetPowerPump.cpp"
        "${PROJECT_SOURCE_DIR}/src/JelleggorbesFojtas.cpp"
        "${PROJECT_SOURCE_DIR}/src/KonstNyomas.cpp"
        "${PROJECT_SOURCE_DIR}/src/Szivattyu.cpp"
        "${PROJECT_SOURCE_DIR}/src/lubksb.cpp"
        "${PROJECT_SOURCE_DIR}/src/ludcmp.cpp"
        "${PROJECT_SOURCE_DIR}/src/epanet_document.cpp"
        "${PROJECT_SOURCE_DIR}/src/epanet_writer.cpp"
    )
    target_compile_features(staci_epanet_node_metadata_test PRIVATE cxx_std_17)
    set_target_properties(staci_epanet_node_metadata_test PROPERTIES CXX_EXTENSIONS OFF)
    target_include_directories(staci_epanet_node_metadata_test PRIVATE
        "${PROJECT_SOURCE_DIR}/include/staci"
    )

    add_test(NAME epanet_node_metadata_storage_and_export
        COMMAND $<TARGET_FILE:staci_epanet_node_metadata_test>
                "${CMAKE_CURRENT_BINARY_DIR}/epanet-node-metadata-native.inp"
    )
    set_tests_properties(epanet_node_metadata_storage_and_export PROPERTIES
        FAIL_REGULAR_EXPRESSION "TCV-1.*omitted;TCV-1.*no lossless EPANET mapping"
    )

    add_test(NAME epanet_demand_metadata_import
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -g
                "${PROJECT_SOURCE_DIR}/tests/epanet_node_metadata.inp"
                -e J1 -p demand
    )
    set_tests_properties(epanet_demand_metadata_import PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "preserved demand components: 3"
    )

    add_test(NAME epanet_initial_quality_metadata_import
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -g
                "${PROJECT_SOURCE_DIR}/tests/epanet_node_metadata.inp"
                -e J1 -p demand
    )
    set_tests_properties(epanet_initial_quality_metadata_import PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "preserved initial quality values: 2"
    )

    add_test(NAME epanet_demand_metadata_snapshot
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -g
                "${PROJECT_SOURCE_DIR}/tests/epanet_node_metadata.inp"
                -e J1 -p demand
    )
    set_tests_properties(epanet_demand_metadata_snapshot PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "value       : 41.04"
    )

    add_test(NAME epanet_reservoir_head_pattern_metadata
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -g
                "${PROJECT_SOURCE_DIR}/tests/epanet_eps_smoke.inp"
                -e EPANET_RESERVOIR_R1 -p head
    )

    add_test(NAME epanet_pipe_minor_loss_import
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -g
                "${PROJECT_SOURCE_DIR}/tests/epanet_node_metadata.inp"
                -e P1 -p minor_loss
    )
    set_tests_properties(epanet_pipe_minor_loss_import PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "value       : 2.5"
    )

    add_test(NAME epanet_pipe_status_override
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -g
                "${PROJECT_SOURCE_DIR}/tests/epanet_node_metadata.inp"
                -e P1 -p status
    )
    set_tests_properties(epanet_pipe_status_override PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "value       : 0"
    )

    add_test(NAME epanet_pipe_check_valve_import
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -g
                "${PROJECT_SOURCE_DIR}/tests/epanet_node_metadata.inp"
                -e P2 -p check_valve
    )
    set_tests_properties(epanet_pipe_check_valve_import PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "value       : 1"
    )

    add_test(NAME epanet_tcv_setting_import
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -g
                "${PROJECT_SOURCE_DIR}/tests/epanet_tcv.inp"
                -e TCV1 -p tcv_setting
    )
    set_tests_properties(epanet_tcv_setting_import PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "value       : 2"
        FAIL_REGULAR_EXPRESSION "TCV was not imported;EPANET valve type.*has no equivalent"
    )

    add_test(NAME epanet_tcv_position_import
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -g
                "${PROJECT_SOURCE_DIR}/tests/epanet_tcv.inp"
                -e TCV1 -p position
    )
    set_tests_properties(epanet_tcv_position_import PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "value       : 50"
    )

    add_test(NAME epanet_tcv_minor_loss_import
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -g
                "${PROJECT_SOURCE_DIR}/tests/epanet_tcv.inp"
                -e TCV1 -p tcv_minor_loss
    )
    set_tests_properties(epanet_tcv_minor_loss_import PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "value       : 0.5"
    )

    add_test(NAME epanet_modify_tcv_setting
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -m
                "${PROJECT_SOURCE_DIR}/tests/epanet_tcv.inp"
                -e TCV1 -p tcv_setting -n 12
                -o "${CMAKE_CURRENT_BINARY_DIR}/epanet-tcv-modified.inp"
    )
    set_tests_properties(epanet_modify_tcv_setting PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "EPANET modification: wrote"
    )

    add_test(NAME epanet_modify_tcv_setting_roundtrip
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -g
                "${CMAKE_CURRENT_BINARY_DIR}/epanet-tcv-modified.inp"
                -e TCV1 -p tcv_setting
    )
    set_tests_properties(epanet_modify_tcv_setting_roundtrip PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        DEPENDS epanet_modify_tcv_setting
        PASS_REGULAR_EXPRESSION "value       : 12"
    )

    add_test(NAME epanet_modify_pipe_minor_loss
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -m
                "${PROJECT_SOURCE_DIR}/tests/epanet_node_metadata.inp"
                -e P1 -p minor_loss -n 3.75
                -o "${CMAKE_CURRENT_BINARY_DIR}/epanet-minor-loss-modified.inp"
    )
    set_tests_properties(epanet_modify_pipe_minor_loss PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "EPANET modification: wrote"
    )

    add_test(NAME epanet_modify_pipe_minor_loss_roundtrip
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -g
                "${CMAKE_CURRENT_BINARY_DIR}/epanet-minor-loss-modified.inp"
                -e P1 -p minor_loss
    )
    set_tests_properties(epanet_modify_pipe_minor_loss_roundtrip PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        DEPENDS epanet_modify_pipe_minor_loss
        PASS_REGULAR_EXPRESSION "value       : 3.75"
    )
    set_tests_properties(epanet_reservoir_head_pattern_metadata PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "preserved reservoir head patterns: 1"
    )

    add_test(NAME epanet_power_pump_speed_import
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -g
                "${PROJECT_SOURCE_DIR}/tests/epanet_pumps.inp"
                -e PWR -p base_speed
    )
    set_tests_properties(epanet_power_pump_speed_import PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "value       : 0.8"
    )

    add_test(NAME epanet_head_pump_metadata_import
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -g
                "${PROJECT_SOURCE_DIR}/tests/epanet_pumps.inp"
                -e PHEAD -p speed_pattern_length
    )
    set_tests_properties(epanet_head_pump_metadata_import PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "value       : 3"
    )

    add_test(NAME epanet_one_point_head_pump_import
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -g
                "${PROJECT_SOURCE_DIR}/tests/epanet_pumps.inp"
                -e PONE -p head_curve_points
    )
    set_tests_properties(epanet_one_point_head_pump_import PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "value       : 1"
    )

    add_test(NAME epanet_numeric_pump_status_import
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -g
                "${PROJECT_SOURCE_DIR}/tests/epanet_pumps.inp"
                -e PONE -p speed
    )
    set_tests_properties(epanet_numeric_pump_status_import PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "value       : 0.9"
    )

    add_test(NAME epanet_power_pump_efficiency_import
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -g
                "${PROJECT_SOURCE_DIR}/tests/epanet_pumps.inp"
                -e PWR -p efficiency_curve_points
    )
    set_tests_properties(epanet_power_pump_efficiency_import PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "value       : 3"
    )

    set(STACI_PUMP_EPS_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/epanet-pump-eps")
    add_test(NAME epanet_pump_speed_pattern_eps
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -z
                "${PROJECT_SOURCE_DIR}/tests/epanet_pumps.inp"
                -o "${STACI_PUMP_EPS_PREFIX}"
    )
    set_tests_properties(epanet_pump_speed_pattern_eps PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "3 hydraulic states, 0 failed states"
    )

    add_test(NAME epanet_pump_speed_pattern_results
        COMMAND "${CMAKE_COMMAND}"
                -DNODES_CSV=${STACI_PUMP_EPS_PREFIX}-nodes.csv
                -P "${PROJECT_SOURCE_DIR}/tests/verify_epanet_pumps.cmake"
    )
    set_tests_properties(epanet_pump_speed_pattern_results PROPERTIES
        DEPENDS epanet_pump_speed_pattern_eps
    )

    add_test(NAME epanet_inp_import
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -l
                "${PROJECT_SOURCE_DIR}/tests/epanet_smoke.inp"
    )
    set_tests_properties(epanet_inp_import PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "imported STACI elements: 3"
    )

    add_test(NAME epanet_export_spr
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> --export-epanet
                "${PROJECT_SOURCE_DIR}/tests/anytown_1med.spr"
                -o "${CMAKE_CURRENT_BINARY_DIR}/anytown-export.inp"
    )
    set_tests_properties(epanet_export_spr PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "EPANET export: wrote"
    )

    add_test(NAME epanet_export_roundtrip
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -l
                "${CMAKE_CURRENT_BINARY_DIR}/anytown-export.inp"
    )
    set_tests_properties(epanet_export_roundtrip PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        DEPENDS epanet_export_spr
        PASS_REGULAR_EXPRESSION "imported STACI elements: 44"
    )

    add_test(NAME epanet_metadata_export
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> --export-epanet
                "${PROJECT_SOURCE_DIR}/tests/epanet_metadata_roundtrip.inp"
                -o "${CMAKE_CURRENT_BINARY_DIR}/epanet-metadata-roundtrip.inp"
    )
    set_tests_properties(epanet_metadata_export PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )

    add_test(NAME epanet_metadata_roundtrip
        COMMAND "${CMAKE_COMMAND}" -E compare_files
                "${PROJECT_SOURCE_DIR}/tests/epanet_metadata_roundtrip.inp"
                "${CMAKE_CURRENT_BINARY_DIR}/epanet-metadata-roundtrip.inp"
    )
    set_tests_properties(epanet_metadata_roundtrip PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        DEPENDS epanet_metadata_export
    )

    add_test(NAME epanet_modify_inp
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -m
                "${PROJECT_SOURCE_DIR}/tests/epanet_smoke.inp"
                -e Pipe1 -p diameter -n 0.2
                -o "${CMAKE_CURRENT_BINARY_DIR}/epanet-modified.inp"
    )
    set_tests_properties(epanet_modify_inp PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "EPANET modification: wrote"
    )

    add_test(NAME epanet_modify_roundtrip
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -g
                "${CMAKE_CURRENT_BINARY_DIR}/epanet-modified.inp"
                -e Pipe1 -p diameter
    )
    set_tests_properties(epanet_modify_roundtrip PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        DEPENDS epanet_modify_inp
        PASS_REGULAR_EXPRESSION "value       : 0.2"
    )

    add_test(NAME epanet_modify_patterned_demand
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -m
                "${PROJECT_SOURCE_DIR}/tests/epanet_smoke.inp"
                -e J1 -p demand -n 7.2
                -o "${CMAKE_CURRENT_BINARY_DIR}/epanet-demand-modified.inp"
    )
    set_tests_properties(epanet_modify_patterned_demand PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "EPANET modification: wrote"
    )

    add_test(NAME epanet_modify_patterned_demand_roundtrip
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> -g
                "${CMAKE_CURRENT_BINARY_DIR}/epanet-demand-modified.inp"
                -e J1 -p demand
    )
    set_tests_properties(epanet_modify_patterned_demand_roundtrip PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        DEPENDS epanet_modify_patterned_demand
        PASS_REGULAR_EXPRESSION "value       : 7.2"
    )

    add_test(NAME epanet_extended_period_simulation
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> --epanet-eps
                "${PROJECT_SOURCE_DIR}/tests/epanet_eps_smoke.inp"
                -o "${CMAKE_CURRENT_BINARY_DIR}/epanet-eps-result"
    )
    set_tests_properties(epanet_extended_period_simulation PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "EPANET EPS complete: 3 hydraulic states, 0 failed states"
    )

    add_test(NAME epanet_extended_period_results
        COMMAND "${CMAKE_COMMAND}"
                -DRESULT_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/epanet-eps-result
                -DEXPECT_HDF5=${STACI_HAS_HDF5}
                -DH5DUMP_EXECUTABLE=${STACI_H5DUMP_EXECUTABLE}
                -P "${PROJECT_SOURCE_DIR}/tests/verify_epanet_eps.cmake"
    )
    set_tests_properties(epanet_extended_period_results PROPERTIES
        DEPENDS epanet_extended_period_simulation
    )

    add_test(NAME epanet_simple_controls_simulation
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> --epanet-eps
                "${PROJECT_SOURCE_DIR}/tests/epanet_controls.inp"
                -o "${CMAKE_CURRENT_BINARY_DIR}/epanet-controls-result"
    )
    set_tests_properties(epanet_simple_controls_simulation PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "EPANET EPS complete: 33 hydraulic states, 0 failed states"
    )

    add_test(NAME epanet_simple_controls_results
        COMMAND "${CMAKE_COMMAND}"
                -DRESULT_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/epanet-controls-result
                -P "${PROJECT_SOURCE_DIR}/tests/verify_epanet_controls.cmake"
    )
    set_tests_properties(epanet_simple_controls_results PROPERTIES
        DEPENDS epanet_simple_controls_simulation
    )

    add_test(NAME epanet_rule_controls_simulation
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> --epanet-eps
                "${PROJECT_SOURCE_DIR}/tests/epanet_rules.inp"
                -o "${CMAKE_CURRENT_BINARY_DIR}/epanet-rules-result"
    )
    set_tests_properties(epanet_rule_controls_simulation PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "EPANET EPS complete: 10 hydraulic states, 0 failed states"
    )

    add_test(NAME epanet_rule_controls_results
        COMMAND "${CMAKE_COMMAND}"
                -DRESULT_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/epanet-rules-result
                -P "${PROJECT_SOURCE_DIR}/tests/verify_epanet_rules.cmake"
    )
    set_tests_properties(epanet_rule_controls_results PROPERTIES
        DEPENDS epanet_rule_controls_simulation
    )

    add_test(NAME epanet_tcv_controls_simulation
        COMMAND $<TARGET_FILE:${PROJECT_NAME}> --epanet-eps
                "${PROJECT_SOURCE_DIR}/tests/epanet_tcv.inp"
                -o "${CMAKE_CURRENT_BINARY_DIR}/epanet-tcv-result"
    )
    set_tests_properties(epanet_tcv_controls_simulation PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        PASS_REGULAR_EXPRESSION "EPANET EPS complete: 4 hydraulic states, 0 failed states"
    )

    add_test(NAME epanet_tcv_controls_results
        COMMAND "${CMAKE_COMMAND}"
                -DRESULT_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/epanet-tcv-result
                -P "${PROJECT_SOURCE_DIR}/tests/verify_epanet_tcv.cmake"
    )
    set_tests_properties(epanet_tcv_controls_results PROPERTIES
        DEPENDS epanet_tcv_controls_simulation
    )
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
