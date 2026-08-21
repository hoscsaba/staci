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

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_EXTENSIONS OFF
)

target_include_directories(${PROJECT_NAME} PRIVATE
    "${PROJECT_SOURCE_DIR}/include/staci"
)

# M_PI is available by default with GCC and Clang, but MSVC requires this
# definition before including <cmath> or <math.h>.
if(MSVC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE _USE_MATH_DEFINES)
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
        target_compile_definitions(${PROJECT_NAME} PRIVATE STACI_HAS_HDF5=1)
        if(TARGET HDF5::HDF5)
            target_link_libraries(${PROJECT_NAME} PRIVATE HDF5::HDF5)
        else()
            target_include_directories(${PROJECT_NAME} PRIVATE ${HDF5_INCLUDE_DIRS})
            target_link_libraries(${PROJECT_NAME} PRIVATE ${HDF5_LIBRARIES})
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
    target_link_libraries(${PROJECT_NAME} PRIVATE SuiteSparse::UMFPACK)
elseif(TARGET SuiteSparse::UMFPACK_static)
    target_link_libraries(${PROJECT_NAME} PRIVATE SuiteSparse::UMFPACK_static)
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

    target_include_directories(${PROJECT_NAME} PRIVATE
        "${SUITESPARSE_INCLUDE_DIR}"
    )
    target_link_libraries(${PROJECT_NAME} PRIVATE
        "${SUITESPARSE_UMFPACK_LIBRARY}"
        "${SUITESPARSE_CHOLMOD_LIBRARY}"
        "${SUITESPARSE_AMD_LIBRARY}"
        "${SUITESPARSE_COLAMD_LIBRARY}"
        "${SUITESPARSE_CONFIG_LIBRARY}"
    )
endif()

include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

include(CTest)
if(BUILD_TESTING)
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
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
