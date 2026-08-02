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

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
