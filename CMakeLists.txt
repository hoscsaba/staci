cmake_minimum_required(VERSION 3.16)
project(staci VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include headers
include_directories(${PROJECT_SOURCE_DIR}/include/staci)
include_directories("/opt/homebrew/include/suitesparse")
link_directories("/opt/homebrew/lib")

# Collect all .cpp files
file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
)

# Separate main file (Staci.cpp)
set(MAIN_FILE "${PROJECT_SOURCE_DIR}/src/Staci.cpp")
set(EXCLUDE_FILES
        "${PROJECT_SOURCE_DIR}/src/staci_calibrate.cpp"
        "${PROJECT_SOURCE_DIR}/src/staci_split.cpp"
)

# Remove excluded files from sources
list(REMOVE_ITEM SRC_FILES ${MAIN_FILE} ${EXCLUDE_FILES})

# Show final list of sources to compile
message(STATUS "Final sources: ${SRC_FILES}")

# Add executable
add_executable(${PROJECT_NAME} ${MAIN_FILE} ${SRC_FILES})

# Add the SuiteSparse include path
include_directories("/opt/homebrew/include/suitesparse")

# Link SuiteSparse libraries
target_link_libraries(${PROJECT_NAME}
        umfpack
        cholmod
        amd
        colamd
        suitesparseconfig
)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")