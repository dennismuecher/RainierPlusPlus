cmake_minimum_required(VERSION 3.15)
project(RAINIER VERSION 2.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS Core Hist Tree RIO Graf MathCore)
include(${ROOT_USE_FILE})

# OpenMP support (optional)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP found - parallel execution enabled")
endif()

# Include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${ROOT_INCLUDE_DIRS}
)

# Source files - NOW INCLUDING PRIORITY 3 FILES
set(SOURCES
    src/Config.cpp
    src/core/Level.cpp
    src/core/Nucleus.cpp
    src/io/LevelFileReader.cpp
    src/models/LevelDensity.cpp
    src/models/SpinCutoff.cpp
    src/models/GammaStrength.cpp
    src/simulation/CascadeEvent.cpp
    src/simulation/DecaySimulator.cpp
    src/io/OutputManager.cpp
)

# Header files
set(HEADERS
    include/Config.h
    include/core/Level.h
    include/core/DiscreteLevel.h
    include/core/ContinuumLevel.h
    include/core/Transition.h
    include/core/Nucleus.h
    include/models/LevelDensity.h
    include/models/SpinCutoff.h
    include/models/GammaStrength.h
    include/simulation/DecaySimulator.h
    include/io/LevelFileReader.h
    include/io/OutputManager.h
    include/utils/PhysicsConstants.h
)

# Create library
add_library(rainier_lib STATIC ${SOURCES} ${HEADERS})
target_link_libraries(rainier_lib ${ROOT_LIBRARIES})

# Main executable
add_executable(rainier src/main.cpp)
target_link_libraries(rainier rainier_lib ${ROOT_LIBRARIES})

# Compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(rainier_lib PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Install
install(TARGETS rainier DESTINATION bin)
install(TARGETS rainier_lib DESTINATION lib)
install(DIRECTORY include/ DESTINATION include/rainier)

# Print configuration summary
message(STATUS "")
message(STATUS "RAINIER 2.0 Configuration:")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  ROOT version: ${ROOT_VERSION}")
message(STATUS "  Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")
