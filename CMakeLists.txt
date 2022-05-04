#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.1)
project(dpm++)


#----------------------------------------------------------------------------
# Find Geant4 package
#
find_package(Geant4 REQUIRED)


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})


#-------------------------------------------------------------------------------
# Set the CXX flags and default build type
#
if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
    if (NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
    endif ()
endif ()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")


#-------------------------------------------------------------------------------
# Set the headers and sources for the data GENeration an SIMulation parts
#
file(GLOB headers_GEN
        ${PROJECT_SOURCE_DIR}/DataInit/inc/*.hh
        ${PROJECT_SOURCE_DIR}/Utils/inc/*.hh
        )
file(GLOB cxx_sources_GEN
        ${PROJECT_SOURCE_DIR}/DataInit/src/*.cc
        ${PROJECT_SOURCE_DIR}/Utils/src/*.cc
        )

file(GLOB headers_SIM
        ${PROJECT_SOURCE_DIR}*.hh
        ${PROJECT_SOURCE_DIR}/Simulation/inc/*.hh
        )
file(GLOB cxx_sources_SIM
        ${PROJECT_SOURCE_DIR}/Simulation/src/*.cc
        )
file(GLOB validation_sources
        ${PROJECT_SOURCE_DIR}/modelValidationTests/src/*.cc
        )


#----------------------------------------------------------------------------
# Locate headers
#
include_directories(${Geant4_INCLUDE_DIR}
        ${PROJECT_SOURCE_DIR}
        ${PROJECT_SOURCE_DIR}/DataInit/inc
        ${PROJECT_SOURCE_DIR}/Utils/inc
        ${PROJECT_SOURCE_DIR}/Simulation/inc/
        ${PROJECT_SOURCE_DIR}/modelValidationTests/inc/
        )

#----------------------------------------------------------------------------
# Build target
#
add_executable(dpm_GenerateData ${PROJECT_SOURCE_DIR}/dpm_GenerateData.cc ${cxx_sources_GEN})
target_link_libraries(dpm_GenerateData ${Geant4_LIBRARIES})

add_executable(dpm_Simulate ${PROJECT_SOURCE_DIR}/dpm_Simulate.cc ${cxx_sources_SIM})
target_link_libraries(dpm_Simulate)

## add model level test for bremsstrahlung
add_executable(test_brem ${PROJECT_SOURCE_DIR}/test_brem.cc ${cxx_sources_SIM})
target_link_libraries(test_brem)


add_library(dpm_validation STATIC ${validation_sources} ${cxx_sources_SIM})
target_link_libraries(dpm_validation)
target_include_directories(dpm_validation PUBLIC
        ${PROJECT_SOURCE_DIR}/modelValidationTests/inc
        ${PROJECT_SOURCE_DIR}/Simulation/inc)

add_executable(modeltests model_validation.cpp)
target_link_libraries(modeltests dpm_validation)

