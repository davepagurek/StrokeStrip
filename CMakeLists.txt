cmake_minimum_required(VERSION 3.3)

project(strokestrip C)

# Set c++11
# https://stackoverflow.com/questions/10851247/how-to-activate-c-11-in-cmake
if (POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW)
endif ()
set (CMAKE_CXX_STANDARD 14)

option(CXX "enable C++ compilation" ON)
if(CXX)
    enable_language(CXX)
endif()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Visual Studio compiler with static runtime libraries
if(MSVC AND MT)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
endif()

find_package(GUROBI REQUIRED)

include_directories(${GUROBI_INCLUDE_DIRS})

# Source files
file(GLOB SOURCE_FILES src/*.cpp src/*.h)
add_executable(${CMAKE_PROJECT_NAME} ${SOURCE_FILES})

if(CXX)
    set(CMAKE_CXX_STANDARD 11)
    target_link_libraries(${CMAKE_PROJECT_NAME} optimized ${GUROBI_CXX_LIBRARY} debug ${GUROBI_CXX_DEBUG_LIBRARY})
endif()

target_link_libraries(${CMAKE_PROJECT_NAME} ${GUROBI_LIBRARY})

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    include(FeatureSummary)
    feature_summary(WHAT ALL)
endif()

# Eigen
find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(${CMAKE_PROJECT_NAME} Eigen3::Eigen)

# Glm
add_subdirectory(external/glm)
target_link_libraries(${CMAKE_PROJECT_NAME} glm)

# Cornucopia
add_subdirectory(external/Cornucopia)
target_link_libraries(${CMAKE_PROJECT_NAME} Cornucopia)

target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE external/)
