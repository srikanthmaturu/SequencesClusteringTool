cmake_minimum_required(VERSION 3.4)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")
include(AppendCompilerFlags)

enable_testing()

project(sequences_clustering_tool C CXX)

append_cxx_compiler_flags("-msse4.2 -std=c++14 -Wall -DNDEBUG -fopenmp" "GCC" CMAKE_CXX_FLAGS)
append_cxx_compiler_flags("-O3 -ffast-math -funroll-loops" "GCC" CMAKE_CXX_FLAGS)
set(CMAKE_CXX_STANDARD 14)

INCLUDE_DIRECTORIES(${CMAKE_HOME_DIRECTORY}/include
        ${CMAKE_HOME_DIRECTORY}/src
        ${CMAKE_HOME_DIRECTORY}/include/xxhash/
        ${CMAKE_HOME_DIRECTORY}/external/FALCONN/src/include/
        ${CMAKE_HOME_DIRECTORY}/external/FALCONN/src/tf_idf_index/
        ${CMAKE_HOME_DIRECTORY}/external/FALCONN/external/eigen
        ${CMAKE_HOME_DIRECTORY}/external/edlib/edlib/include/)

add_library(xxhash include/xxhash/xxhash.c)
add_library(edlib external/edlib/edlib/src/edlib.cpp)

ADD_SUBDIRECTORY(lib)

ADD_EXECUTABLE(sequences_cluster_generator src/main.cpp)

LINK_LIBRARIES(sequences_cluster_generator edlib xxhash)
