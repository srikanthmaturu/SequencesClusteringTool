cmake_minimum_required(VERSION 3.4)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")
include(AppendCompilerFlags)

enable_testing()

project(sequences_clustering_tool C CXX)

append_cxx_compiler_flags("-msse4.2 -std=c++14 -Wall -DNDEBUG -fopenmp" "GCC" CMAKE_CXX_FLAGS)
append_cxx_compiler_flags("-O3 -ffast-math -funroll-loops" "GCC" CMAKE_CXX_FLAGS)
set(CMAKE_CXX_STANDARD 14)


set(CMAKE_PREFIX_PATH  ${CMAKE_PREFIX_PATH} ${CMAKE_HOME_DIRECTORY}/external/seqan/util/cmake/)
set(SEQAN_INCLUDE_PATH "${CMAKE_HOME_DIRECTORY}/external/seqan/include")

# Load the SeqAn module and fail if not found.
find_package (SeqAn REQUIRED)
MESSAGE("SEQAN status ${SEQAN_FOUND}")

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

INCLUDE_DIRECTORIES(${CMAKE_HOME_DIRECTORY}/include
        ${CMAKE_HOME_DIRECTORY}/src
        ${CMAKE_HOME_DIRECTORY}/include/xxhash/
        ${CMAKE_HOME_DIRECTORY}/external/FALCONN/src/include/
        ${CMAKE_HOME_DIRECTORY}/external/FALCONN/src/tf_idf_index/
        ${CMAKE_HOME_DIRECTORY}/external/FALCONN/external/eigen
        ${CMAKE_HOME_DIRECTORY}/external/edlib/edlib/include/
        ${SEQAN_INCLUDE_DIRS})

add_library(xxhash include/xxhash/xxhash.c)
add_library(edlib external/edlib/edlib/src/edlib.cpp)

ADD_SUBDIRECTORY(lib)

FILE(STRINGS ${CMAKE_HOME_DIRECTORY}/exp0.config exp0.lines REGEX "^[^#].*")
FOREACH(line ${exp0.lines})
    MESSAGE("line = ${line}")
    LIST(GET line 0 index_name)
    LIST(GET line 1 index_type)
    MESSAGE("Exp 0 ${index_name}; ${index_type};")
    STRING(REGEX REPLACE "<(.*)>" "\\1" t_part ${index_type})
    STRING(REGEX REPLACE "," ";" t_params ${t_part})
    LIST(GET t_params 0 NGRAM_LENGTH)
    LIST(GET t_params 1 LSH_HASH_TYPE)
    LIST(GET t_params 2 NUMBER_OF_HASH_TABLES)
    LIST(GET t_params 3 NUMBER_OF_HASH_BITS)
    LIST(GET t_params 4 NUMBER_OF_PROBES)
    LIST(GET t_params 5 THRESHOLD)
    LIST(GET t_params 6 POINT_TYPE)
    SET(exec ${index_name}_NL_${NGRAM_LENGTH}_LT_${LSH_HASH_TYPE}_NHT_${NUMBER_OF_HASH_TABLES}_NHB_${NUMBER_OF_HASH_BITS}_NP_${NUMBER_OF_PROBES}_TH_${THRESHOLD}_PT_${POINT_TYPE})
    ADD_EXECUTABLE(${exec} src/main.cpp)
    TARGET_LINK_LIBRARIES(${exec} xxhash edlib ${SEQAN_LIBRARIES})
    STRING(COMPARE EQUAL ${POINT_TYPE} DenseVectorFloat IsDenseVectorFloatType)
    SET_PROPERTY(TARGET ${exec} PROPERTY COMPILE_DEFINITIONS
            INDEX_NAME="${index_name}"
            INDEX_TYPE=${index_type}
            NGRAM_LENGTH=${NGRAM_LENGTH}
            LSH_HASH_TYPE=${LSH_HASH_TYPE}
            NUMBER_OF_HASH_TABLES=${NUMBER_OF_HASH_TABLES}
            NUMBER_OF_HASH_BITS=${NUMBER_OF_HASH_BITS}
            NUMBER_OF_PROBES=${NUMBER_OF_PROBES}
            THRESHOLD=${THRESHOLD}
            POINT_TYPE=${POINT_TYPE}
            )
    if(IsDenseVectorFloatType)
        GET_TARGET_PROPERTY(current_properties ${exec} COMPILE_DEFINITIONS)
        SET_PROPERTY(TARGET ${exec} PROPERTY COMPILE_DEFINITIONS ${current_properties}
                VT_DVF=true)
        message("Set VT_DVF compile definition")
    endif()
    message("Point type ${POINT_TYPE}")
    message("Target ${exec}")
ENDFOREACH()
