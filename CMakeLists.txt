cmake_minimum_required(VERSION 3.20)
project(fylat LANGUAGES C CXX Fortran
        DESCRIPTION "FYLAT — FY-3D MERSI-II Cloud Mask Retrieval System"
        VERSION 3.2.0)

# ============================================================================
# Phase 1: IO Migration — C++ IO libraries
# ============================================================================
# Currently Fortran core is compiled via legacy Makefile (build.sh).
# This CMake build covers the new C++ IO components and will gradually
# replace Fortran modules across phases 1-4.

# --- C++ standard ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# --- Dependencies ---
# HDF5 C++ wrapper (header-only, optional for Phase 1)
# find_package(HighFive QUIET)
# if(HighFive_FOUND)
#     message(STATUS "HighFive found: ${HighFive_INCLUDE_DIRS}")
# else()
#     message(STATUS "HighFive not found — HDF5 IO via Python h5py (Phase 1)")
# endif()

# pybind11 for Python bindings
find_package(pybind11 QUIET)
if(pybind11_FOUND)
    message(STATUS "pybind11 found: ${pybind11_INCLUDE_DIRS}")
else()
    message(STATUS "pybind11 not found — Python bindings deferred")
endif()

# Eigen for linear algebra (NWP interpolation, Phase 2+)
find_package(Eigen3 QUIET)
if(Eigen3_FOUND)
    message(STATUS "Eigen3 found: ${EIGEN3_INCLUDE_DIRS}")
endif()

# ============================================================================
# Phase 1: C++ MERSI L1 HDF5 Reader (HighFive-based, deferred until available)
# ============================================================================
# add_library(fylat_io STATIC
#     src/cpp/io/mersi_l1_reader.cpp
#     src/cpp/io/nwp_grib_reader.cpp
#     src/cpp/io/clm_writer.cpp
# )
# target_include_directories(fylat_io PUBLIC
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp
#     ${HighFive_INCLUDE_DIRS}
# )
# target_link_libraries(fylat_io PUBLIC HighFive::HighFive)

# ============================================================================
# Phase 1: pybind11 Python module (deferred until HighFive available)
# ============================================================================
# if(pybind11_FOUND)
#     pybind11_add_module(fylat_py
#         src/python_bindings/bindings.cpp
#     )
#     target_link_libraries(fylat_py PRIVATE fylat_io)
# endif()

# ============================================================================
# Legacy Fortran build (via external Makefile for now)
# ============================================================================
add_custom_target(fortran_build
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/build.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Building Fortran core via legacy Makefile"
)

# ============================================================================
# Tests
# ============================================================================
enable_testing()
add_custom_target(check
    COMMAND ${CMAKE_COMMAND} -E echo "Phase 1 verification: run 3 orbits"
    COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/run_fylat.py --date 20220803
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

message(STATUS "FYLAT v3.2 — Phase 1 IO migration build configured")
