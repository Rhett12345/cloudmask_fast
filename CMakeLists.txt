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
find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

message(STATUS "HDF5 include dirs: ${HDF5_INCLUDE_DIRS}")
message(STATUS "Python executable: ${Python_EXECUTABLE}")

# Eigen for linear algebra (NWP interpolation, Phase 2+)
find_package(Eigen3 QUIET)
if(Eigen3_FOUND)
    message(STATUS "Eigen3 found: ${EIGEN3_INCLUDE_DIRS}")
endif()

# ============================================================================
# Phase 1: C++ HDF5 IO library and Python bindings
# ============================================================================
add_library(fylat_io STATIC
    src/cpp/io/hdf5_io.cpp
)
target_include_directories(fylat_io PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp
    ${HDF5_INCLUDE_DIRS}
)
target_link_libraries(fylat_io PUBLIC
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_LIBRARIES}
)

pybind11_add_module(fylat_py
    src/python_bindings/bindings.cpp
)
target_link_libraries(fylat_py PRIVATE fylat_io)

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
add_test(
    NAME cpp_io_smoke
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cpp_io_smoke.py
)
set_tests_properties(cpp_io_smoke PROPERTIES
    ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:fylat_py>"
)
add_test(
    NAME nwp_reader_layout
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_nwp_reader_layout.py
)
add_test(
    NAME mersi_io_cpp_backend
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_mersi_io_cpp_backend.py
)
set_tests_properties(mersi_io_cpp_backend PROPERTIES
    ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:fylat_py>:${CMAKE_CURRENT_SOURCE_DIR}/python;FYLAT_IO_BACKEND=cpp"
)

add_custom_target(check
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

message(STATUS "FYLAT v3.2 — Phase 1 IO migration build configured")
