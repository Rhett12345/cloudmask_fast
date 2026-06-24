cmake_minimum_required(VERSION 3.18)
project(fylat VERSION 0.4.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Release build with native optimizations
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")

# ---------------------------------------------------------------------------
# pybind11
# ---------------------------------------------------------------------------
find_package(pybind11 REQUIRED)

# ---------------------------------------------------------------------------
# FYLAT core C++ library (header-only for now, will add .cpp sources later)
# ---------------------------------------------------------------------------
add_library(fylat_cpp INTERFACE)
target_include_directories(fylat_cpp INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# ---------------------------------------------------------------------------
# pybind11 module
# ---------------------------------------------------------------------------
pybind11_add_module(fylat_core
    src/bindings/pybind_module.cpp
    src/core/bitops.cpp
)
target_include_directories(fylat_core PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(fylat_core PRIVATE fylat_cpp)

# ---------------------------------------------------------------------------
# Tests (optional)
# ---------------------------------------------------------------------------
option(BUILD_TESTS "Build C++ unit tests" OFF)
if(BUILD_TESTS)
    enable_testing()
    find_package(GTest QUIET)
    if(GTest_FOUND)
        add_executable(test_bitops tests/test_bitops.cpp src/core/bitops.cpp)
        target_include_directories(test_bitops PRIVATE include)
        target_link_libraries(test_bitops GTest::GTest GTest::Main)
        add_test(NAME bitops COMMAND test_bitops)
    endif()
endif()
