cmake_minimum_required(VERSION 3.16)
project(radcor VERSION 0.4.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

# Parity with the Fortran needs IEEE-exact float arithmetic: no fast-math, no FMA contraction.
add_compile_options(-Wall -Wextra -ffp-contract=off)

option(RADCOR_WITH_ROOT "Build ROOT TTree output if ROOT is found" ON)

add_library(radcor src/Input.cpp src/Run.cpp src/Options.cpp src/RootIO.cpp)
target_include_directories(radcor PUBLIC include)

if(RADCOR_WITH_ROOT)
  find_package(ROOT QUIET COMPONENTS Tree RIO)
  if(ROOT_FOUND)
    message(STATUS "radcor: ROOT ${ROOT_VERSION} found, ROOT output enabled")
    target_compile_definitions(radcor PRIVATE RADCOR_HAS_ROOT)
    target_link_libraries(radcor PUBLIC ROOT::Tree ROOT::RIO)
  else()
    message(STATUS "radcor: ROOT not found, ROOT output disabled")
  endif()
endif()

add_executable(radcor_app apps/radcor.cpp)
set_target_properties(radcor_app PROPERTIES OUTPUT_NAME radcor)
target_link_libraries(radcor_app PRIVATE radcor)

enable_testing()
add_test(NAME parity_fortran COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_parity.sh)
add_test(NAME input_formats COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_formats.sh)
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_FOUND)
  add_test(NAME python_interface COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_python.py)
endif()
find_program(GFORTRAN gfortran)
if(GFORTRAN)
  add_test(NAME fortran_rebuild COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fortran_parity.sh)
endif()

# Example: the library used from C++ (examples/engine_from_csv.cpp)
add_executable(engine_from_csv examples/engine_from_csv.cpp)
target_link_libraries(engine_from_csv PRIVATE radcor)
