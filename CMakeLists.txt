# --------------------------------------------------------------------------- #
# Root level CMake file for caeror
# --------------------------------------------------------------------------- #

cmake_minimum_required(VERSION 3.24)

project(
  caeror
  VERSION 0.0.0
  LANGUAGES CXX C
)
set(BLT_CXX_STD "c++23" CACHE STRING "C++ standard used by BLT")

include(cmake/CaerorPathVariables.cmake)
include(cmake/CaerorVariables.cmake)
include(cmake/CaerorSetupBLT.cmake)
include(cmake/CaerorSetupRAJA.cmake)

add_subdirectory("${CAEROR_SRC_DIR}")
target_compile_features(caeror PUBLIC cxx_std_23)

if(CAEROR_BUILD_TESTS)
  enable_testing()
  add_subdirectory("${CAEROR_TESTS_DIR}")
endif()