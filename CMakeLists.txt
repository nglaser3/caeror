# --------------------------------------------------------------------------- #
# Root level CMake file for caeror
# --------------------------------------------------------------------------- #

cmake_minimum_required(VERSION 3.23)

project(
  caeror
  VERSION 0.0.0
  LANGUAGES CXX C
)
set(BLT_CXX_STD "c++23" CACHE STRING "")

include(cmake/CaerorPathVariables.cmake)
include(cmake/CaerorVariables.cmake)
include(cmake/CaerorSetupBLT.cmake)
include(cmake/CaerorSetupRAJA.cmake)

add_subdirectory("${CAEROR_SRC_DIR}")

if(CAEROR_BUILD_TESTS)
  enable_testing()
  add_subdirectory("${CAEROR_TESTS_DIR}")
endif()