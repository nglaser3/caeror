# --------------------------------------------------------------------------- #
# Configure coverage for caeror
# --------------------------------------------------------------------------- #

option(CAEROR_GENERATE_COVERAGE "Enable coverage" OFF)

if(CAEROR_GENERATE_COVERAGE)
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR "CAEROR_GENERATE_COVERAGE requires GCC")
  endif()

  target_compile_options(caeror PUBLIC --coverage -O0 -g)
  target_link_options(caeror PUBLIC --coverage)
endif()
