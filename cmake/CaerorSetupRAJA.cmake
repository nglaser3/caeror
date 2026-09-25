# --------------------------------------------------------------------------- #
# Setup RAJA for caeror
# --------------------------------------------------------------------------- #

find_package(RAJA QUIET)
if(NOT RAJA_FOUND)
  set(RAJA_ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
  set(RAJA_ENABLE_EXERCISES OFF CACHE BOOL "" FORCE)
  set(RAJA_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
  set(CAMP_ENABLE_TESTS OFF CACHE BOOL "" FORCE)

  add_subdirectory(${CAEROR_TPL_DIR}/raja)
  # raja uses some deprecated std functions, turn off the compile error locally
  target_compile_options(RAJA INTERFACE -Wno-deprecated-declarations)
  target_compile_options(camp INTERFACE -Wno-deprecated-declarations)
endif(NOT RAJA_FOUND)
list(APPEND RCANA_EXTERNAL_DEPS RAJA)