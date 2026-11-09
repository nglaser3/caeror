# --------------------------------------------------------------------------- #
# Setup RAJA for caeror
# --------------------------------------------------------------------------- #

if(NOT TARGET RAJA)
  find_package(RAJA CONFIG QUIET)

  if(NOT RAJA_FOUND)
    if(NOT EXISTS "${CAEROR_TPL_DIR}/raja/CMakeLists.txt")
      message(FATAL_ERROR
        "RAJA was not found and the git submodule has not been checked out."
      )
    endif()

    set(RAJA_ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(RAJA_ENABLE_EXERCISES OFF CACHE BOOL "" FORCE)
    set(RAJA_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
    set(CAMP_ENABLE_TESTS OFF CACHE BOOL "" FORCE)

    add_subdirectory(${CAEROR_TPL_DIR}/raja)

    # raja uses some deprecated std functions, turn off the compile error locally
    target_compile_options(RAJA INTERFACE -Wno-deprecated-declarations)
    target_compile_options(camp INTERFACE -Wno-deprecated-declarations)

  endif()
endif()