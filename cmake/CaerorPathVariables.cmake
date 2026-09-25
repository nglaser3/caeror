# --------------------------------------------------------------------------- #
# Setup global path variables for caeror
# --------------------------------------------------------------------------- #

set(CAEROR_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

cmake_path(APPEND CAEROR_ROOT_DIR "src" OUTPUT_VARIABLE CAEROR_SRC_DIR)
cmake_path(APPEND CAEROR_ROOT_DIR "include" OUTPUT_VARIABLE CAEROR_INCLUDE_DIR)
cmake_path(APPEND CAEROR_ROOT_DIR "tests" OUTPUT_VARIABLE CAEROR_TESTS_DIR)
