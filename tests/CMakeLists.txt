# --------------------------------------------------------------------------- #
# Setup tests for caeror
# --------------------------------------------------------------------------- #

blt_add_executable(
    NAME caeror_tests
    SOURCES
        surfaces/plane.test.cc
        surfaces/axis_aligned_plane.test.cc
        surfaces/sphere.test.cc
        surfaces/axis_aligned_cylinder.test.cc
        surfaces/dispatch.test.cc
    DEPENDS_ON
        caeror
        gtest
)

blt_add_test(
    NAME caeror_tests
    COMMAND caeror_tests
)
