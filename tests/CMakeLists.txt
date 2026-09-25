# --------------------------------------------------------------------------- #
# Setup tests for caeror
# --------------------------------------------------------------------------- #

blt_add_executable(
    NAME caeror_tests
    SOURCES
        surface.test.cc
    DEPENDS_ON
        caeror
        gtest
)

blt_add_test(
    NAME caeror_tests
    COMMAND caeror_tests
)