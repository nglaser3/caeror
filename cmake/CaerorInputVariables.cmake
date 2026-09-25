# --------------------------------------------------------------------------- #
# Setup compile variables for caeror
# --------------------------------------------------------------------------- #

option(CAEROR_BUILD_TESTS "" ON)
if(CAEROR_BUILD_TESTS)
  set(ENABLE_GTEST ON
    CACHE BOOL ""
  )
endif()