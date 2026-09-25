# --------------------------------------------------------------------------- #
# Setup BLT
# --------------------------------------------------------------------------- #
if(NOT DEFINED BLT_SOURCE_DIR)
  set(BLT_SOURCE_DIR
    "${CAEROR_ROOT_DIR}/tpl/blt"
    CACHE PATH ""
  )
endif()

if(NOT EXISTS "${BLT_SOURCE_DIR}/SetupBLT.cmake")
  message(FATAL_ERROR "BLT not found: ${BLT_SOURCE_DIR}")
endif()

include("${BLT_SOURCE_DIR}/SetupBLT.cmake")