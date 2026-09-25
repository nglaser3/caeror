# --------------------------------------------------------------------------- #
# Setup caeror library 
# --------------------------------------------------------------------------- #

blt_add_library(
    NAME
        caeror
    SOURCES
        surface.cc
    HEADERS
        ${CAEROR_INCLUDE_DIR}/caeror/surface.h
    DEPENDS_ON
        RAJA
)

target_include_directories(caeror
    PUBLIC ${CAEROR_INCLUDE_DIR}
)