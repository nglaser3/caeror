# --------------------------------------------------------------------------- #
# Setup caeror library 
# --------------------------------------------------------------------------- #

blt_add_library(
    NAME
        caeror
    SOURCES
        caeror_geometry.cc
        geometry_constructor.cc

    HEADERS
        ${CAEROR_INCLUDE_DIR}/caeror/geometry_constructor.h
    DEPENDS_ON
        RAJA
)

target_include_directories(caeror
    PUBLIC ${CAEROR_INCLUDE_DIR}
)