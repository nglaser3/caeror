# --------------------------------------------------------------------------- #
# Configure local API documentation
# --------------------------------------------------------------------------- #

configure_file(doxygen/Doxyfile.in Doxyfile @ONLY)

add_custom_target(caeror_docs
  COMMAND doxygen "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile"
  WORKING_DIRECTORY "${CAEROR_ROOT_DIR}"
  COMMENT "Generating Caeror API documentation"
  VERBATIM
)
add_dependencies(docs caeror_docs)
