#pragma once

#include <vector>

#include "caeror/raja_layouts.h"
#include "caeror/region.h"

namespace caeror {
/**
 * @brief Container holding information of a cell
 */
struct Cell {
  /// @brief The ID of the cell
  const CellID id_;
  /// @brief Region of the cell
  Region region_;
  /// @brief Material fill (if applicable)
  MaterialID material_fill_{MAXCaerorIndex};
  /// @brief Universe fill (if applicable)
  UniverseID universe_fill_{MAXCaerorIndex};
};

/**
 * @brief Container holding the information of a universe
 */
struct Universe {
  /// @brief The ID of the universe
  const UniverseID id_;
  /// @brief List of cells contained in this universe
  std::vector<CellID> cells_;
};
}  // namespace caeror
