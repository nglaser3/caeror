#pragma once
#include "caeror/raja_layouts.h"

namespace caeror {
/**
 * @brief Simple container for minimum data required to represent a particle in geometry
 */
struct Particle {
  /// @name Spatial position
  /// @{
  double x, y, z;
  /// @}

  /// @name Unit direction
  /// @{
  double ux, uy, uz;
  /// @}

  /// @brief Current Cell
  CellID cell;
};
}  // namespace caeror
