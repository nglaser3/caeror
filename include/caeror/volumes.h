#pragma once

#include <vector>

#include "caeror/raja_layouts.h"
#include "caeror/region.h"

namespace caeror
{
  struct Cell {
    const CellID id_;
    Region region_;
    MaterialID material_fill_;
    UniverseID universe_fill_;
  };

  struct Universe {
    const UniverseID id_;
    std::vector<CellID> cells_;
  };
} // namespace caeror
