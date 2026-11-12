#pragma once
#include "caeror/raja_layouts.h"

namespace caeror {
struct Particle {
  double x, y, z;
  double ux, uy, uz;

  CellID cell;
};
}  // namespace caeror
