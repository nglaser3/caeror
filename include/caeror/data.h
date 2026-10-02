#pragma once
#include <cmath>

#include "caeror/raja_layouts.h"
#include <RAJA/RAJA.hpp>

namespace caeror {
struct Particle {
  double x, y, z;
  double mu, phi;
  double ux, uy, uz;

  CellID cell;

  RAJA_HOST_DEVICE
  void UpdateXYZ() {
    const double sin_theta = sqrt(1.0 - mu * mu);
    x = sin_theta * cos(phi);
    y = sin_theta * sin(phi);
    z = mu;
  }
};
} // namespace caeror