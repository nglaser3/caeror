#pragma once
#include <cmath>

#include <RAJA/RAJA.hpp>
#include "caeror/raja_layouts.h"

namespace caeror
{
struct Point{
  double x, y, z;
};

struct Direction{
  double mu, phi;
  double x, y, z;

  RAJA_HOST_DEVICE
  void UpdateXYZ() {
    const double sin_theta = sqrt(1.0 - mu * mu);
    x = sin_theta * cos(phi);
    y = sin_theta * sin(phi);
    z = mu;
  }
};

struct Particle{
  double x, y, z;
  double mu, phi;
  double ux, uy, uz;
  
  CellID cell;
};
} // namespace caeror