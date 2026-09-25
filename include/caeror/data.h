#pragma once
#include <cmath>

#include <RAJA/RAJA.hpp>

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
} // namespace caeror