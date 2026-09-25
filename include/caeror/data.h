#pragma once
#include <cmath>

struct Point{
  double x, y, z;
};

struct Direction{
  double mu, phi;
  double x, y, z;
  void UpdateXYZ() {
    using std::sqrt, std::sin, std::cos;
    const double sin_theta = sqrt(1.0 - mu * mu);
    x = sin_theta * cos(phi);
    y = sin_theta * sin(phi);
    z = mu;
  }
};