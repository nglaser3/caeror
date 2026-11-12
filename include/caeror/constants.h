#pragma once

#include <limits>

namespace caeror {

// Absolute geometry tolerance in centimeters.
inline constexpr double TOLERANCE = 1e-10;
inline constexpr double INF = std::numeric_limits<double>::infinity();

}  // namespace caeror
