#pragma once

#include "caeror/geometry_constructor.h"

namespace caeror {

// Two infinite slabs: 0 < x < 2 and 2 < x < 5 centimeters.
inline GeometryConstructor CreateSlabs() {
  GeometryConstructor input;
  auto left = input.CreateAxisAlignedPlane(Axis::X, 0.0);
  auto middle = input.CreateAxisAlignedPlane(Axis::X, 2.0);
  auto right = input.CreateAxisAlignedPlane(Axis::X, 5.0);
  auto first = input.CreateCell(+left && -middle, MaterialID{1});
  auto second = input.CreateCell(+middle && -right, MaterialID{2});
  input.CreateUniverse({first, second});
  return input;
}

}  // namespace caeror
