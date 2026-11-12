#pragma once

#include "caeror/geometry_constructor.h"

namespace caeror {

// Kobayashi dog-leg geometry.
// Dimensions are centimeters; bounds are (0,0,0) to (60,100,60).
// Material and cell IDs: shield = 0, void = 1, source = 2.
inline GeometryConstructor CreateKobayashi() {
  GeometryConstructor geometry;

  // x=0, y=0, z=0 are reflective; x=60, y=100, z=60 are vacuum.
  const auto x1 = geometry.CreateAxisAlignedPlane(
      Axis::X, 0.0, BoundaryCondition::Reflective);
  const auto x2 = geometry.CreateAxisAlignedPlane(Axis::X, 10.0);
  const auto x3 = geometry.CreateAxisAlignedPlane(Axis::X, 30.0);
  const auto x4 = geometry.CreateAxisAlignedPlane(Axis::X, 40.0);
  const auto x5 =
      geometry.CreateAxisAlignedPlane(Axis::X, 60.0, BoundaryCondition::Vacuum);
  const auto y1 = geometry.CreateAxisAlignedPlane(
      Axis::Y, 0.0, BoundaryCondition::Reflective);
  const auto y2 = geometry.CreateAxisAlignedPlane(Axis::Y, 10.0);
  const auto y3 = geometry.CreateAxisAlignedPlane(Axis::Y, 50.0);
  const auto y4 = geometry.CreateAxisAlignedPlane(Axis::Y, 60.0);
  const auto y5 = geometry.CreateAxisAlignedPlane(Axis::Y, 100.0,
                                                  BoundaryCondition::Vacuum);
  const auto z1 = geometry.CreateAxisAlignedPlane(
      Axis::Z, 0.0, BoundaryCondition::Reflective);
  const auto z2 = geometry.CreateAxisAlignedPlane(Axis::Z, 10.0);
  const auto z3 = geometry.CreateAxisAlignedPlane(Axis::Z, 30.0);
  const auto z4 = geometry.CreateAxisAlignedPlane(Axis::Z, 40.0);
  const auto z5 =
      geometry.CreateAxisAlignedPlane(Axis::Z, 60.0, BoundaryCondition::Vacuum);

  const auto source_region = +x1 && -x2 && +y1 && -y2 && +z1 && -z2;
  const auto channel_1 = +x1 && -x2 && +y2 && -y3 && +z1 && -z2;
  const auto channel_2 = +x1 && -x3 && +y3 && -y4 && +z1 && -z2;
  const auto channel_3 = +x3 && -x4 && +y3 && -y4 && +z1 && -z3;
  const auto channel_4 = +x3 && -x4 && +y3 && -y5 && +z3 && -z4;
  const auto void_region = channel_1 || channel_2 || channel_3 || channel_4;
  const auto problem_box = +x1 && -x5 && +y1 && -y5 && +z1 && -z5;
  const auto shield_region = problem_box && !void_region && !source_region;

  const auto shield_cell = geometry.CreateCell(shield_region, MaterialID{0});
  const auto void_cell = geometry.CreateCell(void_region, MaterialID{1});
  const auto source_cell = geometry.CreateCell(source_region, MaterialID{2});
  geometry.CreateUniverse({shield_cell, void_cell, source_cell});

  return geometry;
}
}  // namespace caeror
