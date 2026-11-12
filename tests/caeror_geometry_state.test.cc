#include <gtest/gtest.h>

#include "caeror/caeror_geometry.h"
#include "caeror/constants.h"
#include "sample_constructors/kobayashi.h"
#include "sample_constructors/slabs.h"

namespace {
using namespace caeror;

TEST(CaerorGeometryState, Intersection) {
  GeometryConstructor input;
  auto x = input.CreateAxisAlignedPlane(Axis::X, 0.0);
  auto y = input.CreateAxisAlignedPlane(Axis::Y, 0.0);
  input.CreateCell(+x && +y, MaterialID{0});

  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();

  Particle p{};
  p.ux = 1.0;

  p.x = 1.0;
  p.y = 1.0;
  EXPECT_TRUE(state.InCell(p, CellID{0}));

  p.x = -1.0;
  EXPECT_FALSE(state.InCell(p, CellID{0}));

  p.x = 1.0;
  p.y = -1.0;
  EXPECT_FALSE(state.InCell(p, CellID{0}));
}

TEST(CaerorGeometryState, Union) {
  GeometryConstructor input;
  auto x = input.CreateAxisAlignedPlane(Axis::X, 0.0);
  auto y = input.CreateAxisAlignedPlane(Axis::Y, 0.0);
  input.CreateCell(+x || +y, MaterialID{0});
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.ux = 1.0;
  p.x = 1.0;
  p.y = -1.0;
  auto result = state.InCell(p, CellID{0});
  EXPECT_TRUE(result);

  p.x = -1.0;
  p.y = 1.0;
  result = state.InCell(p, CellID{0});
  EXPECT_TRUE(result);

  p.y = -1.0;
  result = state.InCell(p, CellID{0});
  EXPECT_FALSE(result);
}

TEST(CaerorGeometryState, Complement) {
  GeometryConstructor input;
  auto x = input.CreateAxisAlignedPlane(Axis::X, 0.0);
  input.CreateCell(!(+x), MaterialID{0});
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.ux = 1.0;
  p.x = -1.0;
  auto result = state.InCell(p, CellID{0});
  EXPECT_TRUE(result);

  p.x = 1.0;
  result = state.InCell(p, CellID{0});
  EXPECT_FALSE(result);
}

TEST(CaerorGeometryState, NestedBoolean) {
  GeometryConstructor input;
  auto x = input.CreateAxisAlignedPlane(Axis::X, 0.0);
  auto y = input.CreateAxisAlignedPlane(Axis::Y, 0.0);
  auto z = input.CreateAxisAlignedPlane(Axis::Z, 0.0);
  input.CreateCell(+x || (+y && +z), MaterialID{0});
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.ux = 1.0;
  p.x = 1.0;
  p.y = -1.0;
  p.z = -1.0;
  // Inside +x, so the union must be true even though both other sides are
  // false.
  auto result = state.InCell(p, CellID{0});
  EXPECT_TRUE(result);
}

TEST(CaerorGeometryState, FindSlabOrOutside) {
  auto input = CreateSlabs();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  const double positions[] = {-1.0, 1.0, 3.0, 6.0};
  const CellID expected[] = {CellID{MAXCaerorIndex}, CellID{0}, CellID{1},
                             CellID{MAXCaerorIndex}};
  for (int i = 0; i < 4; ++i) {
    Particle p{};
    p.x = positions[i];
    p.ux = 1.0;
    auto result = state.FindCell(p);
    EXPECT_EQ(result, expected[i]);
  }
}

TEST(CaerorGeometryState, EmptyGeometryHasNoCell) {
  GeometryConstructor input;
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.ux = 1.0;
  auto result = state.FindCell(p);
  EXPECT_EQ(result, CellID{MAXCaerorIndex});
}

TEST(CaerorGeometryState, BoundaryMembershipUsesDirection) {
  auto input = CreateSlabs();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  for (double direction : {-1.0, 1.0}) {
    Particle p{};
    p.x = 2.0;
    p.ux = direction;
    p.cell = CellID{MAXCaerorIndex};
    auto first = state.InCell(p, CellID{0});
    EXPECT_EQ(first, direction < 0.0);
    auto second = state.InCell(p, CellID{1});
    EXPECT_EQ(second, direction > 0.0);
    auto result = state.FindCell(p);
    EXPECT_EQ(result, direction < 0.0 ? CellID{0} : CellID{1});
  }
}

TEST(CaerorGeometryState, TangentPlaneRetainsCurrentCell) {
  auto input = CreateSlabs();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  for (CaerorIndexType i = 0; i < 2; ++i) {
    Particle p{};
    p.x = 2.0;
    p.uy = 1.0;
    p.cell = CellID{i};
    auto found = state.FindCell(p);
    EXPECT_EQ(found, p.cell);
    auto result = state.CrossSurface(p, SurfaceID{1});
    EXPECT_EQ(result.cell, p.cell);
    EXPECT_FALSE(result.leaked);
  }
}

TEST(CaerorGeometryState, TangencyDoesNotOverride) {
  GeometryConstructor input;
  auto x = input.CreateAxisAlignedPlane(Axis::X, 0.0);
  auto y = input.CreateAxisAlignedPlane(Axis::Y, 0.0);
  input.CreateCell(+x && +y, MaterialID{0});

  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();

  Particle p{};
  p.x = 0.0;
  p.y = -1.0;
  p.uy = 1.0;
  p.cell = CellID{0};

  auto result = state.InCell(p, CellID{0});
  EXPECT_FALSE(result);
}

TEST(CaerorGeometryState, SlabDistances) {
  auto input = CreateSlabs();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  for (double direction : {-1.0, 1.0}) {
    Particle p{};
    p.x = 0.5;
    p.ux = direction;
    p.cell = CellID{0};
    auto result = state.DistanceToSurface(p);
    EXPECT_DOUBLE_EQ(result.distance, direction < 0.0 ? 0.5 : 1.5);
    EXPECT_EQ(result.surface, direction < 0.0 ? SurfaceID{0} : SurfaceID{1});
  }
}

TEST(CaerorGeometryState, ParallelRayHasNoBoundary) {
  auto input = CreateSlabs();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.x = 1.0;
  p.uy = 1.0;
  p.cell = CellID{0};
  auto result = state.DistanceToCellBoundary(p);
  EXPECT_EQ(result.distance, INF);
  EXPECT_EQ(result.surface, SurfaceID{MAXCaerorIndex});
}

TEST(CaerorGeometryState, StartingOnEntranceSkipsZeroDistance) {
  auto input = CreateSlabs();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.x = 2.0;
  p.ux = 1.0;
  p.cell = CellID{1};
  auto result = state.DistanceToSurface(p);
  EXPECT_DOUBLE_EQ(result.distance, 3.0);
  EXPECT_EQ(result.surface, SurfaceID{2});
  auto boundary = state.DistanceToCellBoundary(p);
  EXPECT_DOUBLE_EQ(boundary.distance, 3.0);
  EXPECT_EQ(boundary.surface, SurfaceID{2});
}

TEST(CaerorGeometryState, StartingOnExitReportsZeroDistance) {
  auto input = CreateSlabs();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.x = 2.0;
  p.ux = 1.0;
  p.cell = CellID{0};
  auto result = state.DistanceToSurface(p);
  EXPECT_DOUBLE_EQ(result.distance, 0.0);
  EXPECT_EQ(result.surface, SurfaceID{1});
  auto boundary = state.DistanceToCellBoundary(p);
  EXPECT_DOUBLE_EQ(boundary.distance, 0.0);
  EXPECT_EQ(boundary.surface, SurfaceID{1});
}

TEST(CaerorGeometryState, UnionSkipsInternalSurface) {
  GeometryConstructor input;
  auto left = input.CreateAxisAlignedPlane(Axis::X, 0.0);
  auto middle = input.CreateAxisAlignedPlane(Axis::X, 2.0);
  auto right = input.CreateAxisAlignedPlane(Axis::X, 5.0);
  input.CreateCell((+left && -middle) || (+middle && -right), MaterialID{0});
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.x = 1.0;
  p.ux = 1.0;
  p.cell = CellID{0};
  auto result = state.DistanceToCellBoundary(p);
  EXPECT_DOUBLE_EQ(result.distance, 4.0);
  EXPECT_EQ(result.surface, right.id_);
  p.x = 2.0;
  auto cross_result = state.CrossSurface(p, middle.id_);
  EXPECT_EQ(cross_result.cell, p.cell);
  EXPECT_FALSE(cross_result.leaked);
}

TEST(CaerorGeometryState, UnionReportsInternalSurface) {
  GeometryConstructor input;
  auto left = input.CreateAxisAlignedPlane(Axis::X, 0.0);
  auto middle = input.CreateAxisAlignedPlane(Axis::X, 2.0);
  auto right = input.CreateAxisAlignedPlane(Axis::X, 5.0);
  input.CreateCell((+left && -middle) || (+middle && -right), MaterialID{0});
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.x = 1.0;
  p.ux = 1.0;
  p.cell = CellID{0};
  auto result = state.DistanceToSurface(p);
  EXPECT_DOUBLE_EQ(result.distance, 1.0);
  EXPECT_EQ(result.surface, middle.id_);
}

TEST(CaerorGeometryState, BoundaryBeyondOverlappingSpheres) {
  GeometryConstructor input;
  auto first = input.CreateSphere(2.0, 0.0, 0.0, 0.0);
  auto second = input.CreateSphere(2.0, 2.0, 0.0, 0.0);
  input.CreateCell(-first || -second, MaterialID{0});
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.x = 0.0;
  p.ux = 1.0;
  p.cell = CellID{0};
  auto result = state.DistanceToCellBoundary(p);
  EXPECT_DOUBLE_EQ(result.distance, 4.0);
  EXPECT_EQ(result.surface, second.id_);
}

TEST(CaerorGeometryState, CurvedSurfaceDistances) {
  GeometryConstructor input;
  auto sphere = input.CreateSphere(2.0, 0.0, 0.0, 0.0);
  auto cylinder = input.CreateAxisAlignedCylinder(Axis::Z, 3.0, 0.0, 0.0);
  input.CreateCell(-sphere, MaterialID{0});
  input.CreateCell(-cylinder, MaterialID{1});
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  // Query each cell directly; these deliberately overlap.
  for (CaerorIndexType i = 0; i < 2; ++i) {
    Particle p{};
    p.x = 0.5;
    p.ux = 1.0;
    p.cell = CellID{i};
    auto inside = state.InCell(p, p.cell);
    EXPECT_TRUE(inside);
    auto result = state.DistanceToSurface(p);
    EXPECT_DOUBLE_EQ(result.distance, i == 0 ? 1.5 : 2.5);
    EXPECT_EQ(result.surface, i == 0 ? sphere.id_ : cylinder.id_);
  }
}

TEST(CaerorGeometryState, ObliquePlaneDistance) {
  GeometryConstructor input;
  auto plane = input.CreatePlane(1.0, 1.0, 0.0, -4.0);
  input.CreateCell(-plane, MaterialID{0});
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.x = 1.0;
  p.y = 1.0;
  p.ux = 0.6;
  p.uy = 0.8;
  p.cell = CellID{0};
  auto result = state.DistanceToSurface(p);
  EXPECT_DOUBLE_EQ(result.distance, 2.0 / 1.4);
  EXPECT_EQ(result.surface, plane.id_);
}

TEST(CaerorGeometryState, TangentSphereIsNotACrossing) {
  GeometryConstructor input;
  auto sphere = input.CreateSphere(1.0, 0.0, 0.0, 0.0);
  auto right = input.CreateAxisAlignedPlane(Axis::X, 3.0);
  input.CreateCell(+sphere && -right, MaterialID{0});
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.x = -2.0;
  p.y = 1.0;
  p.ux = 1.0;
  p.cell = CellID{0};
  auto result = state.DistanceToCellBoundary(p);
  EXPECT_DOUBLE_EQ(result.distance, 5.0);
  EXPECT_EQ(result.surface, right.id_);
}

TEST(CaerorGeometryState, CrossSlabsBothDirections) {
  auto input = CreateSlabs();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  for (double direction : {-1.0, 1.0}) {
    Particle p{};
    p.x = 2.0;
    p.ux = direction;
    p.cell = direction > 0.0 ? CellID{0} : CellID{1};
    auto next = direction > 0.0 ? CellID{1} : CellID{0};
    auto entering = state.InCell(p, next);
    EXPECT_TRUE(entering);
    auto leaving = state.InCell(p, p.cell);
    EXPECT_FALSE(leaving);
    auto result = state.CrossSurface(p, SurfaceID{1});
    EXPECT_EQ(result.cell, next);
    EXPECT_FALSE(result.leaked);
  }
}

TEST(CaerorGeometryState, CrossAwayFromBoundaryKeepsCell) {
  auto input = CreateSlabs();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.x = 1.0;
  p.ux = 1.0;
  p.cell = CellID{0};
  auto result = state.CrossSurface(p, SurfaceID{1});
  EXPECT_EQ(result.cell, p.cell);
  EXPECT_FALSE(result.leaked);
}

TEST(CaerorGeometryState, TransmissionWithoutNeighborLeaks) {
  auto input = CreateSlabs();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.x = 5.0;
  p.ux = 1.0;
  p.cell = CellID{1};
  auto result = state.CrossSurface(p, SurfaceID{2});
  EXPECT_EQ(result.cell, CellID{MAXCaerorIndex});
  EXPECT_TRUE(result.leaked);
}

TEST(CaerorGeometryState, VacuumExitLeaks) {
  GeometryConstructor input;
  auto wall =
      input.CreateAxisAlignedPlane(Axis::X, 2.0, BoundaryCondition::Vacuum);
  input.CreateCell(-wall, MaterialID{0});
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.x = 2.0;
  p.ux = 1.0;
  p.cell = CellID{0};
  auto result = state.CrossSurface(p, wall.id_);
  EXPECT_EQ(result.cell, CellID{MAXCaerorIndex});
  EXPECT_TRUE(result.leaked);
}

TEST(CaerorGeometryState, VacuumRequiresExit) {
  GeometryConstructor input;
  auto wall =
      input.CreateAxisAlignedPlane(Axis::X, 2.0, BoundaryCondition::Vacuum);
  input.CreateCell(-wall, MaterialID{0});
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.x = 1.0;
  p.ux = 1.0;
  p.cell = CellID{0};
  auto result = state.CrossSurface(p, wall.id_);
  EXPECT_EQ(result.cell, p.cell);
  EXPECT_FALSE(result.leaked);

  p.x = 2.0;
  p.ux = -1.0;
  result = state.CrossSurface(p, wall.id_);
  EXPECT_EQ(result.cell, p.cell);
  EXPECT_FALSE(result.leaked);

  p.ux = 0.0;
  p.uy = 1.0;
  result = state.CrossSurface(p, wall.id_);
  EXPECT_EQ(result.cell, p.cell);
  EXPECT_FALSE(result.leaked);
}

TEST(CaerorGeometryState, VacuumWithNeighborTransmits) {
  GeometryConstructor input;
  auto wall =
      input.CreateAxisAlignedPlane(Axis::X, 2.0, BoundaryCondition::Vacuum);
  input.CreateCell(-wall, MaterialID{0});
  input.CreateCell(+wall, MaterialID{1});
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.x = 2.0;
  p.ux = 1.0;
  p.cell = CellID{0};
  // BCs apply only when leaving the geometry, not between neighboring cells.
  auto result = state.CrossSurface(p, wall.id_);
  EXPECT_EQ(result.cell, CellID{1});
  EXPECT_FALSE(result.leaked);
}

TEST(CaerorGeometryState, CrossingChecksEntireCandidateRegion) {
  GeometryConstructor input;
  auto x = input.CreateAxisAlignedPlane(Axis::X, 0.0);
  auto y = input.CreateAxisAlignedPlane(Axis::Y, 0.0);
  input.CreateCell(+x && +y, MaterialID{0});
  input.CreateCell(+x && -y, MaterialID{1});
  input.CreateCell(-x, MaterialID{2});
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.y = -1.0;
  p.ux = 1.0;
  p.cell = CellID{2};
  auto entering = state.InCell(p, CellID{0});
  EXPECT_FALSE(entering);
  auto result = state.CrossSurface(p, x.id_);
  EXPECT_EQ(result.cell, CellID{1});
  EXPECT_FALSE(result.leaked);
}

TEST(CaerorGeometryState, ComplementBoundaryMembership) {
  GeometryConstructor input;
  auto x = input.CreateAxisAlignedPlane(Axis::X, 0.0);
  auto y = input.CreateAxisAlignedPlane(Axis::Y, 0.0);
  input.CreateCell(!(+x && +y), MaterialID{0});
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.y = 1.0;
  p.ux = -1.0;
  auto result = state.InCell(p, CellID{0});
  EXPECT_TRUE(result);
}

TEST(CaerorGeometryState, TraceSlabsWithoutNudging) {
  auto input = CreateSlabs();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.x = 1.0;
  p.ux = 1.0;
  p.cell = state.FindCell(p);
  for (CaerorIndexType i = 0; i < 2; ++i) {
    ASSERT_EQ(p.cell, CellID{i});
    auto material = state.material_fills(p.cell);
    EXPECT_EQ(material, i == 0 ? MaterialID{1} : MaterialID{2});
    auto boundary = state.DistanceToCellBoundary(p);
    ASSERT_DOUBLE_EQ(boundary.distance, i == 0 ? 1.0 : 3.0);
    ASSERT_EQ(boundary.surface, SurfaceID{i + 1});
    p.x += boundary.distance * p.ux;
    auto result = state.CrossSurface(p, boundary.surface);
    EXPECT_EQ(result.leaked, i == 1);
    p.cell = result.cell;
  }
  EXPECT_EQ(p.cell, CellID{MAXCaerorIndex});
  EXPECT_DOUBLE_EQ(p.x, 5.0);
}

TEST(CaerorGeometryState, KobayashiMembership) {
  auto input = CreateKobayashi();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  const double positions[][3] = {{5, 5, 5},    {5, 20, 5},   {20, 55, 5},
                                 {35, 55, 20}, {35, 80, 35}, {50, 80, 50},
                                 {-1, 5, 5},   {61, 5, 5}};
  const CaerorIndexType expected[] = {
      2, 1, 1, 1, 1, 0, MAXCaerorIndex, MAXCaerorIndex};
  for (int i = 0; i < 8; ++i) {
    Particle p{};
    p.x = positions[i][0];
    p.y = positions[i][1];
    p.z = positions[i][2];
    p.uy = 1.0;
    auto result = state.FindCell(p);
    EXPECT_EQ(result, CellID{expected[i]}) << "point=" << i;
  }
}

TEST(CaerorGeometryState, TraceKobayashiAlongY) {
  auto input = CreateKobayashi();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.x = 5.0;
  p.y = 5.0;
  p.z = 5.0;
  p.uy = 1.0;
  p.cell = state.FindCell(p);
  const CaerorIndexType cells[] = {2, 1, 0};
  const CaerorIndexType surfaces[] = {6, 8, 9};
  const double distances[] = {5.0, 50.0, 40.0};
  // y=50 joins two pieces of the same void cell; it is not an external event.
  for (int i = 0; i < 3; ++i) {
    ASSERT_EQ(p.cell, CellID{cells[i]});
    auto material = state.material_fills(p.cell);
    EXPECT_EQ(material, MaterialID{cells[i]});
    auto boundary = state.DistanceToCellBoundary(p);
    ASSERT_DOUBLE_EQ(boundary.distance, distances[i]);
    ASSERT_EQ(boundary.surface, SurfaceID{surfaces[i]});
    p.y += boundary.distance * p.uy;
    auto result = state.CrossSurface(p, boundary.surface);
    EXPECT_EQ(result.leaked, i == 2);
    p.cell = result.cell;
  }
  EXPECT_EQ(p.cell, CellID{MAXCaerorIndex});
  EXPECT_DOUBLE_EQ(p.y, 100.0);
}

TEST(CaerorGeometryState, ReflectObliquePlane) {
  GeometryConstructor input;
  auto wall =
      input.CreatePlane(3.0, 4.0, 0.0, -25.0, BoundaryCondition::Reflective);
  input.CreateCell(-wall, MaterialID{0});

  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();

  Particle p{};
  p.x = 3.0;
  p.y = 4.0;
  p.ux = 1.0;
  p.cell = CellID{0};

  auto result = state.CrossSurface(p, wall.id_);
  EXPECT_EQ(result.cell, CellID{0});
  EXPECT_FALSE(result.leaked);
  EXPECT_DOUBLE_EQ(p.ux, 0.28);
  EXPECT_DOUBLE_EQ(p.uy, -0.96);
  EXPECT_DOUBLE_EQ(p.uz, 0.0);
}

TEST(CaerorGeometryState, ReflectOffSurface) {
  GeometryConstructor input;
  auto left = input.CreateAxisAlignedPlane(Axis::X, 0.0);
  auto right =
      input.CreateAxisAlignedPlane(Axis::X, 2.0, BoundaryCondition::Reflective);
  input.CreateCell(+left && -right, MaterialID{0});

  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();

  Particle p{};
  p.x = 2.0;
  p.ux = 0.6;
  p.uy = 0.8;
  p.cell = CellID{0};

  auto result = state.CrossSurface(p, right.id_);
  EXPECT_EQ(result.cell, p.cell);
  EXPECT_FALSE(result.leaked);
  EXPECT_DOUBLE_EQ(p.ux, -0.6);
  EXPECT_DOUBLE_EQ(p.uy, 0.8);
  EXPECT_DOUBLE_EQ(p.uz, 0.0);
  EXPECT_DOUBLE_EQ(p.x, 2.0);
  EXPECT_DOUBLE_EQ(p.y, 0.0);
  EXPECT_DOUBLE_EQ(p.z, 0.0);
  auto next = state.DistanceToCellBoundary(p);
  EXPECT_DOUBLE_EQ(next.distance, 2.0 / 0.6);
  EXPECT_EQ(next.surface, left.id_);
}

}  // namespace
