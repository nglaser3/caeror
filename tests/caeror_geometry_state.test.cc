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
  // Inside +x, so the union must be true even though both other sides are false.
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
    auto crossed = state.CrossSurface(p, SurfaceID{1});
    EXPECT_EQ(crossed, p.cell);
  }
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
  auto result = state.DistanceToSurface(p);
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
  auto result = state.DistanceToSurface(p);
  EXPECT_DOUBLE_EQ(result.distance, 4.0);
  EXPECT_EQ(result.surface, right.id_);
  p.x = 2.0;
  auto crossed = state.CrossSurface(p, middle.id_);
  EXPECT_EQ(crossed, p.cell);
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
  EXPECT_NEAR(result.distance, 2.0 / 1.4, TOLERANCE);
  EXPECT_EQ(result.surface, plane.id_);
}

TEST(CaerorGeometryState, SphereEntranceFindsFarRoot) {
  GeometryConstructor input;
  auto sphere = input.CreateSphere(2.0, 0.0, 0.0, 0.0);
  input.CreateCell(-sphere, MaterialID{0});
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.x = -2.0;
  p.ux = 1.0;
  p.cell = CellID{0};
  auto result = state.DistanceToSurface(p);
  EXPECT_DOUBLE_EQ(result.distance, 4.0);
  EXPECT_EQ(result.surface, sphere.id_);
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
  auto result = state.DistanceToSurface(p);
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
    auto entering = state.EnteringCell(p, next, SurfaceID{1});
    EXPECT_TRUE(entering);
    auto leaving = state.EnteringCell(p, p.cell, SurfaceID{1});
    EXPECT_FALSE(leaving);
    auto result = state.CrossSurface(p, SurfaceID{1});
    EXPECT_EQ(result, next);
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
  EXPECT_EQ(result, p.cell);
}

TEST(CaerorGeometryState, ExitGeometryHasNoNextCell) {
  auto input = CreateSlabs();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.x = 5.0;
  p.ux = 1.0;
  p.cell = CellID{1};
  auto result = state.CrossSurface(p, SurfaceID{2});
  EXPECT_EQ(result, CellID{MAXCaerorIndex});
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
  auto entering = state.EnteringCell(p, CellID{0}, x.id_);
  EXPECT_FALSE(entering);
  auto result = state.CrossSurface(p, x.id_);
  EXPECT_EQ(result, CellID{1});
}

TEST(CaerorGeometryState, EnteringComplementOfIntersection) {
  GeometryConstructor input;
  auto x = input.CreateAxisAlignedPlane(Axis::X, 0.0);
  auto y = input.CreateAxisAlignedPlane(Axis::Y, 0.0);
  input.CreateCell(!(+x && +y), MaterialID{0});
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  Particle p{};
  p.y = 1.0;
  p.ux = -1.0;
  auto result = state.EnteringCell(p, CellID{0}, x.id_);
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
    EXPECT_EQ(material, i == 0 ? MaterialID{10} : MaterialID{20});
    auto boundary = state.DistanceToSurface(p);
    ASSERT_DOUBLE_EQ(boundary.distance, i == 0 ? 1.0 : 3.0);
    ASSERT_EQ(boundary.surface, SurfaceID{i + 1});
    p.x += boundary.distance * p.ux;
    auto next = state.CrossSurface(p, boundary.surface);
    p.cell = next;
  }
  EXPECT_EQ(p.cell, CellID{MAXCaerorIndex});
  EXPECT_DOUBLE_EQ(p.x, 5.0);
}

TEST(CaerorGeometryState, KobayashiMembership) {
  auto input = CreateKobayashi();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  const double positions[][3] = {{5, 5, 5}, {5, 20, 5}, {20, 55, 5},
                                 {35, 55, 20}, {35, 80, 35}, {50, 80, 50},
                                 {-1, 5, 5}, {61, 5, 5}};
  const CaerorIndexType expected[] = {2, 1, 1, 1, 1, 0, MAXCaerorIndex, MAXCaerorIndex};
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
    auto boundary = state.DistanceToSurface(p);
    ASSERT_DOUBLE_EQ(boundary.distance, distances[i]);
    ASSERT_EQ(boundary.surface, SurfaceID{surfaces[i]});
    p.y += boundary.distance * p.uy;
    auto next = state.CrossSurface(p, boundary.surface);
    p.cell = next;
  }
  EXPECT_EQ(p.cell, CellID{MAXCaerorIndex});
  EXPECT_DOUBLE_EQ(p.y, 100.0);
}

} // namespace
