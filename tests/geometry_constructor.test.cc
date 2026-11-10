#include <vector>

#include <gtest/gtest.h>

#include "caeror/geometry_constructor.h"

namespace {
using namespace caeror;

TEST(GeometryConstructor, StartsEmpty) {
  GeometryConstructor geometry;
  EXPECT_TRUE(geometry.GetSurfaces().empty());
  EXPECT_TRUE(geometry.GetSurfaceTypes().empty());
  EXPECT_TRUE(geometry.GetCells().empty());
  EXPECT_TRUE(geometry.GetUniverses().empty());
}

TEST(GeometryConstructor, MixedSurfaces) {
  GeometryConstructor geometry;
  const auto plane = geometry.CreatePlane(1.0, -2.0, 3.0, -4.0);
  EXPECT_EQ(plane.id_, SurfaceID{0});
  EXPECT_DOUBLE_EQ(camp::get<0>(plane.data_), 1.0);
  EXPECT_DOUBLE_EQ(camp::get<1>(plane.data_), -2.0);
  EXPECT_DOUBLE_EQ(camp::get<2>(plane.data_), 3.0);
  EXPECT_DOUBLE_EQ(camp::get<3>(plane.data_), -4.0);

  const auto aligned = geometry.CreateAxisAlignedPlane(Axis::Y, -5.0);
  EXPECT_EQ(aligned.id_, SurfaceID{1});
  EXPECT_EQ(camp::get<0>(aligned.data_), Axis::Y);
  EXPECT_DOUBLE_EQ(camp::get<1>(aligned.data_), -5.0);

  const auto sphere = geometry.CreateSphere(6.0, 7.0, -8.0, 9.0);
  EXPECT_EQ(sphere.id_, SurfaceID{2});
  EXPECT_DOUBLE_EQ(camp::get<0>(sphere.data_), 6.0);
  EXPECT_DOUBLE_EQ(camp::get<1>(sphere.data_), 7.0);
  EXPECT_DOUBLE_EQ(camp::get<2>(sphere.data_), -8.0);
  EXPECT_DOUBLE_EQ(camp::get<3>(sphere.data_), 9.0);

  const auto cylinder = geometry.CreateAxisAlignedCylinder(Axis::Z, 10.0, -11.0, 12.0);
  EXPECT_EQ(cylinder.id_, SurfaceID{3});
  EXPECT_EQ(camp::get<0>(cylinder.data_), Axis::Z);
  EXPECT_DOUBLE_EQ(camp::get<1>(cylinder.data_), 10.0);
  EXPECT_DOUBLE_EQ(camp::get<2>(cylinder.data_), -11.0);
  EXPECT_DOUBLE_EQ(camp::get<3>(cylinder.data_), 12.0);

  const auto &surfaces = geometry.GetSurfaces();
  ASSERT_EQ(surfaces.size(), 4);

  const std::vector<SurfaceType> expected_types = {
      SurfaceType::Plane, SurfaceType::AxisAlignedPlane,
      SurfaceType::Sphere, SurfaceType::AxisAlignedCylinder};
  EXPECT_EQ(geometry.GetSurfaceTypes(), expected_types);

  for (size_t i = 0; i < surfaces.size(); ++i) 
    EXPECT_EQ(surfaces[i]->id_, SurfaceID{i});
}

TEST(GeometryConstructor, MaterialCells) {
  GeometryConstructor geometry;

  const auto a = geometry.CreateAxisAlignedPlane(Axis::X, -1.0);
  const auto b = geometry.CreateAxisAlignedPlane(Axis::X, 1.0);

  const auto region = +a && -b;
  const auto first = geometry.CreateCell(region, MaterialID{0});
  const auto second = geometry.CreateCell(!region, MaterialID{1});

  const auto &cells = geometry.GetCells();
  ASSERT_EQ(cells.size(), 2);
  ASSERT_NE(cells[0], nullptr);
  ASSERT_NE(cells[1], nullptr);

  const std::vector<RPNToken> inside = {0, 1, RPN_NOT, RPN_AND};
  EXPECT_EQ(first.id_, CellID{0});
  EXPECT_EQ(first.region_.logic_, inside);
  EXPECT_EQ(first.region_.stack_size, 2);
  EXPECT_EQ(first.material_fill_, MaterialID{0});
  EXPECT_EQ(first.universe_fill_, UniverseID{MAXCaerorIndex});

  const std::vector<RPNToken> outside = {0, 1, RPN_NOT, RPN_AND, RPN_NOT};
  EXPECT_EQ(second.id_, CellID{1});
  EXPECT_EQ(second.region_.logic_, outside);
  EXPECT_EQ(second.region_.stack_size, 2);
  EXPECT_EQ(second.material_fill_, MaterialID{1});
  EXPECT_EQ(second.universe_fill_, UniverseID{MAXCaerorIndex});
}

TEST(GeometryConstructor, UniversesAndUniverseFill) {
  GeometryConstructor geometry;

  const auto plane = geometry.CreateAxisAlignedPlane(Axis::X, 0.0);

  const auto first = geometry.CreateCell(+plane, MaterialID{0});
  const auto second = geometry.CreateCell(-plane, MaterialID{1});
  const auto universe = geometry.CreateUniverse({first, second});

  const auto parent_cell = geometry.CreateCell(+plane, universe);
  const auto root_universe = geometry.CreateUniverse({parent_cell});

  const auto &universes = geometry.GetUniverses();
  ASSERT_EQ(universes.size(), 2);
  ASSERT_NE(universes[0], nullptr);
  ASSERT_NE(universes[1], nullptr);

  const std::vector<CellID> expected = {CellID{0}, CellID{1}};
  EXPECT_EQ(universe.id_, UniverseID{0});
  EXPECT_EQ(universe.cells_, expected);

  const std::vector<CellID> expected_outer = {CellID{2}};
  EXPECT_EQ(root_universe.id_, UniverseID{1});
  EXPECT_EQ(universes[1]->id_, UniverseID{1});
  EXPECT_EQ(root_universe.cells_, expected_outer);
  EXPECT_EQ(universes[1]->cells_, expected_outer);

  const auto &cells = geometry.GetCells();
  ASSERT_EQ(cells.size(), 3);
  const std::vector<RPNToken> expected_region = {0};
  EXPECT_EQ(parent_cell.id_, CellID{2});
  EXPECT_EQ(parent_cell.material_fill_, MaterialID{MAXCaerorIndex});
  EXPECT_EQ(parent_cell.universe_fill_, UniverseID{0});
  EXPECT_EQ(parent_cell.region_.logic_, expected_region);
  EXPECT_EQ(parent_cell.region_.stack_size, 1);

  EXPECT_EQ(cells[0]->material_fill_, MaterialID{0});
  EXPECT_EQ(cells[1]->material_fill_, MaterialID{1});
  EXPECT_EQ(cells[0]->universe_fill_, UniverseID{MAXCaerorIndex});
  EXPECT_EQ(cells[1]->universe_fill_, UniverseID{MAXCaerorIndex});
}

} // namespace
