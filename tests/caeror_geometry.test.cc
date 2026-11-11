#include <algorithm>
#include <vector>
#include <type_traits>

#include <gtest/gtest.h>

#include "caeror/caeror_geometry.h"
#include "caeror/caeror_geometry_state.h"
#include "sample_constructors/kobayashi.h"

namespace {
using namespace caeror;

TEST(CaerorGeometry, KobayashiSurfaces) {
  auto input = CreateKobayashi();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();

  ASSERT_EQ(state.surfaces.get_layout().get_dim_size<0>(), 15);
  ASSERT_EQ(state.surfaces.get_layout().get_dim_size<1>(), 4);
  const double intercepts[] = {
      0, 10, 30, 40, 60, 0, 10, 50, 60, 100, 0, 10, 30, 40, 60};
  for (CaerorIndexType i = 0; i < 15; ++i) {
    const SurfaceID id{i};
    auto type = state.surfaces(id, SurfaceDataIndex{0});
    EXPECT_EQ(type, static_cast<double>(SurfaceType::AxisAlignedPlane));
    auto axis = state.surfaces(id, SurfaceDataIndex{1});
    EXPECT_EQ(axis, static_cast<double>(i / 5)); // kinda hacky, but x was made first then y then z lol
    auto intercept = state.surfaces(id, SurfaceDataIndex{2});
    EXPECT_DOUBLE_EQ(intercept, intercepts[i]);
  }
}

TEST(CaerorGeometry, KobayashiCellLogic) {
  auto input = CreateKobayashi();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();

  ASSERT_EQ(state.logic.get_layout().get_dim_size<0>(), 3);
  ASSERT_EQ(state.logic.get_layout().get_dim_size<1>(), 91);

  ASSERT_EQ(state.logic_sizes.get_layout().get_dim_size<0>(), 3);
  const size_t lengths[] = {91, 59, 14};
  for (CaerorIndexType i = 0; i < 3; ++i) {
    const CellID id{i};
    auto size = state.logic_sizes(id);
    ASSERT_EQ(size, lengths[i]);
    const auto &expected = input.GetCells()[i]->region_.logic_;
    ASSERT_EQ(size, expected.size());
    for (CaerorIndexType j = 0; j < size; ++j) {
      auto token = state.logic(id, RPNTokenIndex{j});
      EXPECT_EQ(token, expected[j]);
    }
  }
}

TEST(CaerorGeometry, KobayashiMatsAndUniverses) {
  auto input = CreateKobayashi();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();

  ASSERT_EQ(state.material_fills.get_layout().get_dim_size<0>(), 3);
  ASSERT_EQ(state.universe_fills.get_layout().get_dim_size<0>(), 3);
  ASSERT_EQ(state.universe_sizes.get_layout().get_dim_size<0>(), 1);
  ASSERT_EQ(state.universe_cells.get_layout().get_dim_size<0>(), 1);
  ASSERT_EQ(state.universe_cells.get_layout().get_dim_size<1>(), 3);

  auto size = state.universe_sizes(UniverseID{0});
  ASSERT_EQ(size, 3);

  for (CaerorIndexType i = 0; i < 3; ++i) {
    auto material = state.material_fills(CellID{i});
    EXPECT_EQ(material, MaterialID{i});
    auto universe = state.universe_fills(CellID{i});
    EXPECT_EQ(universe, UniverseID{MAXCaerorIndex});
    auto cell = state.universe_cells(UniverseID{0}, LocalCellID{i});
    EXPECT_EQ(cell, CellID{i});
  }
}

TEST(CaerorGeometry, KobayashiSurfaceCells) {
  auto input = CreateKobayashi();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();

  ASSERT_EQ(state.surface_cells.get_layout().get_dim_size<0>(), 15u);
  ASSERT_EQ(state.surface_num_cells.get_layout().get_dim_size<0>(), 15u);
  const std::vector<CaerorIndexType> expected[] = {
      {0, 1, 2},
      {0, 1, 2},
      {0, 1},
      {0, 1}, 
      {0},
      {0, 2},
      {0, 1, 2},
      {0, 1},
      {0, 1},
      {0, 1},
      {0, 1, 2},
      {0, 1, 2},
      {0, 1},
      {0, 1},
      {0}
    };
  for (CaerorIndexType i = 0; i < 15; ++i) {
    const SurfaceID id{i};
    auto size = state.surface_num_cells(id);
    ASSERT_LE(size, state.surface_cells.get_layout().get_dim_size<1>());
    std::vector<CaerorIndexType> cells;
    for (CaerorIndexType j = 0; j < size; ++j) {
      auto cell = state.surface_cells(id, LocalCellID{j});
      cells.push_back(*cell);
    }
    std::sort(cells.begin(), cells.end());
    EXPECT_EQ(cells, expected[i]);
  }
}

TEST(CaerorGeometry, EmptyGeometry) {
  GeometryConstructor input;
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  EXPECT_EQ(state.surfaces.get_layout().get_dim_size<0>(), 0);
  EXPECT_EQ(state.logic.get_layout().get_dim_size<0>(), 0);
  EXPECT_EQ(state.universe_sizes.get_layout().get_dim_size<0>(), 0);
}

// ============================================================================
// Memory stuff
// Created with help of GPT-6 Astra
// ============================================================================
TEST(CaerorGeometry, GeometryOutlivesInput) {
  auto geometry = [] {
    auto input = CreateKobayashi();
    return CaerorGeometry(input);
  }();
  auto state = geometry.GetState<RAJA::resources::Host>();
  EXPECT_EQ(state.material_fills(CellID{1}), MaterialID{1});
}

// Helper struct for faking a resource, used to count allocations to verify geometry owns their lifetime.
struct CountingHostResource {
  inline static int allocations = 0;
  RAJA::resources::Host host;

  template <typename T> T *allocate(size_t count) {
    auto *pointer = host.allocate<T>(count);
    ++allocations;
    return pointer;
  }
  void deallocate(void *pointer) {
    host.deallocate(pointer);
    --allocations;
  }
  void memcpy(void *destination, const void *source, size_t bytes) {
    host.memcpy(destination, source, bytes);
  }
  void wait() {host.wait();}
};

TEST(CaerorGeometry, GeometryReleasesAllocations) {
  auto input = CreateKobayashi();
  {
    CaerorGeometry geometry(input);
    geometry.GetState<CountingHostResource>();
    EXPECT_EQ(CountingHostResource::allocations, 9);
  }
  EXPECT_EQ(CountingHostResource::allocations, 0);
}

TEST(CaerorGeometry, AdditionalStatesDontCorrupt) {
  auto input = CreateKobayashi();
  {
    CaerorGeometry geometry(input);
    auto first = geometry.GetState(CountingHostResource{});
    auto second = geometry.GetState(CountingHostResource{});
    EXPECT_EQ(CountingHostResource::allocations, 18);
    EXPECT_NE(first.surfaces.get_data(), second.surfaces.get_data());
    auto material = first.material_fills(CellID{2});
    EXPECT_EQ(material, MaterialID{2});
    material = second.material_fills(CellID{2});
    EXPECT_EQ(material, MaterialID{2});

  }
  EXPECT_EQ(CountingHostResource::allocations, 0);
}

TEST(CaerorGeometry, CopiedStateSharesStorage) {
  auto input = CreateKobayashi();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  {
    auto copy = state;
    EXPECT_EQ(copy.surfaces.get_data(), state.surfaces.get_data());
    auto material = copy.material_fills(CellID{2});
    EXPECT_EQ(material, MaterialID{2});
  }
  auto material = state.material_fills(CellID{2});
  EXPECT_EQ(material, MaterialID{2});
}

TEST(CaerorGeometry, CaptureStateInKernel) {
  auto input = CreateKobayashi();
  CaerorGeometry geometry(input);
  auto state = geometry.GetState<RAJA::resources::Host>();
  CaerorIndexType output[3]{};
  auto *results = output;
  RAJA::forall<RAJA::seq_exec>(RAJA::RangeSegment(0, 3), [=](RAJA::Index_type i) {
    auto local = state;
    results[i] = *local.material_fills(CellID{static_cast<CaerorIndexType>(i)});
    results[i] += *state.universe_cells(UniverseID{0}, LocalCellID{0});
  });
  for (CaerorIndexType i = 0; i < 3; ++i) EXPECT_EQ(output[i], i);
}
} // namespace
