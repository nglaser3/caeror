#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <vector>

#include <RAJA/RAJA.hpp>
#include <camp/tuple.hpp>

#include "caeror/caeror_geometry_state.h"
#include "caeror/geometry_constructor.h"
#include "caeror/raja_layouts.h"
#include "caeror/stack.h"

namespace caeror {

class CaerorGeometry {
public:
  CaerorGeometry(const GeometryConstructor &input) {
    GenerateSurfacesData(input);
    GenerateCellsData(input);
    GenerateUniverseData(input);
  }

  CaerorGeometry(const CaerorGeometry &) = delete;
  CaerorGeometry &operator=(const CaerorGeometry &) = delete;

  ~CaerorGeometry() {
    for (auto &release : release_resource_data_) release();

    delete[] surfaces_data_;

    delete[] rpn_logic_data_;
    delete[] rpn_logic_sizes_;
    delete[] surface_cell_ids_;
    delete[] surface_cell_sizes_;
    delete[] cell_mat_fill_ids_;
    delete[] cell_uni_fill_ids_;

    delete[] uni_cell_ids_;
    delete[] uni_cell_sizes_;
  }

  void GenerateSurfacesData(const GeometryConstructor &input) {
    const auto &surface_types = input.GetSurfaceTypes();
    num_surfaces_ = surface_types.size();
    max_surf_data_size_ = 0;
    for (size_t surf_id = 0; surf_id < surface_types.size(); surf_id++) {
      auto stype = surface_types[surf_id];
      size_t data_size = 0;
      switch (stype) {
      case SurfaceType::Plane:
        data_size = SurfaceDataReq<SurfaceType::Plane>();
        break;
      case SurfaceType::AxisAlignedPlane:
        data_size = SurfaceDataReq<SurfaceType::AxisAlignedPlane>();
        break;
      case SurfaceType::AxisAlignedCylinder:
        data_size = SurfaceDataReq<SurfaceType::AxisAlignedCylinder>();
        break;
      case SurfaceType::Sphere:
        data_size = SurfaceDataReq<SurfaceType::Sphere>();
        break;
      }
      if (data_size > max_surf_data_size_)
        max_surf_data_size_ = data_size;
    }

    const auto &surfaces = input.GetSurfaces();
    surfaces_data_ = new double[max_surf_data_size_ * surfaces.size()]{};
    for (size_t surf_id = 0; surf_id < surface_types.size(); surf_id++) {
      auto stype = surface_types[surf_id];
      auto offset = max_surf_data_size_ * surf_id;
      surfaces_data_[offset + SURFACE_TYPE_INDEX] = static_cast<double>(stype);

      const auto &surf = *surfaces[surf_id];
      surfaces_data_[offset + SURFACE_BC_INDEX] = static_cast<double>(surf.bc_);
      switch (stype) {
      case SurfaceType::Plane: {
        FillSurfaceData<SurfaceType::Plane, PlaneData::A>(surf, offset);
        FillSurfaceData<SurfaceType::Plane, PlaneData::B>(surf, offset);
        FillSurfaceData<SurfaceType::Plane, PlaneData::C>(surf, offset);
        FillSurfaceData<SurfaceType::Plane, PlaneData::D>(surf, offset);
        break;
      }
      case SurfaceType::AxisAlignedPlane: {
        FillSurfaceData<SurfaceType::AxisAlignedPlane,
                        AxisAlignedPlaneData::Axis>(surf, offset);
        FillSurfaceData<SurfaceType::AxisAlignedPlane,
                        AxisAlignedPlaneData::Intercept>(surf, offset);
        break;
      }
      case SurfaceType::AxisAlignedCylinder: {
        FillSurfaceData<SurfaceType::AxisAlignedCylinder,
                        AxisAlignedCylinderData::Axis>(surf, offset);
        FillSurfaceData<SurfaceType::AxisAlignedCylinder,
                        AxisAlignedCylinderData::Radius>(surf, offset);
        FillSurfaceData<SurfaceType::AxisAlignedCylinder,
                        AxisAlignedCylinderData::Center1>(surf, offset);
        FillSurfaceData<SurfaceType::AxisAlignedCylinder,
                        AxisAlignedCylinderData::Center2>(surf, offset);
        break;
      }
      case SurfaceType::Sphere: {
        FillSurfaceData<SurfaceType::Sphere, SphereData::Radius>(surf, offset);
        FillSurfaceData<SurfaceType::Sphere, SphereData::XCenter>(surf, offset);
        FillSurfaceData<SurfaceType::Sphere, SphereData::YCenter>(surf, offset);
        FillSurfaceData<SurfaceType::Sphere, SphereData::ZCenter>(surf, offset);
        break;
      }
      }
    }
  }

  template <SurfaceType Type, auto Index>
  void FillSurfaceData(const SurfaceBase &s, size_t offset) {
    const auto &surface = CastSurface<Type>(s);
    constexpr size_t index = static_cast<size_t>(Index);
    surfaces_data_[offset + SURFACE_DATA_OFFSET + index] =
        static_cast<double>(camp::get<index>(surface.data_));
  }

  void GenerateCellsData(const GeometryConstructor &input) {
    const auto &cells = input.GetCells();
    num_cells_ = cells.size();
    rpn_logic_sizes_ = new size_t[cells.size()]{};
    cell_mat_fill_ids_ = new MaterialID[cells.size()];
    cell_uni_fill_ids_ = new UniverseID[cells.size()];
    surface_cell_sizes_ = new size_t[num_surfaces_]{};
    std::vector<std::vector<CaerorIndexType>> cell_surfaces(num_cells_);
    max_logic_size_ = 0;
    max_surfaces_ = 0;
    max_cells_ = 0;
    for (const auto &cell : cells) {
      const auto cell_id = *cell->id_;
      cell_uni_fill_ids_[cell_id] = cell->universe_fill_;
      cell_mat_fill_ids_[cell_id] = cell->material_fill_;
      const auto &logic = cell->region_.logic_;
      rpn_logic_sizes_[cell_id] = logic.size();
      max_logic_size_ = std::max(max_logic_size_, logic.size());
      assert(cell->region_.stack_size <= MAX_STACK_DEPTH);
      auto &surfaces = cell_surfaces[cell_id];
      for (auto token : logic) {
        if (token != RPN_OR && token != RPN_AND && token != RPN_NOT) {
          assert(token < num_surfaces_);
          surfaces.push_back(token);
        }
      }
      std::sort(surfaces.begin(), surfaces.end());
      surfaces.erase(std::unique(surfaces.begin(), surfaces.end()), surfaces.end());
      max_surfaces_ = std::max(max_surfaces_, surfaces.size());
      for (auto surface : surfaces)
        max_cells_ = std::max(max_cells_, ++surface_cell_sizes_[surface]);
    }

    rpn_logic_data_ = new RPNToken[max_logic_size_ * num_cells_]{};
    surface_cell_ids_ = new CellID[max_cells_ * num_surfaces_]{};
    std::vector<size_t> current_surf_index(num_surfaces_, 0);
    for (const auto &cell : cells) {
      const auto cell_id = *cell->id_;
      const auto &logic = cell->region_.logic_;
      std::copy(logic.begin(), logic.end(),
                rpn_logic_data_ + max_logic_size_ * cell_id);
      for (auto surface : cell_surfaces[cell_id]) {
        const auto local_index = current_surf_index[surface]++;
        surface_cell_ids_[max_cells_ * surface + local_index] = cell->id_;
      }
    }
  }

  void GenerateUniverseData(const GeometryConstructor &input) {
    const auto &universes = input.GetUniverses();
    num_universes_ = universes.size();
    uni_cell_sizes_ = new size_t[universes.size()]{};
    max_uni_cells_ = 0;
    for (const auto &universe : universes) {
      const auto uni_id = *(universe->id_);
      auto num_cells = universe->cells_.size();
      uni_cell_sizes_[uni_id] = num_cells;
      if (num_cells > max_uni_cells_)
        max_uni_cells_ = num_cells;
    }

    uni_cell_ids_ = new CellID[max_uni_cells_ * universes.size()]{};
    for (const auto &universe : universes) {
      const auto uni_id = *(universe->id_);
      const auto offset = max_uni_cells_ * uni_id;
      for (size_t local_index = 0; local_index < universe->cells_.size();
           local_index++) {
        const auto cell_id = universe->cells_[local_index];
        uni_cell_ids_[offset + local_index] = cell_id;
      }
    }
  }

  template <typename RajaResource>
  CaerorGeometryState GetState(RajaResource resource = {}) & {
    CaerorGeometryData data{};
    data.num_surfaces = num_surfaces_;
    data.surface_width = max_surf_data_size_;
    data.max_cells = max_cells_;
    data.num_cells = num_cells_;
    data.logic_width = max_logic_size_;
    data.num_universes = num_universes_;
    data.universe_width = max_uni_cells_;
    auto copy = [&]<typename T>(T *&destination, const T *source, size_t count) {
      destination = resource.template allocate<T>(count);
      release_resource_data_.emplace_back([resource, destination]() mutable {
        resource.deallocate(destination);
      });
      resource.memcpy(destination, source, count * sizeof(T));
    };
    copy(data.surfaces_data, surfaces_data_, num_surfaces_ * max_surf_data_size_);
    copy(data.surface_cell_ids, surface_cell_ids_, num_surfaces_ * max_cells_);
    copy(data.surface_cell_sizes, surface_cell_sizes_, num_surfaces_);
    copy(data.rpn_logic_data, rpn_logic_data_, num_cells_ * max_logic_size_);
    copy(data.rpn_logic_sizes, rpn_logic_sizes_, num_cells_);
    copy(data.cell_mat_fill_ids, cell_mat_fill_ids_, num_cells_);
    copy(data.cell_uni_fill_ids, cell_uni_fill_ids_, num_cells_);
    copy(data.uni_cell_ids, uni_cell_ids_, num_universes_ * max_uni_cells_);
    copy(data.uni_cell_sizes, uni_cell_sizes_, num_universes_);
    resource.wait();
    return CaerorGeometryState(data);
  }

  template <typename RajaResource>
  CaerorGeometryState GetState(RajaResource resource = {}) && = delete;

private:
  /// @brief Host cleanup functions for arrays allocated through GetState
  std::vector<std::function<void()>> release_resource_data_;

  /// @brief Number of surfaces
  size_t num_surfaces_;
  /// @brief 2d array for surface data and type (0th index)
  double *surfaces_data_;
  /// @brief second dimension size for surface data
  size_t max_surf_data_size_;

  /// @brief Number of cells
  size_t num_cells_;
  /// @brief 2d array for region logic
  RPNToken *rpn_logic_data_;
  /// @brief 1d array for logic size of cells
  size_t *rpn_logic_sizes_;
  /// @brief second dimension size for rpn logic
  size_t max_logic_size_;
  /// @brief 2d array mapping surfaces to owning cells
  CellID *surface_cell_ids_;
  /// @brief 1d array for number of cells that own each surface
  size_t *surface_cell_sizes_;
  /// @brief max number of cell owners for a single surface
  size_t max_cells_;
  /// @brief max number of surfaces on cell
  size_t max_surfaces_;
  /// @brief 1d array mapping cell index to material id
  MaterialID *cell_mat_fill_ids_;
  ///@brief 1d array mapping cell index to universe id
  UniverseID *cell_uni_fill_ids_;

  /// @brief Number of universes
  size_t num_universes_;
  /// @brief 2d array for cells in universes
  CellID *uni_cell_ids_;
  /// @brief 1d array for length of universe cells
  size_t *uni_cell_sizes_;
  /// @brief second dimension size for universe cell ids
  size_t max_uni_cells_;
};

} // namespace caeror
