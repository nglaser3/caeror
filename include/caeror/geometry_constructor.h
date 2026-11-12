#pragma once

#include <RAJA/RAJA.hpp>
#include <camp/tuple.hpp>
#include <memory>
#include <utility>
#include <vector>

#include "caeror/raja_layouts.h"
#include "caeror/region.h"
#include "caeror/surface.h"
#include "caeror/volumes.h"

namespace caeror {
class GeometryConstructor {
 public:
  // Plane equation: a*x + b*y + c*z + d = 0.
  Plane CreatePlane(double a, double b, double c, double d,
                    BoundaryCondition bc = BoundaryCondition::Transmission) {
    return CreateSurface<SurfaceType::Plane>(camp::make_tuple(a, b, c, d), bc);
  }

  AxisAlignedPlane CreateAxisAlignedPlane(
      Axis axis, double intercept,
      BoundaryCondition bc = BoundaryCondition::Transmission) {
    return CreateSurface<SurfaceType::AxisAlignedPlane>(
        camp::make_tuple(axis, intercept), bc);
  }

  Sphere CreateSphere(double radius, double x, double y, double z,
                      BoundaryCondition bc = BoundaryCondition::Transmission) {
    return CreateSurface<SurfaceType::Sphere>(camp::make_tuple(radius, x, y, z),
                                              bc);
  }

  AxisAlignedCylinder CreateAxisAlignedCylinder(
      Axis axis, double radius, double center1, double center2,
      BoundaryCondition bc = BoundaryCondition::Transmission) {
    return CreateSurface<SurfaceType::AxisAlignedCylinder>(
        camp::make_tuple(axis, radius, center1, center2), bc);
  }

  template <SurfaceType Type, typename... DataTypes>
  Surface<Type, DataTypes...> CreateSurface(
      camp::tuple<DataTypes...>&& data,
      BoundaryCondition bc = BoundaryCondition::Transmission) {
    const auto id = SurfaceID{surfaces_.size()};
    auto surf_ptr =
        std::make_unique<Surface<Type, DataTypes...>>(id, std::move(data), bc);
    surface_types_.push_back(surf_ptr->surface_type_);
    surfaces_.push_back(std::move(surf_ptr));
    return static_cast<const Surface<Type, DataTypes...>&>(*surfaces_.back());
  }

  template <typename FillIDType>
  Cell CreateCell(const Region region, const FillIDType);

  Universe CreateUniverse(std::vector<Cell> cells);

  auto& GetSurfaces() const { return surfaces_; }
  auto& GetSurfaceTypes() const { return surface_types_; }
  auto& GetCells() const { return cells_; }
  auto& GetUniverses() const { return universes_; }

 private:
  std::vector<std::unique_ptr<SurfaceBase>> surfaces_;
  std::vector<SurfaceType> surface_types_;
  std::vector<std::unique_ptr<Cell>> cells_;
  std::vector<std::unique_ptr<Universe>> universes_;
};

template <>
inline Cell GeometryConstructor::CreateCell<MaterialID>(
    const Region region, const MaterialID fill_id) {
  auto id = cells_.size();
  auto cell_ptr = std::make_unique<Cell>(CellID{id}, region);
  cell_ptr->material_fill_ = fill_id;
  cells_.push_back(std::move(cell_ptr));
  return *cells_.back();
}

template <>
inline Cell GeometryConstructor::CreateCell<Universe>(const Region region,
                                                      const Universe fill) {
  auto id = cells_.size();
  auto cell_ptr = std::make_unique<Cell>(CellID{id}, region);
  cell_ptr->universe_fill_ = fill.id_;
  cells_.push_back(std::move(cell_ptr));
  return *cells_.back();
}
}  // namespace caeror
