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
/**
 * @brief Class for constructing a geometry input
 */
class GeometryConstructor {
 public:
  /**
   * @brief Create a Plane surface
   * Plane equation: a*x + b*y + c*z + d = 0.
   */
  Plane CreatePlane(double a, double b, double c, double d,
                    BoundaryCondition bc = BoundaryCondition::Transmission) {
    return CreateSurface<SurfaceType::Plane>(camp::make_tuple(a, b, c, d), bc);
  }

  /**
   * @brief Create an axis aligned plane surface
   */
  AxisAlignedPlane CreateAxisAlignedPlane(
      Axis axis, double intercept,
      BoundaryCondition bc = BoundaryCondition::Transmission) {
    return CreateSurface<SurfaceType::AxisAlignedPlane>(
        camp::make_tuple(axis, intercept), bc);
  }

  /**
   * @brief Create a sphere surface
   */
  Sphere CreateSphere(double radius, double x, double y, double z,
                      BoundaryCondition bc = BoundaryCondition::Transmission) {
    return CreateSurface<SurfaceType::Sphere>(camp::make_tuple(radius, x, y, z),
                                              bc);
  }

  /**
   * @brief Create an axis aligned cylinder surface
   */
  AxisAlignedCylinder CreateAxisAlignedCylinder(
      Axis axis, double radius, double center1, double center2,
      BoundaryCondition bc = BoundaryCondition::Transmission) {
    return CreateSurface<SurfaceType::AxisAlignedCylinder>(
        camp::make_tuple(axis, radius, center1, center2), bc);
  }

  /**
   * @brief Helper function for creating a surface instance
   */
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

  /**
   * @brief Create a cell from a region and fill
   */
  template <typename FillIDType>
  Cell CreateCell(const Region region, const FillIDType);

  /**
   * @brief Create a universe from a list of cells
   */
  Universe CreateUniverse(std::vector<Cell> cells);

  /**
   * @brief Getter for all surfaces in geometry
   */
  auto& GetSurfaces() const { return surfaces_; }
  
  /**
   * @brief Getter for all surface types in geometry
   */
  auto& GetSurfaceTypes() const { return surface_types_; }

  /**
   * @brief Getter for all cells in geometry
   */
  auto& GetCells() const { return cells_; }

  /**
   * @brief Getter for all universes in geometry
   */
  auto& GetUniverses() const { return universes_; }

 private:
  /// @brief List of all surfaces in geometry
  std::vector<std::unique_ptr<SurfaceBase>> surfaces_;
  /// @brief List of all surfaces' types in geometry
  std::vector<SurfaceType> surface_types_;
  /// @brief List of all cells 
  std::vector<std::unique_ptr<Cell>> cells_;
  /// @brief List of all universes
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
