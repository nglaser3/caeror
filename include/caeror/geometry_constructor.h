#pragma once

#include <vector>
#include <memory>

#include <RAJA/RAJA.hpp>

#include "caeror/raja_layouts.h"
#include "caeror/region.h"
#include "caeror/surface.h"
#include "caeror/volumes.h"

namespace caeror {
  class GeometryConstructor {
    public:
      template <SurfaceType Type, typename... DataTypes> 
      auto CreateSurface(camp::tuple<DataTypes...>&& data) {
        auto id = surfaces_.size();
        auto surf_ptr = std::make_unique<Surface<Type, DataTypes...>>(id, data);
        surface_types_.push_back(surf_ptr->surface_type_);
        surfaces_.push_back(std::move(surf_ptr));
        return static_cast<Surface<Type, DataTypes...>>(*surfaces_.back());
      }

      template <typename FillIDType>
      Cell CreateCell(const Region region, const FillIDType);

      Universe CreateUniverse(std::vector<Cell> cells);
    private:
      std::vector<std::unique_ptr<SurfaceBase>> surfaces_;
      std::vector<SurfaceType> surface_types_;
      std::vector<std::unique_ptr<Cell>> cells_;
      std::vector<std::unique_ptr<Universe>> universes_; 
  };

  template <>
  Cell GeometryConstructor::CreateCell<MaterialID>(const Region region, const MaterialID fill_id) {
    auto id = cells_.size();
    auto cell_ptr = std::make_unique<Cell>(CellID{id}, region);
    cell_ptr->material_fill_ = fill_id;
    cells_.push_back(std::move(cell_ptr));
    return *cells_.back();
  }

  template <>
  Cell GeometryConstructor::CreateCell<Universe>(const Region region, const Universe fill) {
    auto id = cells_.size();
    auto cell_ptr = std::make_unique<Cell>(CellID{id}, region);
    cell_ptr->universe_fill_ = fill.id_;
    cells_.push_back(std::move(cell_ptr));
    return *cells_.back();
  }
} // namespace caeror
