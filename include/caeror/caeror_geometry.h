#pragma once

#include <camp/tuple.hpp>
#include <RAJA/RAJA.hpp>

#include "caeror/geometry_constructor.h"
#include "caeror/raja_layouts.h"

namespace caeror {

  class CaerorGeometry {
    public:
      CaerorGeometry(const GeometryConstructor& input) {
        GenerateSurfacesData(input);
        GenerateCellsData(input);
        GenerateUniverseData(input);
      }

      ~CaerorGeometry() {
        delete[] surfaces_data_;

        delete[] rpn_logic_data_;
        delete[] rpn_logic_sizes_;
        delete[] cell_mat_fill_ids_;
        delete[] cell_uni_fill_ids_;

        delete[] uni_cell_ids_;
        delete[] uni_cell_sizes_;
      }

      void GenerateSurfacesData(const GeometryConstructor& input){
        const auto& surface_types = input.GetSurfaceTypes();
        max_surf_data_size_ = 0;
        for (size_t surf_id = 0; surf_id < surface_types.size(); surf_id++) {
          auto stype = surface_types[surf_id];
          size_t data_size;
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
          if (data_size > max_surf_data_size_) max_surf_data_size_ = data_size;
        }

        const auto& surfaces = input.GetSurfaces();
        surfaces_data_ = new double[max_surf_data_size_ * surfaces.size()];
        for(size_t surf_id = 0; surf_id < surface_types.size(); surf_id++) {
          auto stype = surface_types[surf_id];
          auto offset = max_surf_data_size_ * surf_id;
          surfaces_data_[offset] = static_cast<double>(stype);

          const auto& surf = *surfaces[surf_id];
          switch (stype) {
            case SurfaceType::Plane: {
              FillSurfaceData<SurfaceType::Plane, PlaneData::A>(surf, offset); 
              FillSurfaceData<SurfaceType::Plane, PlaneData::B>(surf, offset); 
              FillSurfaceData<SurfaceType::Plane, PlaneData::C>(surf, offset); 
              FillSurfaceData<SurfaceType::Plane, PlaneData::D>(surf, offset); 
              break;
            }
            case SurfaceType::AxisAlignedPlane: {
              FillSurfaceData<SurfaceType::AxisAlignedPlane, AxisAlignedPlaneData::Axis>(surf, offset); 
              FillSurfaceData<SurfaceType::AxisAlignedPlane, AxisAlignedPlaneData::Intercept>(surf, offset); 
              break;
            }
            case SurfaceType::AxisAlignedCylinder: {
              FillSurfaceData<SurfaceType::AxisAlignedCylinder, AxisAlignedCylinderData::Axis>(surf, offset); 
              FillSurfaceData<SurfaceType::AxisAlignedCylinder, AxisAlignedCylinderData::Radius>(surf, offset); 
              FillSurfaceData<SurfaceType::AxisAlignedCylinder, AxisAlignedCylinderData::Center1>(surf, offset); 
              FillSurfaceData<SurfaceType::AxisAlignedCylinder, AxisAlignedCylinderData::Center2>(surf, offset); 
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
      void FillSurfaceData(const SurfaceBase& s, size_t offset) {
        const auto& surface = CastSurface<Type>(s);
        constexpr size_t index = static_cast<size_t>(Index);
        surfaces_data_[offset + 1 + index] = static_cast<double>(camp::get<index>(surface.data_));
      }

      void GenerateCellsData(const GeometryConstructor& input) {
        const auto& cells = input.GetCells();
        rpn_logic_sizes_ = new size_t[cells.size()];
        cell_mat_fill_ids_ = new CaerorIndexType[cells.size()];
        cell_uni_fill_ids_ = new CaerorIndexType[cells.size()];
        max_logic_size_ = 0;
        for (const auto& cell : cells) {
          const auto cell_id = *(cell->id_);
          cell_uni_fill_ids_[cell_id] = *(cell->universe_fill_);
          cell_mat_fill_ids_[cell_id] = *(cell->material_fill_);
          size_t logic_size = 0;
          if (*(cell->material_fill_) != MAXCaerorIndex) 
            logic_size = cell->region_.logic_.size();
          rpn_logic_sizes_[cell_id] = logic_size;
          if (logic_size > max_logic_size_) max_logic_size_ = logic_size;
        }

        rpn_logic_data_ = new CaerorTokenType[max_logic_size_ * cells.size()];
        for (const auto& cell : cells) {
          const auto cell_id = *(cell->id_);
          if (cell_mat_fill_ids_[cell_id] == MAXCaerorIndex) continue;
          auto rpn_offset = max_logic_size_ * (*(cell->id_));
          const auto& logic = cell->region_.logic_;
          for (size_t local_index = 0; local_index < logic.size(); local_index++) {
            rpn_logic_data_[rpn_offset + local_index] = *logic[local_index];
          }
        }
      }

      void GenerateUniverseData(const GeometryConstructor& input) {
        const auto& universes = input.GetUniverses();
        uni_cell_sizes_ = new size_t[universes.size()];
        max_uni_cells_ = 0;
        for (const auto& universe : universes) {
          const auto uni_id = *(universe->id_);
          auto num_cells = universe->cells_.size();
          uni_cell_sizes_[uni_id] = num_cells;
          if (num_cells > max_uni_cells_) max_uni_cells_ = num_cells;
        }

        uni_cell_ids_ = new CaerorIndexType[max_uni_cells_ * universes.size()];
        for (const auto& universe : universes) {
          const auto uni_id = *(universe->id_);
          const auto offset = max_uni_cells_ * uni_id;
          for (size_t local_index = 0; local_index < universe->cells_.size(); local_index++) {
            const auto cell_id = *(universe->cells_[local_index]);
            uni_cell_ids_[offset + local_index] = cell_id;
          }
        }
      }

    private:
      /// @brief 2d array for surface data and type (0th index)
      double* surfaces_data_;
      /// @brief second dimension size for surface data
      size_t max_surf_data_size_; 

      /// @brief 2d array for region logic
      CaerorTokenType* rpn_logic_data_;
      /// @brief 1d array for logic size of cells
      size_t* rpn_logic_sizes_;
      /// @brief second dimension size for rpn logic
      size_t max_logic_size_;
      /// @brief 1d array mapping cell index to material id
      CaerorIndexType* cell_mat_fill_ids_;
      ///@brief 1d array mapping cell index to universe id
      CaerorIndexType* cell_uni_fill_ids_;

      /// @brief 2d array for cells in universes
      CaerorIndexType* uni_cell_ids_;
      /// @brief 1d array for length of universe cells
      size_t* uni_cell_sizes_;
      /// @brief second dimension size for universe cell ids
      size_t max_uni_cells_;
  };
} // namespace caeror
