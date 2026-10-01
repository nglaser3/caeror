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
        delete[] surface_data_offsets_;

        delete[] rpn_logic_data_;
        delete[] cell_logic_offsets_;
        delete[] cell_mat_fill_ids_;
        delete[] cell_uni_fill_ids_;

        delete[] uni_cell_ids_;
        delete[] universe_offsets_;
      }

      void GenerateSurfacesData(const GeometryConstructor& input){
        const auto& surface_types = input.GetSurfaceTypes();
        size_t surf_data_size = 0;
        surface_data_offsets_ = new size_t[surface_types.size() + 1];
        surface_data_offsets_[0] = 0;
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
          surface_data_offsets_[surf_id + 1] = surface_data_offsets_[surf_id] + data_size;
          surf_data_size += data_size;
        }

        surfaces_data_ = new double[surf_data_size];
        const auto& surfaces = input.GetSurfaces();
        for(size_t surf_id = 0; surf_id < surface_types.size(); surf_id++) {
          auto stype = surface_types[surf_id];
          auto offset = surface_data_offsets_[surf_id];
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
        cell_logic_offsets_ = new size_t[cells.size() + 1];
        cell_logic_offsets_[0] = 0;
        cell_mat_fill_ids_ = new CaerorIndexType[cells.size()];
        cell_uni_fill_ids_ = new CaerorIndexType[cells.size()];
        size_t rpn_token_size = 0;
        for (const auto& cell : cells) {
          const auto cell_id = *(cell->id_);
          cell_uni_fill_ids_[cell_id] = *(cell->universe_fill_);
          cell_mat_fill_ids_[cell_id] = *(cell->material_fill_);
          size_t logic_size = 0;
          if (*(cell->material_fill_) != MAXCaerorIndex) 
            logic_size = cell->region_.logic_.size();
          cell_logic_offsets_[cell_id + 1] = cell_logic_offsets_[cell_id] + logic_size;
          rpn_token_size += logic_size;
        }

        rpn_logic_data_ = new CaerorTokenType[rpn_token_size];
        for (const auto& cell : cells) {
          const auto cell_id = *(cell->id_);
          if (cell_mat_fill_ids_[cell_id] == MAXCaerorIndex) continue;
          auto rpn_offset = cell_logic_offsets_[*(cell->id_)];
          const auto& logic = cell->region_.logic_;
          for (size_t local_index = 0; local_index < logic.size(); local_index++) {
            rpn_logic_data_[rpn_offset + local_index] = *logic[local_index];
          }
        }
      }

      void GenerateUniverseData(const GeometryConstructor& input) {
        const auto& universes = input.GetUniverses();
        universe_offsets_ = new size_t[universes.size() + 1];
        universe_offsets_[0] = 0;
        size_t total_num_uni_cells = 0;
        for (const auto& universe : universes) {
          const auto uni_id = *(universe->id_);
          auto num_cells = universe->cells_.size();
          universe_offsets_[uni_id + 1] = universe_offsets_[uni_id] + num_cells;
          total_num_uni_cells += num_cells;
        }

        uni_cell_ids_ = new CaerorIndexType[total_num_uni_cells];
        for (const auto& universe : universes) {
          const auto uni_id = *(universe->id_);
          const auto offset = universe_offsets_[uni_id];
          for (size_t local_index = 0; local_index < universe->cells_.size(); local_index++) {
            const auto cell_id = *(universe->cells_[local_index]);
            uni_cell_ids_[offset + local_index] = cell_id;
          }
        }
      }

    private:
      /// @brief 1d array for surface data and type (0th index)
      double* surfaces_data_;
      /// @brief 1d array for surface data offsets
      size_t* surface_data_offsets_;

      /// @brief 1d array for region logic
      CaerorTokenType* rpn_logic_data_;
      /// @brief 1d array for cell offsets in rpn logic
      size_t* cell_logic_offsets_;
      /// @brief 1d array mapping cell index to material id
      CaerorIndexType* cell_mat_fill_ids_;
      ///@brief 1d array mapping cell index to universe id
      CaerorIndexType* cell_uni_fill_ids_;

      /// @brief 1d array for cells in universes
      CaerorIndexType* uni_cell_ids_;
      /// @brief 1d array for universe offsets into cells
      size_t* universe_offsets_;
  };
} // namespace caeror
