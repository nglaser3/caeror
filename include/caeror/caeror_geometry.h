#pragma once

#include <camp/tuple.hpp>
#include <RAJA/RAJA.hpp>

#include "caeror/caeror_geometry_state.h"
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
        delete[] surface_cell_ids_;
        delete[] surface_cell_sizes_;
        delete[] cell_mat_fill_ids_;
        delete[] cell_uni_fill_ids_;

        delete[] uni_cell_ids_;
        delete[] uni_cell_sizes_;
      }

      void GenerateSurfacesData(const GeometryConstructor& input){
        const auto& surface_types = input.GetSurfaceTypes();
        num_surfaces_ = surface_types.size();
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
        num_cells_ = cells.size();
        rpn_logic_sizes_ = new size_t[cells.size()]{};
        cell_mat_fill_ids_ = new CaerorIndexType[cells.size()];
        cell_uni_fill_ids_ = new CaerorIndexType[cells.size()];
        surface_cell_sizes_ = new size_t[num_surfaces_]{};
        auto cell_sizes = new size_t[num_cells_]{};
        for (const auto& cell : cells) {
          const auto cell_id = *(cell->id_);
          cell_uni_fill_ids_[cell_id] = *(cell->universe_fill_);
          cell_mat_fill_ids_[cell_id] = *(cell->material_fill_);
          size_t logic_size = 0;
          if (*(cell->material_fill_) != MAXCaerorIndex) {
            const auto& logic = cell->region_.logic_;
            rpn_logic_sizes_[cell_id] = logic_size;
            assert(cell->region_.stack_size < MAX_STACK_DEPTH);
            for (const auto& token : logic) {
              if(token != RPN_OR && token != RPN_AND && token != RPN_NOT) {
                auto surf_id = static_cast<CaerorIndexType>(token);
                surface_cell_sizes_[surf_id]++;
                cell_sizes[cell_id]++;
              }
            }
          }
        }

        max_logic_size_ = *std::max_element(rpn_logic_sizes_, rpn_logic_sizes_+num_cells_);
        max_cells_ = *std::max_element(surface_cell_sizes_, surface_cell_sizes_ + num_surfaces_);
        max_surfaces_ = *std::max_element(cell_sizes, cell_sizes + num_cells_);
        delete[] cell_sizes;

        rpn_logic_data_ = new RPNToken[max_logic_size_ * cells.size()];
        surface_cell_ids_ = new CaerorIndexType[max_cells_ * num_surfaces_];
        auto current_surf_index = new size_t[num_surfaces_]{};
        for (const auto& cell : cells) {
          const auto cell_id = *(cell->id_);
          if (cell_mat_fill_ids_[cell_id] == MAXCaerorIndex) continue;
          auto rpn_offset = max_logic_size_ * (*(cell->id_));
          const auto& logic = cell->region_.logic_;
          for (size_t local_index = 0; local_index < logic.size(); local_index++) {
            rpn_logic_data_[rpn_offset + local_index] = logic[local_index];
          }
          for (const auto& token : logic) {
            if(token != RPN_OR && token != RPN_AND && token != RPN_NOT) {
              const auto& surf_id = static_cast<CaerorIndexType>(token);
              auto c_index = current_surf_index[surf_id]++;
              auto index = max_cells_ * surf_id + c_index;
              surface_cell_ids_[index] = cell_id;
            }
          }
        }

        delete[] current_surf_index;
      }

      void GenerateUniverseData(const GeometryConstructor& input) {
        const auto& universes = input.GetUniverses();
        num_universes_ = universes.size();
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

      template <typename RajaResource>
      RAJA_HOST_DEVICE
      CaerorGeometryState<RajaResource> GetState() const {
        CaerorGeometryState<RajaResource> state{};
        state.resource = RajaResource{};

        auto surf_size = num_surfaces_ * max_surf_data_size_ * sizeof(double);
        state.surfaces_data_ = state.resource.allocate(surf_size);
        state.resource.memcpy(state.surfaces_data_, surfaces_data_, surf_size);
        state.surfaces = {state.surfaces_data_, num_surfaces_, max_surf_data_size_};

        auto surf_cell_size = max_cells_ * num_surfaces_ * sizeof(CaerorIndexType);
        state.surface_cell_ids_ = state.resource.allocate(surf_cell_size);
        state.resource.memcpy(state.surface_cell_ids, surface_cell_ids_, surf_cell_size);
        state.surface_cells = {state.surface_cell_ids_, num_surfaces_, max_cells_};

        state.surface_cell_sizes_ = state.resource.allocate(num_surfaces_ * sizeof(size_t));
        state.resource.memcpy(state.surface_cell_sizes_, surface_cell_sizes_, num_surfaces_ * sizeof(size_t));
        state.surface_num_cells = {state.surface_cell_sizes_, num_surfaces_};

        auto logic_size = num_cells_ * max_logic_size_ * sizeof(RPNToken);
        state.rpn_logic_data_ = state.resource.allocate(logic_size);
        state.resource.memcpy(state.rpn_logic_data_, rpn_logic_data_, logic_size);
        state.logic = {state.rpn_logic_data_, num_cells_, max_logic_size_};

        state.rpn_logic_sizes_ = state.resource.allocate(num_cells_*sizeof(size_t));
        state.resource.memcpy(state.rpn_logic_sizes_, rpn_logic_sizes_, num_cells_*sizeof(size_t));
        state.logic_sizes = {state.rpn_logic_sizes_, num_cells_};

        auto cell_size = num_cells_ * sizeof(CaerorIndexType);
        state.cell_mat_fill_ids_ = state.resource.allocate(cell_size);
        state.resource.memcpy(state.cell_mat_fill_ids_, cell_mat_fill_ids_, cell_size);
        state.material_fills = {state.cell_mat_fill_ids_, cell_size};

        state.cell_uni_fill_ids_ = state.resource.allocate(cell_size);
        state.resource.memcpy(state.cell_uni_fill_ids_, cell_uni_fill_ids_, cell_size);
        state.universe_fills = {state.cell_uni_fill_ids_, cell_size};

        auto unicell_size = num_universes_ * max_uni_cells_ * sizeof(CaerorIndexType);
        state.uni_cell_ids_ = state.resource.allocate(unicell_size);
        state.resource.memcpy(state.uni_cell_ids_, uni_cell_ids_, unicell_size);
        state.universe_cells = {state.uni_cell_ids_, num_universes_, max_uni_cells_};

        auto uni_size = num_universes_ * sizeof(size_t);
        state.uni_cell_sizes_ = state.resource.allocate(uni_size);
        state.resource.memcpy(state.uni_cell_sizes_, uni_cell_sizes_, uni_size);
        state.universe_sizes = {state.uni_cell_sizes_, num_universes_};
        return state;
      } 
      
    private:
      /// @brief Number of surfaces
      size_t num_surfaces_;
      /// @brief 2d array for surface data and type (0th index)
      double* surfaces_data_;
      /// @brief second dimension size for surface data
      size_t max_surf_data_size_; 

      /// @brief Number of cells
      size_t num_cells_;
      /// @brief 2d array for region logic
      RPNToken* rpn_logic_data_;
      /// @brief 1d array for logic size of cells
      size_t* rpn_logic_sizes_;
      /// @brief second dimension size for rpn logic
      size_t max_logic_size_;
      /// @brief 2d array mapping surfaces to owning cells
      CaerorIndexType* surface_cell_ids_;
      /// @brief 1d array for number of cells that own each surface
      size_t* surface_cell_sizes_;
      /// @brief max number of cell owners for a single surface
      size_t max_cells_;
      /// @brief max number of surfaces on cell
      size_t max_surfaces_;
      /// @brief 1d array mapping cell index to material id
      CaerorIndexType* cell_mat_fill_ids_;
      ///@brief 1d array mapping cell index to universe id
      CaerorIndexType* cell_uni_fill_ids_;

      /// @brief Number of universes
      size_t num_universes_;
      /// @brief 2d array for cells in universes
      CaerorIndexType* uni_cell_ids_;
      /// @brief 1d array for length of universe cells
      size_t* uni_cell_sizes_;
      /// @brief second dimension size for universe cell ids
      size_t max_uni_cells_;
  };
} // namespace caeror
