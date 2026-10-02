#pragma once

#include "caeror/data.h"
#include "caeror/raja_layouts.h"
#include "caeror/surface.h"

namespace caeror {

  struct DistanceResult{
    double distance;
    SurfaceID surface;
  };

  template <typename RajaResource>
  struct CaerorGeometryState {
    RajaResource resource;

    CaerorView<const double, SurfaceDataLayout> surfaces;

    CaerorView<const CellID, SurfaceCellLayout> surface_cells;

    CaerorView<const size_t, SurfaceLayout> surface_num_cells;

    CaerorView<const RPNToken, RPNTokenLayout> logic;

    CaerorView<const size_t, CellLayout> logic_sizes;

    CaerorView<const MaterialID, CellLayout> material_fills;

    CaerorView<const UniverseID, CellLayout> universe_fills;

    CaerorView<const CellID, UniverseCellLayout> universe_cells;

    CaerorView<const size_t, UniverseLayout> universe_sizes;

    RAJA_HOST_DEVICE
    CellID FindCell(const Particle& p) const {
        for (const auto& cell : RAJA::range<CellID>(0, logic.get_layout().template size<0>())) {
          if (InCell(p, cell)) return cell;
        }
        return CellID{MAXCaerorIndex};
    }

    RAJA_HOST_DEVICE
    bool InCell(const Particle& p, const CellID c) const;

    RAJA_HOST_DEVICE
    DistanceResult DistanceToSurface(const Particle& p) const;

    RAJA_HOST_DEVICE
    void CrossSurface(Particle& p, const SurfaceID& s) const;

    private:
      friend class CaerorGeometry;

      double* surfaces_data_;
      CaerorIndexType* surface_cell_ids_;
      size_t* surface_cell_sizes_;

      RPNToken* rpn_logic_data_;
      size_t* rpn_logic_sizes_;
      CaerorIndexType* cell_mat_fill_ids_;
      CaerorIndexType* cell_uni_fill_ids_;
      
      CaerorIndexType* uni_cell_ids_;
      size_t* uni_cell_sizes_;
  };
} // namespace caeror
