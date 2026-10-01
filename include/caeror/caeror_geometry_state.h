#pragma once

#include "caeror/raja_layouts.h"

namespace caeror {
  template <typename RajaResource>
  struct CaerorGeometryState {
    RajaResource resource;

    using SurfaceView = CaerorView<const double, SurfaceDataLayout>;
    SurfaceView surfaces;

    CaerorView<const RPNToken, RPNTokenLayout> logic;

    CaerorView<const size_t, CellLayout> logic_sizes;

    CaerorView<const CaerorIndexType, CellLayout> material_fills;

    CaerorView<const CaerorIndexType, CellLayout> universe_fills;

    CaerorView<const CaerorIndexType, UniverseCellLayout> universe_cells;

    CaerorView<const size_t, UniverseLayout> universe_sizes;

    private:
      friend class CaerorGeometry;

      double* surfaces_data_;

      RPNToken* rpn_logic_data_;
      size_t* rpn_logic_sizes_;
      CaerorIndexType* cell_mat_fill_ids_;
      CaerorIndexType* cell_uni_fill_ids_;
      
      CaerorIndexType* uni_cell_ids_;
      size_t* uni_cell_sizes_;
  };
  
} // namespace caeror
