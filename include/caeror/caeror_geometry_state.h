#pragma once

#include <cmath>
#include <cstddef>

#include "caeror/data.h"
#include "caeror/raja_layouts.h"
#include "caeror/stack.h"
#include "caeror/surface.h"

namespace caeror {

struct DistanceResult {
  double distance;
  SurfaceID surface;
};

struct CaerorGeometryData {
  double *surfaces_data = nullptr;
  CellID *surface_cell_ids = nullptr;
  size_t *surface_cell_sizes = nullptr;
  RPNToken *rpn_logic_data = nullptr;
  size_t *rpn_logic_sizes = nullptr;
  MaterialID *cell_mat_fill_ids = nullptr;
  UniverseID *cell_uni_fill_ids = nullptr;
  CellID *uni_cell_ids = nullptr;
  size_t *uni_cell_sizes = nullptr;
  size_t num_surfaces = 0;
  size_t surface_width = 0;
  size_t max_cells = 0;
  size_t num_cells = 0;
  size_t logic_width = 0;
  size_t num_universes = 0;
  size_t universe_width = 0;
};

struct CaerorGeometryState {
  CaerorView<double, SurfaceDataLayout> surfaces;

  CaerorView<CellID, SurfaceCellLayout> surface_cells;

  CaerorView<size_t, SurfaceLayout> surface_num_cells;

  CaerorView<RPNToken, RPNTokenLayout> logic;

  CaerorView<size_t, CellLayout> logic_sizes;

  CaerorView<MaterialID, CellLayout> material_fills;

  CaerorView<UniverseID, CellLayout> universe_fills;

  CaerorView<CellID, UniverseCellLayout> universe_cells;

  CaerorView<size_t, UniverseLayout> universe_sizes;

  RAJA_HOST_DEVICE
  explicit CaerorGeometryState(const CaerorGeometryData &data)
      : surfaces(data.surfaces_data, data.num_surfaces, data.surface_width),
        surface_cells(data.surface_cell_ids, data.num_surfaces, data.max_cells),
        surface_num_cells(data.surface_cell_sizes, data.num_surfaces),
        logic(data.rpn_logic_data, data.num_cells, data.logic_width),
        logic_sizes(data.rpn_logic_sizes, data.num_cells),
        material_fills(data.cell_mat_fill_ids, data.num_cells),
        universe_fills(data.cell_uni_fill_ids, data.num_cells),
        universe_cells(data.uni_cell_ids, data.num_universes, data.universe_width),
        universe_sizes(data.uni_cell_sizes, data.num_universes) {}

  RAJA_HOST_DEVICE
  CaerorGeometryState(const CaerorGeometryState &) = default;

  CaerorGeometryState &operator=(const CaerorGeometryState &) = delete;

  RAJA_HOST_DEVICE
  CellID FindCell(const Particle &p) const {
    for (const auto &cell :
         RAJA::range<CellID>(0, logic.get_layout().template get_dim_size<0>())) {
      if (InCell(p, cell))
        return cell;
    }
    return CellID{MAXCaerorIndex};
  }

  RAJA_HOST_DEVICE
  bool InCell(const Particle &p, const CellID c) const {
    bool stack_data[MAX_STACK_DEPTH];
    Stack stack(stack_data);
    for (const auto t : RAJA::range<RPNTokenIndex>(0, logic_sizes(c))) {
      const auto &token = logic(c, t);
      switch (token) {
      case RPN_NOT:
        stack.Push(!stack.Pop());
        break;
      case RPN_AND: {
        const bool lhs = stack.Pop();
        const bool rhs = stack.Pop();
        stack.Push(lhs && rhs);
        break;
      }
      case RPN_OR: {
        const bool lhs = stack.Pop();
        const bool rhs = stack.Pop();
        stack.Push(lhs || rhs);
        break;
      }
      default:
        const auto &s = SurfaceID{token};
        SurfaceDataIndex index{0};
        const auto s_type =
            static_cast<SurfaceType>(surfaces(s, index++));
        const auto& sense = DirectionalSense(s_type, p, &surfaces(s, index));
        stack.Push(sense != SenseResult::Negative);
        break;
      }
    }
    return stack.Pop();
  };

  RAJA_HOST_DEVICE
  DistanceResult DistanceToSurface(const Particle &p) const {
    const auto &c = p.cell;
    SurfaceID surf{MAXCaerorIndex};
    double min_distance = INFINITY;
    for (const auto t : RAJA::range<RPNTokenIndex>(0, logic_sizes(c))) {
      const auto &token = logic(c, t);
      switch (token) {
      case RPN_NOT:
      case RPN_AND:
      case RPN_OR:
        break;
      default:
        const auto &s = SurfaceID{token};
        SurfaceDataIndex index{0};
        const auto s_type =
            static_cast<SurfaceType>(surfaces(s, index++));
        const auto &result = Intersection(s_type, p, &surfaces(s, index));
        if (result.valid_ && result.distance_ < min_distance) {
          surf = s;
          min_distance = result.distance_;
        }
        break;
      }
    }
    return {min_distance, surf};
  };

  RAJA_HOST_DEVICE
  CellID CrossSurface(const Particle &p, const SurfaceID &s) const {
    const auto stype =
        static_cast<SurfaceType>(surfaces(s, SurfaceDataIndex{0}));
    if (SenseResult::On != Sense(stype, p, &surfaces(s, SurfaceDataIndex{1})))
      // not actually on a boundary, no cross possible
      return p.cell;

    // loop over cells attached to the surface
    for (const auto lc : RAJA::range<LocalCellID>(0, surface_num_cells(s))) {
      const auto &c = surface_cells(s, lc);
      if (EnteringCell(p, c, s))
        return c;
    }
    return CellID{MAXCaerorIndex};
  };

  RAJA_HOST_DEVICE
  bool EnteringCell(const Particle &p, const CellID &c,
                    const SurfaceID &s) const {
    RPNTokenIndex t{0};
    while (*t < logic_sizes(c)) {
      const auto &token = logic(c, t);
      switch (token) {
      case RPN_NOT:
      case RPN_AND:
      case RPN_OR:
        break;
      default:
        const auto &surf = SurfaceID{token};
        if (s != surf)
          break;
        SurfaceDataIndex index{0};
        const auto s_type =
            static_cast<SurfaceType>(surfaces(s, index++));
        const auto &surf_sense =
            (logic(c, t + 1) == RPN_NOT ? SenseResult::Negative
                                        : SenseResult::Positive);
        const auto going_into_surface =
            AlignedNormal(s_type, p, surf_sense, &surfaces(s, index));
        if (going_into_surface)
          return true;
        break;
      }
      t++;
    }
    return false;
  }
};
} // namespace caeror
