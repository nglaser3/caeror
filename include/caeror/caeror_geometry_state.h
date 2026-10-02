#pragma once

#include "caeror/data.h"
#include "caeror/raja_layouts.h"
#include "caeror/stack.h"
#include "caeror/surface.h"

namespace caeror {

struct DistanceResult {
  double distance;
  SurfaceID surface;
};

template <typename RajaResource> struct CaerorGeometryState {
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
  CellID FindCell(const Particle &p) const {
    for (const auto &cell :
         RAJA::range<CellID>(0, logic.get_layout().template size<0>())) {
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
      case RPN_AND:
        stack.Push(stack.Pop() && stack.Pop());
        break;
      case RPN_OR:
        stack.Push(stack.Pop() || stack.Pop());
        break;
      default:
        const auto &s = SurfaceID{token};
        SurfaceDataIndex index{0};
        const auto &s_type = surfaces(s, index++);
        const auto &sense = Sense(s_type, p, &surfaces(s, index));
        stack.Push(sense == SenseResult::Positive);
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
        const auto &s_type = surfaces(s, index++);
        const auto &result = Intersection(s_type, p, &surfaces(s, index));
        if (result.distance < min_distance) {
          surf = s;
          min_distance = result.distance;
        }
        break;
      }
    }
    return {min_distance, surf};
  };

  RAJA_HOST_DEVICE
  CellID CrossSurface(const Particle &p, const SurfaceID &s) const;

private:
  friend class CaerorGeometry;

  double *surfaces_data_;
  CaerorIndexType *surface_cell_ids_;
  size_t *surface_cell_sizes_;

  RPNToken *rpn_logic_data_;
  size_t *rpn_logic_sizes_;
  CaerorIndexType *cell_mat_fill_ids_;
  CaerorIndexType *cell_uni_fill_ids_;

  CaerorIndexType *uni_cell_ids_;
  size_t *uni_cell_sizes_;
};
} // namespace caeror
