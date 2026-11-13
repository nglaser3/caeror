#pragma once

#include <cmath>
#include <cstddef>

#include "caeror/data.h"
#include "caeror/raja_layouts.h"
#include "caeror/stack.h"
#include "caeror/surface.h"

namespace caeror {

/// @brief Container for storing the result from a distance calculation
struct DistanceResult {
  /// @brief Distance to the next surface
  double distance;
  /// @brief Surface hit at distance
  SurfaceID surface;
};

/// @brief Container for storing the result from a surface crossing
struct SurfaceCrossingResult {
  /// @brief Cell the incident particle is in after the surface crossing event
  CellID cell;
  /// @brief Whether the incident particle leaked from the system
  bool leaked;
};

/// @cond
/// @brief Helper struct for constructing CaerorGeometryState
struct CaerorGeometryData {
  double* surfaces_data = nullptr;
  CellID* surface_cell_ids = nullptr;
  size_t* surface_cell_sizes = nullptr;
  RPNToken* rpn_logic_data = nullptr;
  size_t* rpn_logic_sizes = nullptr;
  MaterialID* cell_mat_fill_ids = nullptr;
  UniverseID* cell_uni_fill_ids = nullptr;
  CellID* uni_cell_ids = nullptr;
  size_t* uni_cell_sizes = nullptr;
  size_t num_surfaces = 0;
  size_t surface_width = 0;
  size_t max_cells = 0;
  size_t num_cells = 0;
  size_t logic_width = 0;
  size_t num_universes = 0;
  size_t universe_width = 0;
};
/// @endcond

/**
 * @brief queryable state of CaerorGeometry, where ray-tracing operations exist
 */
struct CaerorGeometryState {
  /// @brief Array containing surface data
  CaerorView<double, SurfaceDataLayout> surfaces;

  /// @brief Array containing the surface - cell mapping
  CaerorView<CellID, SurfaceCellLayout> surface_cells;

  /// @brief Array containing the number of cells each surface belongs to
  CaerorView<size_t, SurfaceLayout> surface_num_cells;

  /// @brief Array containing the region logic for each region
  CaerorView<RPNToken, RPNTokenLayout> logic;

  /// @brief Array containing the size of the logic data for each cell
  CaerorView<size_t, CellLayout> logic_sizes;

  /// @brief Array containing the material IDs for each cell
  CaerorView<MaterialID, CellLayout> material_fills;

  /// @brief Array containing the universe IDs for each cell
  CaerorView<UniverseID, CellLayout> universe_fills;

  /// @brief Array containing the cells within each universe
  CaerorView<CellID, UniverseCellLayout> universe_cells;

  /// @brief Array containing the number of cells in each universe
  CaerorView<size_t, UniverseLayout> universe_sizes;

  /**
   * @brief Construct a CaerorGeometryState from a data instance
   */
  RAJA_HOST_DEVICE
  explicit CaerorGeometryState(const CaerorGeometryData& data)
      : surfaces(data.surfaces_data, data.num_surfaces, data.surface_width),
        surface_cells(data.surface_cell_ids, data.num_surfaces, data.max_cells),
        surface_num_cells(data.surface_cell_sizes, data.num_surfaces),
        logic(data.rpn_logic_data, data.num_cells, data.logic_width),
        logic_sizes(data.rpn_logic_sizes, data.num_cells),
        material_fills(data.cell_mat_fill_ids, data.num_cells),
        universe_fills(data.cell_uni_fill_ids, data.num_cells),
        universe_cells(data.uni_cell_ids, data.num_universes,
                       data.universe_width),
        universe_sizes(data.uni_cell_sizes, data.num_universes) {}

  /**
   * @brief Copy constructor
   */
  RAJA_HOST_DEVICE
  CaerorGeometryState(const CaerorGeometryState&) = default;

  /**
   * @brief Copy assignment
   */
  CaerorGeometryState& operator=(const CaerorGeometryState&) = delete;

  /**
   * @brief Routine to find an incident particle in the geometry
   */
  RAJA_HOST_DEVICE
  CellID FindCell(const Particle& p) const {
    for (const auto& cell : RAJA::range<CellID>(
             0, logic.get_layout().template get_dim_size<0>())) {
      if (InCell(p, cell)) return cell;
    }
    return CellID{MAXCaerorIndex};
  }

  /**
   * @brief Routine to determine if an incident particle is in the specified
   * cell
   */
  RAJA_HOST_DEVICE
  bool InCell(const Particle& p, const CellID c) const {
    bool stack_data[MAX_STACK_DEPTH];
    Stack stack(stack_data);
    for (const auto t : RAJA::range<RPNTokenIndex>(0, logic_sizes(c))) {
      const auto& token = logic(c, t);
      switch (token) {
        case RPN_NOT:
          stack.Push(!stack.Pop());
          break;
        case RPN_AND: {
          const bool rhs = stack.Pop();
          const bool lhs = stack.Pop();
          stack.Push(lhs && rhs);
          break;
        }
        case RPN_OR: {
          const bool rhs = stack.Pop();
          const bool lhs = stack.Pop();
          stack.Push(lhs || rhs);
          break;
        }
        default:
          const auto& s = SurfaceID{token};
          const auto s_type = static_cast<SurfaceType>(
              surfaces(s, SurfaceDataIndex{SURFACE_TYPE_INDEX}));
          const auto& sense = DirectionalSense(
              s_type, p, &surfaces(s, SurfaceDataIndex{SURFACE_DATA_OFFSET}));
          if (sense == SenseResult::On) return InCellTangent(p, c);
          stack.Push(sense == SenseResult::Positive);
          break;
      }
    }
    return stack.Pop();
  };

  /**
   * @brief Routine to determine the distance to the next surface (including
   * interior surfaces)
   */
  RAJA_HOST_DEVICE
  DistanceResult DistanceToSurface(const Particle& p) const {
    const auto& c = p.cell;
    SurfaceID surf{MAXCaerorIndex};
    double min_distance = INFINITY;
    for (const auto t : RAJA::range<RPNTokenIndex>(0, logic_sizes(c))) {
      const auto& token = logic(c, t);
      switch (token) {
        case RPN_NOT:
        case RPN_AND:
        case RPN_OR:
          break;
        default:
          const auto& s = SurfaceID{token};
          const SurfaceDataIndex index{SURFACE_DATA_OFFSET};
          const auto s_type = static_cast<SurfaceType>(
              surfaces(s, SurfaceDataIndex{SURFACE_TYPE_INDEX}));
          auto result = Intersection<false>(s_type, p, &surfaces(s, index));
          if (result.valid_ && result.distance_ == 0.0 && InCell(p, c))
            result = Intersection<true>(s_type, p, &surfaces(s, index));
          if (result.valid_ && result.distance_ < min_distance) {
            surf = s;
            min_distance = result.distance_;
          }
          break;
      }
    }
    return {min_distance, surf};
  };

  /**
   * @brief Routine to find the distance to the nearest surface (not including
   * interior surfaces). This is used to find the distance needed to exit a
   * cell.
   */
  RAJA_HOST_DEVICE
  DistanceResult DistanceToCellBoundary(const Particle& p) const {
    Particle temp = p;
    double distance = 0.0;
    while (true) {
      const auto hit = DistanceToSurface(temp);
      if (hit.surface == SurfaceID{MAXCaerorIndex}) return hit;
      temp.x += hit.distance * temp.ux;
      temp.y += hit.distance * temp.uy;
      temp.z += hit.distance * temp.uz;
      distance += hit.distance;
      if (!InCell(temp, p.cell)) return {distance, hit.surface};
    }
  }

  /**
   * @brief Routine for a particle to cross the specified surface
   */
  RAJA_HOST_DEVICE
  SurfaceCrossingResult CrossSurface(Particle& p, const SurfaceID& s) const {
    const auto stype = static_cast<SurfaceType>(
        surfaces(s, SurfaceDataIndex{SURFACE_TYPE_INDEX}));
    if (SenseResult::On !=
        Sense(stype, p, &surfaces(s, SurfaceDataIndex{SURFACE_DATA_OFFSET})))
      // not actually on a boundary, no cross possible
      return {p.cell, false};

    // if particle is still in the current cell, crossing the surface doesnt do
    // anything
    if (InCell(p, p.cell)) return {p.cell, false};

    // loop over cells attached to the surface
    for (const auto lc : RAJA::range<LocalCellID>(0, surface_num_cells(s))) {
      const auto& c = surface_cells(s, lc);
      if (InCell(p, c)) return {c, false};
    }

    // at a boundary, apply boundary condition
    return ApplyBoundaryCondition(p, s);
  };

 private:
  /**
   * @brief Apply a boundary condition to a particle at specified surface
   */
  RAJA_HOST_DEVICE
  SurfaceCrossingResult ApplyBoundaryCondition(Particle& p,
                                               const SurfaceID& s) const {
    const auto bc_type = static_cast<BoundaryCondition>(
        surfaces(s, SurfaceDataIndex{SURFACE_BC_INDEX}));
    switch (bc_type) {
      case BoundaryCondition::Reflective: {
        const auto stype = static_cast<SurfaceType>(
            surfaces(s, SurfaceDataIndex{SURFACE_TYPE_INDEX}));
        const auto n = Normal(
            stype, p, &surfaces(s, SurfaceDataIndex{SURFACE_DATA_OFFSET}));
        const double dot = n.ux * p.ux + n.uy * p.uy + n.uz * p.uz;
        p.ux -= 2.0 * dot * n.ux;
        p.uy -= 2.0 * dot * n.uy;
        p.uz -= 2.0 * dot * n.uz;
        return {p.cell, false};
      }
      case BoundaryCondition::Vacuum:
        return {CellID{MAXCaerorIndex}, true};
      case BoundaryCondition::Transmission:
        return {CellID{MAXCaerorIndex}, true};
      default:
        return {CellID{MAXCaerorIndex}, true};
    }
  }

  /**
   * @brief Evaluate stack if the particle is perfectly tangent to surface
   */
  RAJA_HOST_DEVICE
  bool InCellTangent(const Particle& p, const CellID c) const {
    SenseResult stack_data[MAX_STACK_DEPTH];
    Stack stack(stack_data);
    for (const auto t : RAJA::range<RPNTokenIndex>(0, logic_sizes(c))) {
      const auto& token = logic(c, t);
      switch (token) {
        case RPN_NOT: {
          const auto value = stack.Pop();
          if (value == SenseResult::Positive)
            stack.Push(SenseResult::Negative);
          else if (value == SenseResult::Negative)
            stack.Push(SenseResult::Positive);
          else
            stack.Push(SenseResult::On);
          break;
        }
        case RPN_AND: {
          const auto rhs = stack.Pop();
          const auto lhs = stack.Pop();
          if (lhs == SenseResult::Negative || rhs == SenseResult::Negative)
            stack.Push(SenseResult::Negative);
          else if (lhs == SenseResult::On || rhs == SenseResult::On)
            stack.Push(SenseResult::On);
          else
            stack.Push(SenseResult::Positive);
          break;
        }
        case RPN_OR: {
          const auto rhs = stack.Pop();
          const auto lhs = stack.Pop();
          if (lhs == SenseResult::Positive || rhs == SenseResult::Positive)
            stack.Push(SenseResult::Positive);
          else if (lhs == SenseResult::On || rhs == SenseResult::On)
            stack.Push(SenseResult::On);
          else
            stack.Push(SenseResult::Negative);
          break;
        }
        default:
          const auto& s = SurfaceID{token};
          const SurfaceDataIndex index{SURFACE_DATA_OFFSET};
          const auto s_type = static_cast<SurfaceType>(
              surfaces(s, SurfaceDataIndex{SURFACE_TYPE_INDEX}));
          const auto& sense = DirectionalSense(s_type, p, &surfaces(s, index));
          stack.Push(sense);
          break;
      }
    }
    const auto result = stack.Pop();
    return result == SenseResult::On ? c == p.cell
                                     : result == SenseResult::Positive;
  };
};
}  // namespace caeror
