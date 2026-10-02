#pragma once

namespace caeror {
  // ==========================================================================
  // Index definitions
  // ==========================================================================
  using CaerorIndexType = uint64_t;
  using RPNToken = uint64_t;
  constexpr CaerorIndexType MAXCaerorIndex = std::numeric_limits<CaerorIndexType>::max();

  RAJA_INDEX_VALUE_T(SurfaceID, CaerorIndexType, "Surface Index");
  RAJA_INDEX_VALUE_T(SurfaceDataIndex, CaerorIndexType, "Surface Data Index");
  RAJA_INDEX_VALUE_T(MaterialID, CaerorIndexType, "Material Index");
  RAJA_INDEX_VALUE_T(CellID, CaerorIndexType, "Cell Index");
  RAJA_INDEX_VALUE_T(RPNTokenIndex, CaerorIndexType, "Reverse Polish Notation Token");
  RAJA_INDEX_VALUE_T(UniverseID, CaerorIndexType, "Universe Index");
  RAJA_INDEX_VALUE_T(LocalCellID, CaerorIndexType, "Local Cell Index");

  constexpr RPNToken RPN_NOT = MAXCaerorIndex;
  constexpr RPNToken RPN_AND = MAXCaerorIndex - 1;
  constexpr RPNToken RPN_OR = MAXCaerorIndex - 2;

  // ==========================================================================
  // Layout definitions
  // ==========================================================================
  template <typename ...IndexTypes>
  using CaerorLayout = RAJA::TypedLayout<CaerorIndexType, RAJA::tuple<IndexTypes...>>;

  using SurfaceDataLayout = CaerorLayout<SurfaceID, SurfaceDataIndex>;
  using SurfaceCellLayout = CaerorLayout<SurfaceID, LocalCellID>;
  using SurfaceLayout = CaerorLayout<SurfaceID>;
  using CellLayout = CaerorLayout<CellID>;
  using RPNTokenLayout = CaerorLayout<CellID, RPNTokenIndex>;
  using UniverseLayout = CaerorLayout<UniverseID>;
  using UniverseCellLayout = CaerorLayout<UniverseID, LocalCellID>;

  template <typename DataType, typename ...IndexTypes>
  using CaerorView = RAJA::View<const DataType, CaerorLayout<IndexTypes...>>;

  template <typename DataType, typename ...IndexTypes>
  using CaerorWritableView = RAJA::View<DataType, CaerorLayout<IndexTypes...>>;

} // namespace caeror
