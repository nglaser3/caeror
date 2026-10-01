#pragma once

namespace caeror {
  // ==========================================================================
  // Index definitions
  // ==========================================================================
  using CaerorIndexType = uint64_t;
  using CaerorTokenType = int64_t;
  constexpr CaerorIndexType MAXCaerorIndex = std::numeric_limits<CaerorIndexType>::max();

  RAJA_INDEX_VALUE_T(SurfaceID, CaerorIndexType, "Surface Index");
  RAJA_INDEX_VALUE_T(MaterialID, CaerorIndexType, "Material Index");
  RAJA_INDEX_VALUE_T(CellID, CaerorIndexType, "Cell Index");
  RAJA_INDEX_VALUE_T(UniverseID, CaerorIndexType, "Universe Index");
  RAJA_INDEX_VALUE_T(RPNToken, CaerorIndexType, "Reverse Polish Notation Token");
  
  constexpr CaerorIndexType RPN_NOT = MAXCaerorIndex;
  constexpr CaerorIndexType RPN_AND = MAXCaerorIndex - 1;
  constexpr CaerorIndexType RPN_OR = MAXCaerorIndex - 2;

  // ==========================================================================
  // Layout definitions
  // ==========================================================================
  template <typename ...IndexTypes>
  using CaerorLayout = RAJA::TypedLayout<CaerorIndexType, RAJA::tuple<IndexTypes...>>;

  template <typename DataType, typename ...IndexTypes>
  using CaerorView = RAJA::View<const DataType, CaerorLayout<IndexTypes...>>;

  template <typename DataType, typename ...IndexTypes>
  using CaerorWritableView = RAJA::View<DataType, CaerorLayout<IndexTypes...>>;

} // namespace caeror
