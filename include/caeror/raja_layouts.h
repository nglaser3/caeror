#pragma once

namespace caeror {
  // ==========================================================================
  // Index definitions
  // ==========================================================================
  using CaerorIndexType = uint64_t;

  RAJA_INDEX_VALUE_T(SurfaceID, CaerorIndexType, "Surface Index");
  RAJA_INDEX_VALUE_T(MaterialID, CaerorIndexType, "Material Index");
  RAJA_INDEX_VALUE_T(CellID, CaerorIndexType, "Cell Index");
  RAJA_INDEX_VALUE_T(UniverseID, CaerorIndexType, "Universe Index");

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
