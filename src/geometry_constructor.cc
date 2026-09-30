#include "caeror/geometry_constructor.h"

namespace caeror
{
  Universe GeometryConstructor::CreateUniverse(std::vector<Cell> cells) {
    auto id = universes_.size();
    std::vector<CellID> cell_ids;
    for (const auto& cell : cells) cell_ids.push_back(cell.id_);
    auto uni_ptr = std::make_unique<Universe>(UniverseID{id}, cell_ids);
    universes_.push_back(std::move(uni_ptr));
    return *universes_.back();
  }
  
} // namespace caeror
