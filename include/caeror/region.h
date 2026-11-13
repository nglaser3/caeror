#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>

#include "caeror/raja_layouts.h"
#include "caeror/surface.h"

namespace caeror {
/// @cond
/// @brief Struct for storing a region's logic
struct Region {
  const std::vector<RPNToken> logic_;
  const size_t stack_size;
};
/// @endcond

/**
 * @brief Create a positive halfspace from a surface
 */
inline Region operator+(const SurfaceBase& surf) {
  return Region{{RPNToken{*surf.id_}}, 1};
}

/**
 * @brief Create a negative halfspace from a surface
 */
inline Region operator-(const SurfaceBase& surf) {
  return Region{{RPNToken{*surf.id_}, RPNToken{RPN_NOT}}, 1};
}

/**
 * @brief Create an intersection of two regions
 */
inline Region operator&&(const Region& r1, const Region& r2) {
  std::vector<RPNToken> logic;
  logic.reserve(r1.logic_.size() + r2.logic_.size() + 1);
  logic.insert(logic.end(), r1.logic_.begin(), r1.logic_.end());
  logic.insert(logic.end(), r2.logic_.begin(), r2.logic_.end());
  logic.emplace_back(RPN_AND);
  return Region{logic, std::max(r1.stack_size, 1 + r2.stack_size)};
}

/**
 * @brief Create a union of two regions
 */
inline Region operator||(const Region& r1, const Region& r2) {
  std::vector<RPNToken> logic;
  logic.reserve(r1.logic_.size() + r2.logic_.size() + 1);
  logic.insert(logic.end(), r1.logic_.begin(), r1.logic_.end());
  logic.insert(logic.end(), r2.logic_.begin(), r2.logic_.end());
  logic.emplace_back(RPN_OR);
  return Region{logic, std::max(r1.stack_size, 1 + r2.stack_size)};
}

/**
 * @brief Create the complement of a region
 */
inline Region operator!(const Region& r1) {
  auto logic = r1.logic_;
  logic.emplace_back(RPN_NOT);
  return Region{logic, r1.stack_size};
}
}  // namespace caeror
