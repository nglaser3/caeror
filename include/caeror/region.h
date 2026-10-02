#pragma once

#include <vector>

#include "caeror/raja_layouts.h"
#include "caeror/surface.h"

namespace caeror {
struct Region {
  const std::vector<RPNToken> logic_;
  const size_t stack_size;
};

Region operator+(const SurfaceBase &surf) {
  return Region{{RPNToken{*surf.id_}}, 1};
}

Region operator-(const SurfaceBase &surf) {
  return Region{{RPNToken{*surf.id_}, RPNToken{RPN_NOT}}, 1};
}

Region operator&&(const Region &r1, const Region &r2) {
  std::vector<RPNToken> logic;
  logic.reserve(r1.logic_.size() + r2.logic_.size() + 1);
  logic.insert(logic.end(), r1.logic_.begin(), r1.logic_.end());
  logic.insert(logic.end(), r2.logic_.begin(), r2.logic_.end());
  logic.emplace_back(RPN_AND);
  return Region{logic, std::max(r1.stack_size, 1 + r2.stack_size)};
}

Region operator||(const Region &r1, const Region &r2) {
  std::vector<RPNToken> logic;
  logic.reserve(r1.logic_.size() + r2.logic_.size() + 1);
  logic.insert(logic.end(), r1.logic_.begin(), r1.logic_.end());
  logic.insert(logic.end(), r2.logic_.begin(), r2.logic_.end());
  logic.emplace_back(RPN_OR);
  return Region{logic, std::max(r1.stack_size, 1 + r2.stack_size)};
}

Region operator!(const Region &r1) {
  auto logic = r1.logic_;
  logic.emplace_back(RPN_NOT);
  return Region{logic, r1.stack_size};
}
} // namespace caeror
