#pragma once

#include <string>

#include "caeror/surface.h"

namespace caeror
{
  struct Region {
    const std::string logic_;
  };

  Region operator+(const SurfaceBase& surf) {
    return Region{"+" + std::to_string(*surf.id_)};
  }

  Region operator-(const SurfaceBase& surf) {
    return Region{"-" + std::to_string(*surf.id_)};
  }

  Region operator&&(const Region& r1, const Region& r2) {
    std::string logic = "(" + r1.logic_ + ")&(" + r2.logic_ + ")";
    return Region{logic};
  }

  Region operator||(const Region& r1, const Region& r2) {
    std::string logic = "(" + r1.logic_ + ")|(" + r2.logic_ + ")";
    return Region{logic};
  }

  Region operator!(const Region& r1) {
    std::string logic = "!(" + r1.logic_ + ")";
    return Region{logic};
  }
} // namespace caeror
