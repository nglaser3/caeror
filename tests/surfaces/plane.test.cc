#include <cmath>

#include <gtest/gtest.h>

#include "caeror/constants.h"
#include "caeror/surface.h"

namespace {
using namespace caeror;

TEST(Plane, Sense) {
  // All coefficients differ so swapping coordinates changes the answer.
  for (double scale : {1.0, 2.0, -1.0}) {
    const double data[] = {scale, 2.0 * scale, 3.0 * scale, -6.0 * scale};
    Particle p{};
    p.ux = 1.0;
    auto result = Sense<SurfaceType::Plane>(p, data);
    EXPECT_EQ(result, scale > 0.0 ? SenseResult::Negative : SenseResult::Positive);

    p.x = 1.0;
    p.y = 1.0;
    p.z = 1.0;
    result = Sense<SurfaceType::Plane>(p, data);
    EXPECT_EQ(result, SenseResult::On);

    p.z = 2.0;
    result = Sense<SurfaceType::Plane>(p, data);
    EXPECT_EQ(result, scale > 0.0 ? SenseResult::Positive : SenseResult::Negative);
  }
}

TEST(Plane, Intersection) {
  struct Case {
    double x, y, z, ux, uy, uz, distance;
    bool valid;
  };
  // Plane x + 2y + 3z = 6; every direction is a unit vector.
  const Case cases[] = {
      {0, 0, 0, 1, 0, 0, 6, true},
      {0, 0, 0, 0, 1, 0, 3, true},
      {0, 0, 0, 0, 0, 1, 2, true},
      {8, 0, 0, -1, 0, 0, 2, true},
      {0, 0, 0, -1, 0, 0, INF, false},
      {8, 0, 0, 1, 0, 0, INF, false},
      {1, 1, 1, 1, 0, 0, 0, true},
      {1, 1, 1, -1, 0, 0, 0, true},
      {0, 0, 0, 0, 0.6, 0.8, 5.0 / 3.0, true},
  };
  for (double scale : {1.0, 2.0, -1.0}) {
    const double data[] = {scale, 2.0 * scale, 3.0 * scale, -6.0 * scale};
    for (const auto &c : cases) {
      Particle p{};
      p.x = c.x; p.y = c.y; p.z = c.z;
      p.ux = c.ux; p.uy = c.uy; p.uz = c.uz;
      auto result = Intersection<SurfaceType::Plane, false>(p, data);
      EXPECT_EQ(result.valid_, c.valid);
      EXPECT_DOUBLE_EQ(result.distance_, c.distance);
    }
  }
}

TEST(Plane, ParallelIntersection) {
  // This unit direction is exactly tangent to x + y = 2.
  const double data[] = {1.0, 1.0, 0.0, -2.0};
  for (double x : {0.0, 2.0, 4.0}) {
    Particle p{};
    p.x = x;
    p.ux = std::sqrt(0.5);
    p.uy = -p.ux;
    auto result = Intersection<SurfaceType::Plane, false>(p, data);
    EXPECT_FALSE(result.valid_);
    EXPECT_EQ(result.distance_, INF);
  }
}

TEST(Plane, AlignedNormal) {
  struct Case {
    double ux, uy, uz;
    bool positive, negative;
  };
  const Case cases[] = {
      {1, 0, 0, true, false}, {0, 1, 0, true, false},
      {0, 0, 1, false, true}, {-1, 0, 0, false, true},
      {0, -1, 0, false, true}, {0, 0, -1, true, false},
      {0, 0.6, 0.8, false, true},
      {0, 0.8, 0.6, true, false},
      {std::sqrt(0.5), -std::sqrt(0.5), 0, false, false},
  };
  for (double d : {-4.0, 0.0, 4.0}) {
    const double data[] = {1.0, 1.0, -1.0, d};
    for (const auto &c : cases) {
      Particle p{};
      p.ux = c.ux; p.uy = c.uy; p.uz = c.uz;
      auto result = AlignedNormal<SurfaceType::Plane>(p, SenseResult::Positive, data);
      EXPECT_EQ(result, c.positive);
      result = AlignedNormal<SurfaceType::Plane>(p, SenseResult::Negative, data);
      EXPECT_EQ(result, c.negative);
      result = AlignedNormal<SurfaceType::Plane>(p, SenseResult::On, data);
      EXPECT_FALSE(result);
    }
  }
}
} // namespace
