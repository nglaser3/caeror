#include <cmath>

#include <gtest/gtest.h>

#include "caeror/constants.h"
#include "caeror/surface.h"

namespace {
using namespace caeror;

TEST(Sphere, Sense) {
  for (double radius : {1.0, 2.0}) {
    for (double shift : {0.0, 4.0}) {
      const double data[] = {radius, shift, -shift, 2.0 * shift};
      for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
        for (double offset : {-2.0, -1.0, 0.0, 1.0, 2.0}) {
          Particle p{};
          p.x = shift; p.y = -shift; p.z = 2.0 * shift;
          p.ux = 1.0;
          if (axis == Axis::X) p.x += offset * radius;
          if (axis == Axis::Y) p.y += offset * radius;
          if (axis == Axis::Z) p.z += offset * radius;
          const auto expected = std::abs(offset) < 1.0 ? SenseResult::Negative :
              std::abs(offset) == 1.0 ? SenseResult::On : SenseResult::Positive;
          auto result = Sense<SurfaceType::Sphere>(p, data);
          EXPECT_EQ(result, expected);
        }
      }
    }
  }
}

TEST(Sphere, Intersection) {
  struct Case {
    double x, y, z, ux, uy, uz, distance;
    bool valid;
  };
  // Positions and distances below are for a unit sphere at the origin.
  const Case cases[] = {
      {-3, 0, 0, 1, 0, 0, 2, true}, // Nearest of two forward roots.
      {0, -3, 0, 0, 1, 0, 2, true},
      {0, 0, -3, 0, 0, 1, 2, true},
      {3, 0, 0, -1, 0, 0, 2, true},
      {0, 0, 0, 1, 0, 0, 1, true},  // Inside: use the positive root.
      {0.5, 0, 0, -1, 0, 0, 1.5, true},
      {3, 0, 0, 1, 0, 0, INF, false}, // Both roots behind.
      {-3, 2, 0, 1, 0, 0, INF, false}, // No real roots.
      {-3, 1, 0, 1, 0, 0, 3, true}, // Tangency ahead and behind.
      {3, 1, 0, 1, 0, 0, INF, false},
      {1, 0, 0, 1, 0, 0, 0, true},  // On surface: outward, inward, tangent.
      {1, 0, 0, -1, 0, 0, 0, true},
      {1, 0, 0, 0, 1, 0, 0, true},
      {0, 0, 0, 0.6, 0.8, 0, 1, true},
  };
  for (double radius : {1.0, 2.0}) {
    for (double shift : {0.0, 4.0}) {
      const double data[] = {radius, shift, -shift, 2.0 * shift};
      for (const auto &c : cases) {
        Particle p{};
        p.x = shift + radius * c.x;
        p.y = -shift + radius * c.y;
        p.z = 2.0 * shift + radius * c.z;
        p.ux = c.ux; p.uy = c.uy; p.uz = c.uz;
        auto result = Intersection<SurfaceType::Sphere, false>(p, data);
        EXPECT_EQ(result.valid_, c.valid);
        EXPECT_DOUBLE_EQ(result.distance_, radius * c.distance);
      }
    }
  }
}

TEST(Sphere, SkipZeroFindsFarRoot) {
  const double data[] = {2.0, 0.0, 0.0, 0.0};
  Particle p{};
  p.x = -2.0;
  p.ux = 1.0;
  auto result = Intersection<SurfaceType::Sphere, true>(p, data);
  EXPECT_TRUE(result.valid_);
  EXPECT_DOUBLE_EQ(result.distance_, 4.0);
}

TEST(Sphere, ExitReportsZeroDistance) {
  const double data[] = {2.0, 0.0, 0.0, 0.0};
  Particle p{};
  p.x = 2.0;
  p.ux = 1.0;
  auto result = Intersection<SurfaceType::Sphere, false>(p, data);
  EXPECT_TRUE(result.valid_);
  EXPECT_DOUBLE_EQ(result.distance_, 0.0);
}

TEST(Sphere, AlignedNormal) {
  for (double shift : {0.0, 4.0}) {
    const double data[] = {2.0, shift, -shift, 2.0 * shift};
    for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
      for (double side : {-1.0, 1.0}) {
        for (double direction : {-1.0, -0.6, 0.0, 0.6, 1.0}) {
          Particle p{};
          p.x = shift; p.y = -shift; p.z = 2.0 * shift;
          const double tangent = std::sqrt(1.0 - direction * direction);
          if (axis == Axis::X) { p.x += 2.0 * side; p.ux = direction; p.uy = tangent; }
          if (axis == Axis::Y) { p.y += 2.0 * side; p.uy = direction; p.uz = tangent; }
          if (axis == Axis::Z) { p.z += 2.0 * side; p.uz = direction; p.ux = tangent; }
          auto result = AlignedNormal<SurfaceType::Sphere>(p, SenseResult::Positive, data);
          EXPECT_EQ(result, side * direction > 0.0);
          result = AlignedNormal<SurfaceType::Sphere>(p, SenseResult::Negative, data);
          EXPECT_EQ(result, side * direction < 0.0);
          result = AlignedNormal<SurfaceType::Sphere>(p, SenseResult::On, data);
          EXPECT_FALSE(result);
        }
      }
    }
  }
}
TEST(Sphere, Normal) {
  const double data[] = {5.0, 1.0, 2.0, 7.0};
  Particle p{};
  p.x = 4.0;
  p.y = -2.0;
  p.z = 7.0;
  auto raw = Normal<SurfaceType::Sphere>(p, data);
  EXPECT_DOUBLE_EQ(raw.ux, 3.0);
  EXPECT_DOUBLE_EQ(raw.uy, -4.0);
  EXPECT_DOUBLE_EQ(raw.uz, 0.0);
  auto result = Normal(SurfaceType::Sphere, p, data);
  EXPECT_DOUBLE_EQ(result.ux, 0.6);
  EXPECT_DOUBLE_EQ(result.uy, -0.8);
  EXPECT_DOUBLE_EQ(result.uz, 0.0);
}

} // namespace
