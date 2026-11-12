#include <gtest/gtest.h>

#include <cmath>

#include "caeror/constants.h"
#include "caeror/surface.h"

namespace {
using namespace caeror;

// Map the two radial coordinates and the axial coordinate onto X, Y, or Z.
Particle MakeParticle(Axis axis, double r1, double r2, double axial, double u1,
                      double u2, double ua) {
  Particle p{};
  switch (axis) {
    case Axis::X:
      p.x = axial;
      p.y = r1;
      p.z = r2;
      p.ux = ua;
      p.uy = u1;
      p.uz = u2;
      break;
    case Axis::Y:
      p.x = r1;
      p.y = axial;
      p.z = r2;
      p.ux = u1;
      p.uy = ua;
      p.uz = u2;
      break;
    case Axis::Z:
      p.x = r1;
      p.y = r2;
      p.z = axial;
      p.ux = u1;
      p.uy = u2;
      p.uz = ua;
      break;
  }
  return p;
}

TEST(AxisAlignedCylinder, Sense) {
  for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
    for (double radius : {1.0, 2.0}) {
      const double data[] = {static_cast<double>(axis), radius, 3.0, -4.0};
      for (double axial : {-10.0, 0.0, 10.0}) {
        for (double offset : {-2.0, -1.0, 0.0, 1.0, 2.0}) {
          const auto expected = std::abs(offset) < 1.0 ? SenseResult::Negative
                                : std::abs(offset) == 1.0
                                    ? SenseResult::On
                                    : SenseResult::Positive;
          auto p =
              MakeParticle(axis, 3.0 + offset * radius, -4.0, axial, 1, 0, 0);
          auto result = Sense<SurfaceType::AxisAlignedCylinder>(p, data);
          EXPECT_EQ(result, expected);
          p = MakeParticle(axis, 3.0, -4.0 + offset * radius, axial, 0, 1, 0);
          result = Sense<SurfaceType::AxisAlignedCylinder>(p, data);
          EXPECT_EQ(result, expected);
        }
      }
    }
  }
}

TEST(AxisAlignedCylinder, Intersection) {
  struct Case {
    double r1, r2, u1, u2, ua, distance;
    bool valid;
  };
  const Case cases[] = {
      {-3, 0, 1, 0, 0, 2, true},  // Nearest forward root.
      {3, 0, -1, 0, 0, 2, true},
      {0, -3, 0, 1, 0, 2, true},
      {0, 3, 0, -1, 0, 2, true},
      {0, 0, 1, 0, 0, 1, true},  // Inside: positive root.
      {0.5, 0, -1, 0, 0, 1.5, true},
      {3, 0, 1, 0, 0, INF, false},   // Moving away.
      {-3, 2, 1, 0, 0, INF, false},  // No real roots.
      {-3, 1, 1, 0, 0, 3, true},     // Tangency ahead and behind.
      {3, 1, 1, 0, 0, INF, false},
      {1, 0, 1, 0, 0, 0, true},  // On surface: outward, inward, tangent.
      {1, 0, -1, 0, 0, 0, true},
      {1, 0, 0, 1, 0, 0, true},
      {0, 0, 0.6, 0.8, 0, 1, true},
      {-4, 0, 0.6, 0, 0.8, 5,
       true},  // Axial component changes travel distance.
  };
  for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
    for (double radius : {1.0, 2.0}) {
      for (double shift : {0.0, 4.0}) {
        const double data[] = {static_cast<double>(axis), radius, shift,
                               -shift};
        for (const auto& c : cases) {
          auto p = MakeParticle(axis, shift + radius * c.r1,
                                -shift + radius * c.r2, 7.0, c.u1, c.u2, c.ua);
          auto result =
              Intersection<SurfaceType::AxisAlignedCylinder, false>(p, data);
          EXPECT_EQ(result.valid_, c.valid);
          EXPECT_DOUBLE_EQ(result.distance_, radius * c.distance);
        }
      }
    }
  }
}

TEST(AxisAlignedCylinder, SkipZeroFindsFarRoot) {
  const double data[] = {static_cast<double>(Axis::Z), 2.0, 0.0, 0.0};
  Particle p{};
  p.x = -2.0;
  p.ux = 1.0;
  auto result = Intersection<SurfaceType::AxisAlignedCylinder, true>(p, data);
  EXPECT_TRUE(result.valid_);
  EXPECT_DOUBLE_EQ(result.distance_, 4.0);
}

TEST(AxisAlignedCylinder, ParallelIntersection) {
  for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
    const double data[] = {static_cast<double>(axis), 1.0, 0.0, 0.0};
    for (double radius : {0.0, 1.0, 2.0}) {
      for (double direction : {-1.0, 1.0}) {
        auto p = MakeParticle(axis, radius, 0, 0, 0, 0, direction);
        auto result =
            Intersection<SurfaceType::AxisAlignedCylinder, false>(p, data);
        EXPECT_FALSE(result.valid_);
        // Distance is not used when valid_ is false.
      }
    }
  }
}

TEST(AxisAlignedCylinder, AlignedNormal) {
  for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
    const double data[] = {static_cast<double>(axis), 2.0, 3.0, -4.0};
    for (double side : {-1.0, 1.0}) {
      for (double direction : {-1.0, -0.6, 0.0, 0.6, 1.0}) {
        for (bool second_coordinate : {false, true}) {
          const double tangent = std::sqrt(1.0 - direction * direction);
          auto p = second_coordinate ? MakeParticle(axis, 3, -4 + 2 * side, 7,
                                                    0, direction, tangent)
                                     : MakeParticle(axis, 3 + 2 * side, -4, 7,
                                                    direction, 0, tangent);
          auto result = AlignedNormal<SurfaceType::AxisAlignedCylinder>(
              p, SenseResult::Positive, data);
          EXPECT_EQ(result, side * direction > 0.0);
          result = AlignedNormal<SurfaceType::AxisAlignedCylinder>(
              p, SenseResult::Negative, data);
          EXPECT_EQ(result, side * direction < 0.0);
          result = AlignedNormal<SurfaceType::AxisAlignedCylinder>(
              p, SenseResult::On, data);
          EXPECT_FALSE(result);
        }
      }
    }
  }
}
TEST(AxisAlignedCylinder, Normal) {
  for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
    const double data[] = {static_cast<double>(axis), 5.0, 1.0, -2.0};
    auto p = MakeParticle(axis, 4.0, -6.0, 7.0, 1.0, 0.0, 0.0);
    auto raw = Normal<SurfaceType::AxisAlignedCylinder>(p, data);
    EXPECT_DOUBLE_EQ(raw.ux, axis == Axis::X ? 0.0 : 3.0);
    EXPECT_DOUBLE_EQ(raw.uy, axis == Axis::X   ? 3.0
                             : axis == Axis::Y ? 0.0
                                               : -4.0);
    EXPECT_DOUBLE_EQ(raw.uz, axis == Axis::Z ? 0.0 : -4.0);
    auto result = Normal(SurfaceType::AxisAlignedCylinder, p, data);
    EXPECT_DOUBLE_EQ(result.ux, raw.ux / 5.0);
    EXPECT_DOUBLE_EQ(result.uy, raw.uy / 5.0);
    EXPECT_DOUBLE_EQ(result.uz, raw.uz / 5.0);
  }
}

}  // namespace
