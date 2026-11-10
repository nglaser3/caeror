#include <cmath>

#include <gtest/gtest.h>

#include "caeror/surface.h"
#include "caeror/constants.h"

namespace {
using namespace caeror;

// Put the requested position and direction on the selected axis.
Particle MakeParticle(Axis axis, double position, double direction) {
  Particle p{};
  const double tangent = std::sqrt(1.0 - direction * direction);
  switch (axis) {
  case Axis::X:
    p.x = position;
    p.ux = direction;
    p.uy = tangent;
    break;
  case Axis::Y:
    p.y = position;
    p.uy = direction;
    p.uz = tangent;
    break;
  case Axis::Z:
    p.z = position;
    p.uz = direction;
    p.ux = tangent;
    break;
  }
  return p;
}

TEST(AxisAlignedPlane, Sense) {
  for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
    for (double intercept : {-4.0, 0.0, 4.0}) {
      const double data[] = {static_cast<double>(axis), intercept};
      for (double direction : {-1.0, 0.0, 1.0}) {
        auto p = MakeParticle(axis, intercept - 2.0, direction);
        auto result = Sense<SurfaceType::AxisAlignedPlane>(p, data);
        EXPECT_EQ(result, SenseResult::Negative);

        p = MakeParticle(axis, intercept, direction);
        result = Sense<SurfaceType::AxisAlignedPlane>(p, data);
        EXPECT_EQ(result, SenseResult::On);

        p = MakeParticle(axis, intercept + 2.0, direction);
        result = Sense<SurfaceType::AxisAlignedPlane>(p, data);
        EXPECT_EQ(result, SenseResult::Positive);
      }
    }
  }
}

TEST(AxisAlignedPlane, Intersection) {
  struct Case {
    double offset, direction, distance;
    bool valid;
  };
  const Case cases[] = {
      {-3.0, 1.0, 3.0, true},       // Approaching from either side.
      {3.0, -1.0, 3.0, true},
      {-3.0, -1.0, INF, false}, // Moving away.
      {3.0, 1.0, INF, false},
      {-3.0, 0.0, INF, false},  // Parallel, including on the plane.
      {3.0, 0.0, INF, false},
      {0.0, 0.0, INF, false},
      {0.0, -0.0, INF, false},
      {0.0, 1.0, 0.0, true},        // Starting on the plane.
      {0.0, -1.0, 0.0, true},
      {-3.0, 0.6, 5.0, true},       // Oblique unit directions.
      {3.0, -0.6, 5.0, true},
      {-1.0, TOLERANCE, 1.0 / TOLERANCE, true}, // Nearly parallel, approaching.
      {1.0, -TOLERANCE, 1.0 / TOLERANCE, true},
      {-1.0, -TOLERANCE, INF, false},
      {1.0, TOLERANCE, INF, false},
  };
  for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
    for (double intercept : {-4.0, 0.0, 4.0}) {
      const double data[] = {static_cast<double>(axis), intercept};
      for (const auto &c : cases) {
        const auto p = MakeParticle(axis, intercept + c.offset, c.direction);
        const auto result = Intersection<SurfaceType::AxisAlignedPlane>(p, data);
        EXPECT_EQ(result.valid_, c.valid);
        EXPECT_DOUBLE_EQ(result.distance_, c.distance);
      }
    }
  }
}

TEST(AxisAlignedPlane, IgnoresTransverseCoordinates) {
  for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
    const double data[] = {static_cast<double>(axis), 4.0};
    auto p = MakeParticle(axis, 1.0, 1.0);
    if (axis != Axis::X) p.x = 11.0;
    if (axis != Axis::Y) p.y = -17.0;
    if (axis != Axis::Z) p.z = 23.0;

    auto sense = Sense<SurfaceType::AxisAlignedPlane>(p, data);
    EXPECT_EQ(sense, SenseResult::Negative);
    const auto result = Intersection<SurfaceType::AxisAlignedPlane>(p, data);
    EXPECT_TRUE(result.valid_);
    EXPECT_DOUBLE_EQ(result.distance_, 3.0);
    auto aligned = AlignedNormal<SurfaceType::AxisAlignedPlane>(
        p, SenseResult::Positive, data);
    EXPECT_TRUE(aligned);
    aligned = AlignedNormal<SurfaceType::AxisAlignedPlane>(
        p, SenseResult::Negative, data);
    EXPECT_FALSE(aligned);
  }
}

TEST(AxisAlignedPlane, AlignedNormal) {
  for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
    for (double intercept : {-4.0, 0.0, 4.0}) {
      const double data[] = {static_cast<double>(axis), intercept};
      for (double direction : {-1.0, -0.6, -TOLERANCE, -0.0,
                               0.0, TOLERANCE, 0.6, 1.0}) {
        for (double offset : {-2.0, 0.0, 2.0}) {
          const auto p = MakeParticle(axis, intercept + offset, direction);
          auto result = AlignedNormal<SurfaceType::AxisAlignedPlane>(
              p, SenseResult::Positive, data);
          EXPECT_EQ(result, direction > 0.0);

          result = AlignedNormal<SurfaceType::AxisAlignedPlane>(
              p, SenseResult::Negative, data);
          EXPECT_EQ(result, direction < 0.0);

          result = AlignedNormal<SurfaceType::AxisAlignedPlane>(
              p, SenseResult::On, data);
          EXPECT_FALSE(result);
        }
      }
    }
  }
}
} // namespace
