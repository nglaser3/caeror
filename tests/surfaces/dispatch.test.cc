#include <gtest/gtest.h>

#include "caeror/surface.h"

namespace {
using namespace caeror;

template <SurfaceType Type>
void CheckDispatch(const double *data) {
  for (double position : {-3.0, 0.0, 1.0, 3.0}) {
    for (double direction : {-1.0, 0.0, 1.0}) {
      Particle p{};
      p.x = position;
      p.ux = direction;
      p.uy = direction == 0.0 ? 1.0 : 0.0;

      auto sense = Sense(Type, p, data);
      auto expected_sense = Sense<Type>(p, data);
      EXPECT_EQ(sense, expected_sense);

      auto hit = Intersection<false>(Type, p, data);
      auto expected_hit = Intersection<Type, false>(p, data);
      EXPECT_EQ(hit.valid_, expected_hit.valid_);
      EXPECT_DOUBLE_EQ(hit.distance_, expected_hit.distance_);

      for (auto side : {SenseResult::Negative, SenseResult::On,
                        SenseResult::Positive}) {
        auto aligned = AlignedNormal(Type, p, side, data);
        auto expected_aligned = AlignedNormal<Type>(p, side, data);
        EXPECT_EQ(aligned, expected_aligned);
      }
    }
  }
}

TEST(SurfaceDispatch, Plane) {
  const double data[] = {1.0, 2.0, 3.0, -2.0};
  CheckDispatch<SurfaceType::Plane>(data);
}

TEST(SurfaceDispatch, AxisAlignedPlane) {
  const double data[] = {static_cast<double>(Axis::X), 1.0};
  CheckDispatch<SurfaceType::AxisAlignedPlane>(data);
}

TEST(SurfaceDispatch, Sphere) {
  const double data[] = {2.0, 1.0, -1.0, 0.0};
  CheckDispatch<SurfaceType::Sphere>(data);
}

TEST(SurfaceDispatch, AxisAlignedCylinder) {
  const double data[] = {static_cast<double>(Axis::Z), 2.0, 1.0, -1.0};
  CheckDispatch<SurfaceType::AxisAlignedCylinder>(data);
}

TEST(SurfaceDispatch, DirectionalSense) {
  const double data[] = {static_cast<double>(Axis::X), 2.0};
  Particle p{};

  p.x = 1.0;
  p.ux = 1.0;
  auto result = DirectionalSense(SurfaceType::AxisAlignedPlane, p, data);
  EXPECT_EQ(result, SenseResult::Negative);

  p.x = 3.0;
  p.ux = -1.0;
  result = DirectionalSense(SurfaceType::AxisAlignedPlane, p, data);
  EXPECT_EQ(result, SenseResult::Positive);

  p.x = 2.0;
  p.ux = 1.0;
  result = DirectionalSense(SurfaceType::AxisAlignedPlane, p, data);
  EXPECT_EQ(result, SenseResult::Positive);

  p.ux = -1.0;
  result = DirectionalSense(SurfaceType::AxisAlignedPlane, p, data);
  EXPECT_EQ(result, SenseResult::Negative);

  p.ux = 0.0;
  p.uy = 1.0;
  result = DirectionalSense(SurfaceType::AxisAlignedPlane, p, data);
  EXPECT_EQ(result, SenseResult::On);
}
} // namespace
