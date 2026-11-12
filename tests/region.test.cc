#include "caeror/region.h"

#include <gtest/gtest.h>

#include <vector>

namespace {
using namespace caeror;

TEST(Region, PositiveHalfSpace) {
  const SurfaceBase surface{SurfaceID{1}};
  auto region = +surface;
  const std::vector<RPNToken> expected = {1};
  EXPECT_EQ(region.logic_, expected);
  EXPECT_EQ(region.stack_size, 1);
}

TEST(Region, NegativeHalfSpace) {
  const SurfaceBase surface{SurfaceID{1}};
  auto region = -surface;
  const std::vector<RPNToken> expected = {1, RPN_NOT};
  EXPECT_EQ(region.logic_, expected);
  EXPECT_EQ(region.stack_size, 1);
}

TEST(Region, Intersection) {
  const SurfaceBase a{SurfaceID{0}}, b{SurfaceID{1}};
  auto region = +a && -b;
  const std::vector<RPNToken> expected = {0, 1, RPN_NOT, RPN_AND};
  EXPECT_EQ(region.logic_, expected);
  EXPECT_EQ(region.stack_size, 2);
}

TEST(Region, Union) {
  const SurfaceBase a{SurfaceID{0}}, b{SurfaceID{1}};
  auto region = -a || +b;
  const std::vector<RPNToken> expected = {0, RPN_NOT, 1, RPN_OR};
  EXPECT_EQ(region.logic_, expected);
  EXPECT_EQ(region.stack_size, 2);
}

TEST(Region, Negation) {
  const SurfaceBase a{SurfaceID{0}}, b{SurfaceID{1}};
  auto region = !(+a && -b);
  const std::vector<RPNToken> expected = {0, 1, RPN_NOT, RPN_AND, RPN_NOT};
  EXPECT_EQ(region.logic_, expected);
  EXPECT_EQ(region.stack_size, 2);

  auto twice = !region;
  const std::vector<RPNToken> expected_twice = {0,       1,       RPN_NOT,
                                                RPN_AND, RPN_NOT, RPN_NOT};
  EXPECT_EQ(twice.logic_, expected_twice);
  EXPECT_EQ(twice.stack_size, 2);
}

TEST(Region, NestedExpressions) {
  const SurfaceBase a{SurfaceID{0}}, b{SurfaceID{1}}, c{SurfaceID{2}},
      d{SurfaceID{3}};

  auto left = ((+a && +b) || +c) && +d;
  const std::vector<RPNToken> expected_left = {0,      1, RPN_AND, 2,
                                               RPN_OR, 3, RPN_AND};
  EXPECT_EQ(left.logic_, expected_left);
  EXPECT_EQ(left.stack_size, 2);

  auto right = +a || (+b && (+c || +d));
  const std::vector<RPNToken> expected_right = {0,      1,       2,     3,
                                                RPN_OR, RPN_AND, RPN_OR};
  EXPECT_EQ(right.logic_, expected_right);
  EXPECT_EQ(right.stack_size, 4);

  auto balanced = (+a || +b) && (+c || +d);
  const std::vector<RPNToken> expected_balanced = {0, 1,      RPN_OR, 2,
                                                   3, RPN_OR, RPN_AND};
  EXPECT_EQ(balanced.logic_, expected_balanced);
  EXPECT_EQ(balanced.stack_size, 3);

  auto deeper_left = (+a && (+b && +c)) || +d;
  const std::vector<RPNToken> expected_deeper_left = {
      0, 1, 2, RPN_AND, RPN_AND, 3, RPN_OR};
  EXPECT_EQ(deeper_left.logic_, expected_deeper_left);
  EXPECT_EQ(deeper_left.stack_size, 3);
}
}  // namespace
