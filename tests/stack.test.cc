#include "caeror/stack.h"

#include <gtest/gtest.h>

namespace {
using namespace caeror;

TEST(Stack, Empty) {
  int data[4]{};
  Stack<int> stack(data);
  EXPECT_TRUE(stack.Empty());

  stack.Push(1);
  EXPECT_FALSE(stack.Empty());

  EXPECT_EQ(stack.Pop(), 1);
  EXPECT_TRUE(stack.Empty());
}

TEST(Stack, PushPop) {
  int data[4]{};
  Stack<int> stack(data);
  stack.Push(1);
  stack.Push(2);
  stack.Push(3);

  EXPECT_EQ(stack.Pop(), 3);
  stack.Push(4);
  EXPECT_EQ(stack.Pop(), 4);
  EXPECT_EQ(stack.Pop(), 2);
  EXPECT_EQ(stack.Pop(), 1);

  EXPECT_TRUE(stack.Empty());
}

TEST(Stack, Clear) {
  int data[4]{};
  Stack<int> stack(data);
  stack.Clear();
  EXPECT_TRUE(stack.Empty());

  stack.Push(20);
  stack.Clear();
  EXPECT_TRUE(stack.Empty());
}

TEST(Stack, BooleanValues) {
  bool data[4]{};
  Stack<bool> stack(data);
  stack.Push(true);
  stack.Push(false);
  stack.Push(true);

  EXPECT_TRUE(stack.Pop());
  EXPECT_FALSE(stack.Pop());
  EXPECT_TRUE(stack.Pop());

  EXPECT_TRUE(stack.Empty());
}
}  // namespace
