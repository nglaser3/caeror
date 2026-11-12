#pragma once

#include <RAJA/RAJA.hpp>
#include <cstddef>

#include "caeror/raja_layouts.h"

namespace caeror {

constexpr size_t MAX_STACK_DEPTH = 16;

/**
 * @brief Helper class for evaluating a stack
 */
template <typename DataType>
class Stack {
 public:
  /**
   * @brief create a stack instance from provided data
   */
  RAJA_HOST_DEVICE
  Stack(DataType* data) : eval_stack_(data), size_() {}

  /**
   * @brief Push value to top of stack
   */
  RAJA_HOST_DEVICE
  void Push(DataType value) { eval_stack_[size_++] = value; }

  /**
   * @brief Pop value from top of stack
   */
  RAJA_HOST_DEVICE
  DataType Pop() { return eval_stack_[--size_]; }

  /**
   * @brief Clear the stack
   */
  RAJA_HOST_DEVICE
  void Clear() { size_ = 0; }

  /**
   * @brief Determine if stack is empty
   */
  RAJA_HOST_DEVICE
  bool Empty() const { return size_ == 0; }

 private:
  /// @brief Array of data
  DataType* eval_stack_;
  /// @brief Current size of stack
  CaerorIndexType size_;
};
}  // namespace caeror
