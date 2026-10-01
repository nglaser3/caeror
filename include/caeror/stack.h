#pragma once

#include <RAJA/RAJA.hpp>

#include "caeror/raja_layouts.h"

namespace caeror {
  
  template <typename DataType>
  class Stack{
    public:
      RAJA_HOST_DEVICE
      Stack(DataType* data)
        : eval_stack_(data),
          size_(){}

      RAJA_HOST_DEVICE
      void Push(DataType value) {
        eval_stack_[size_++] = value;
      }

      RAJA_HOST_DEVICE
      DataType Pop() {
        return eval_stack_[--size_];
      }

      RAJA_HOST_DEVICE
      void Clear() {
        size_ = 0;
      }

      RAJA_HOST_DEVICE
      bool Empty() const {return size_ == 0;}

    private:
      DataType* eval_stack_;
      CaerorIndexType size_;
  };
} // namespace caeror
