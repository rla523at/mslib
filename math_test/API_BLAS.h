#pragma once
#include "msmath/BLAS.h"
#include "gtest/gtest.h"

namespace blas_test_API
{
inline void compare_considering_4ULP(const std::vector<double>& vec1, const std::vector<double>& vec2)
{
  for (size_t i = 0; i < vec1.size(); ++i)
  {
    EXPECT_DOUBLE_EQ(vec1[i], vec2[i]);
  }
}
} // namespace blas_test_API