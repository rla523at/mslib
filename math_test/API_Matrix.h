#pragma once
#include "msmath/Matrix.h"
#include "gtest/gtest.h"

namespace matrix_test_API
{
inline void compare_considering_4ULP(const ms::math::Matrix_Const_Wrapper& cmw1, const ms::math::Matrix_Const_Wrapper& cmw2)
{
  const auto [num_rows, num_columns] = cmw1.size();

  for (int i = 0; i < num_rows; i++)
  {
    for (int j = 0; j < num_columns; j++)
    {
      EXPECT_DOUBLE_EQ(cmw1.at(i, j), cmw2.at(i, j));
    }
  }
}

inline void compare_considering_epsilon(const ms::math::Matrix_Const_Wrapper& cmw1, const ms::math::Matrix_Const_Wrapper& cmw2, const double epsilon)
{
  const auto [num_rows, num_columns] = cmw1.size();

  for (int i = 0; i < num_rows; i++)
  {
    for (int j = 0; j < num_columns; j++)
    {
      EXPECT_NEAR(cmw1.at(i, j), cmw2.at(i, j), epsilon);
    }
  }
}
} // namespace matrix_test_API