#pragma once
#include "API_Matrix.h"

namespace ms::math
{

TEST(Matrix, cwrapper1)
{
  const Matrix m       = {1, 2, {1, 2}};
  auto         m_cwrap = m.view();

  EXPECT_EQ(m.at(0, 1), 2);
}
TEST(Matrix, conversion2)
{
  Matrix m        = {1, 2, {1, 2}};
  auto   m_wrap   = m.wrap();
  m_wrap.at(0, 1) = 5;

  EXPECT_EQ(m.at(0, 1), 5);
}

TEST(Matrix, operator_addition_assign_1)
{
  constexpr auto      num_rows    = 1;
  constexpr auto      num_columns = 2;
  std::vector<double> val1        = {1, 2};
  std::vector<double> val2        = {2, 2};

  Matrix            m   = {num_rows, num_columns, std::move(val1)};
  const Matrix_View cmw = {num_rows, num_columns, val2};
  m += cmw;

  const Matrix ref = {1, 2, {3, 4}};
  EXPECT_EQ(m, ref);
}
TEST(Matrix, operator_addition_assign_2)
{
  constexpr auto      num_rows    = 2;
  constexpr auto      num_columns = 1;
  std::vector<double> val1        = {1, 2};
  std::vector<double> val2        = {2, 2};

  Matrix            m   = {num_rows, num_columns, std::move(val1)};
  const Matrix_View cmw = {num_rows, num_columns, val2};
  m += cmw;

  const Matrix ref = {2, 1, {3, 4}};
  EXPECT_EQ(m, ref);
}
TEST(Matrix, operator_addition_assign_3)
{
  constexpr auto      num_rows    = 2;
  constexpr auto      num_columns = 2;
  std::vector<double> val1        = {1, 2, 3, 4};
  std::vector<double> val2        = {1, 3, 4, 5};

  Matrix            m   = {num_rows, num_columns, std::move(val1)};
  const Matrix_View cmw = {num_rows, num_columns, val2};
  m += cmw;

  Matrix ref(2, 2, {2, 5, 7, 9});
  EXPECT_EQ(m, ref);
}
TEST(Matrix, operator_addition_assign_6)
{
  constexpr auto      num_rows1    = 2;
  constexpr auto      num_columns1 = 5;
  std::vector<double> val1         = {1.2345, 2.346345, 6.3262345, 8.5674567, 6.23452345, 2.53462, 6.432452345, 2.345345, 1.234563245, 7.3245345};

  constexpr auto      num_rows2    = 5;
  constexpr auto      num_columns2 = 2;
  std::vector<double> val2         = {1.234234, 2.3462345, 345.324, 2.6345345, 634523.5, 2345345.3, 23453.345, 234534.6, 234523.5, 623452.1};

  Matrix            m    = {num_rows1, num_columns1, std::move(val1)};
  const Matrix_View cmw  = {num_rows2, num_columns2, val2};
  const auto        cmwT = cmw.transpose();
  m += cmwT;

  Matrix ref(2, 5, {2.468734, 347.670345, 634529.8262345, 23461.9124567, 234529.73452345, 4.8808545, 9.066986845, 2345347.645345, 234535.834563245, 623459.4245345});
  matrix_test_API::compare_considering_4ULP(m, ref);
}
TEST(Matrix, operator_addition_assign_7)
{
  constexpr auto      num_rows1    = 3;
  constexpr auto      num_columns1 = 5;
  std::vector<double> val1         = {1.2345, 2.346345, 6.3262345, 8.5674567, 6.23452345, 2.53462, 6.432452345, 2.345345, 1.234563245, 7.3245345, 789.45978, 74.5789123, 74.23541, 4.7894113, 7894.5134};

  constexpr auto      num_rows2    = 5;
  constexpr auto      num_columns2 = 3;
  std::vector<double> val2         = {1.234234, 2.3462345, 789456.0, 345.324, 2.6345345, 74.48651, 634523.5, 2345345.3, 710.1846, 23453.345, 234534.6, 12.5487, 234523.5, 623452.1, 421.7456};

  Matrix            m    = {num_rows1, num_columns1, std::move(val1)};
  const Matrix_View cmw  = {num_rows2, num_columns2, val2};
  const auto        cmwT = cmw.transpose();
  m += cmwT;

  Matrix ref(3, 5, {2.468734, 347.670345, 634529.8262345, 23461.9124567, 234529.73452345, 4.8808545, 9.066986845, 2345347.645345, 234535.834563245, 623459.4245345, 790245.45978, 149.0654223, 784.42001, 17.3381113, 8316.259});
  matrix_test_API::compare_considering_4ULP(m, ref);
}

TEST(blas, mpm)
{
  constexpr auto      num_rows1    = 2;
  constexpr auto      num_columns1 = 2;
  std::vector<double> val1         = {1, 2, 3, 4};

  constexpr auto      num_rows2    = 2;
  constexpr auto      num_columns2 = 2;
  std::vector<double> val2         = {1, 2, 3, 4};

  constexpr auto      num_rows3          = 2;
  constexpr auto      num_columns3       = 2;
  constexpr auto      leading_dimension3 = 3;
  std::vector<double> val3(num_rows3 * leading_dimension3);

  const auto cmw1   = Matrix_View{num_rows1, num_columns1, val1};
  const auto cmw2   = Matrix_View{num_rows2, num_columns2, val2};
  auto       result = Matrix_Wrap(num_rows3, num_columns3, val3, leading_dimension3);

  ms::math::blas::mpm(result, cmw1, cmw2);

  std::vector<double> ref = {2, 4, 0, 6, 8, 0};
  EXPECT_EQ(val3, ref);
}

} // namespace ms::math
