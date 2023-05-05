#pragma once
#include "msmath/Vector.h"
#include "API_Matrix.h"

namespace ms::math
{
// for google test cout message

template <size_t dim>
std::ostream& operator<<(std::ostream& os, const Vector<dim>& vec)
{
  return os << vec.to_string();
}

// for google test cout message

TEST(Matrix_Const_Wrapper, operator_equal_1)
{
  constexpr auto      num_row           = 3;
  constexpr auto      num_column        = 2;
  std::vector<double> val1              = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  constexpr auto      leading_dimension = 3;

  const Matrix_Const_Wrapper m = {num_row, num_column, val1.data(), leading_dimension, Transpose_Type::transposed};
  //		m =	¦£ 1 2 ¦¤T	 =	¦¢	1 4 7	¦¢
  //				¦¢ 4 5 ¦¢				¦¢	2 5 8	¦¢
  //				¦¦ 7 8 ¦¥

  constexpr auto      ref_num_row    = 2;
  constexpr auto      ref_num_column = 3;
  std::vector<double> ref_val        = {1, 4, 7, 2, 5, 8};

  const Matrix_Const_Wrapper ref = {ref_num_row, ref_num_column, ref_val.data()};

  EXPECT_TRUE(m == ref);
}
TEST(Matrix_Const_Wrapper, operator_addition_1)
{
  constexpr auto      num_rows    = 1;
  constexpr auto      num_columns = 2;
  std::vector<double> val1        = {1, 2};
  std::vector<double> val2        = {2, 2};

  const Matrix_Const_Wrapper cmw1   = {num_rows, num_columns, val1.data()};
  const Matrix_Const_Wrapper cmw2   = {num_rows, num_columns, val2.data()};
  const auto                 result = cmw1 + cmw2;

  const Matrix ref = {1, 2, {3, 4}};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_addition_2)
{
  constexpr auto      num_rows    = 2;
  constexpr auto      num_columns = 1;
  std::vector<double> val1        = {1, 2};
  std::vector<double> val2        = {2, 2};

  const Matrix_Const_Wrapper cmw1   = {num_rows, num_columns, val1.data()};
  const Matrix_Const_Wrapper cmw2   = {num_rows, num_columns, val2.data()};
  const auto                 result = cmw1 + cmw2;

  const Matrix ref = {2, 1, {3, 4}};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_addition_3)
{
  constexpr auto      num_rows    = 2;
  constexpr auto      num_columns = 2;
  std::vector<double> val1        = {1, 2, 3, 4};
  std::vector<double> val2        = {1, 3, 4, 5};

  const Matrix_Const_Wrapper cmw1   = {num_rows, num_columns, val1.data()};
  const Matrix_Const_Wrapper cmw2   = {num_rows, num_columns, val2.data()};
  const auto                 result = cmw1 + cmw2;

  Matrix ref(2, 2, {2, 5, 7, 9});
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_addition_4)
{
  constexpr auto      num_rows1    = 2;
  constexpr auto      num_columns1 = 3;
  std::vector<double> val1         = {1, 2, 3, 4, 5, 6};

  constexpr auto      num_rows2    = 3;
  constexpr auto      num_columns2 = 2;
  std::vector<double> val2         = {1, 3, 4, 5, 5, 6};

  const Matrix_Const_Wrapper cmw1 = {num_rows1, num_columns1, val1.data()};
  const Matrix_Const_Wrapper cmw2 = {num_rows2, num_columns2, val2.data()};

  const auto cmw3 = cmw2.transposed_const_matrix_wrapper();

  const auto result = cmw1 + cmw3;

  Matrix ref(2, 3, {2, 6, 8, 7, 10, 12});
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_addition_5)
{
  constexpr auto      num_rows    = 2;
  constexpr auto      num_columns = 2;
  std::vector<double> val1        = {1, 2, 3, 4};
  std::vector<double> val2        = {1, 3, 4, 5};

  const Matrix_Const_Wrapper cmw1 = {num_rows, num_columns, val1.data()};
  const Matrix_Const_Wrapper cmw2 = {num_rows, num_columns, val2.data()};

  const auto cmw1T = cmw1.transposed_const_matrix_wrapper();
  const auto cmw2T = cmw2.transposed_const_matrix_wrapper();

  const auto result = cmw1T + cmw2T;

  Matrix ref(2, 2, {2, 7, 5, 9});
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_addition_6)
{
  constexpr auto      num_rows1    = 2;
  constexpr auto      num_columns1 = 5;
  std::vector<double> val1         = {1.2345, 2.346345, 6.3262345, 8.5674567, 6.23452345, 2.53462, 6.432452345, 2.345345, 1.234563245, 7.3245345};

  constexpr auto      num_rows2    = 5;
  constexpr auto      num_columns2 = 2;
  std::vector<double> val2         = {1.234234, 2.3462345, 345.324, 2.6345345, 634523.5, 2345345.3, 23453.345, 234534.6, 234523.5, 623452.1};

  const Matrix_Const_Wrapper cmw1   = {num_rows1, num_columns1, val1.data()};
  const Matrix_Const_Wrapper cmw2   = {num_rows2, num_columns2, val2.data()};
  const auto                 cmw2T  = cmw2.transposed_const_matrix_wrapper();
  const auto                 result = cmw1 + cmw2T;

  Matrix ref(2, 5, {2.468734, 347.670345, 634529.8262345, 23461.9124567, 234529.73452345, 4.8808545, 9.066986845, 2345347.645345, 234535.834563245, 623459.4245345});
  matrix_test_API::compare_considering_4ULP(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_addition_7)
{
  constexpr auto      num_rows1    = 3;
  constexpr auto      num_columns1 = 5;
  std::vector<double> val1         = {1.2345, 2.346345, 6.3262345, 8.5674567, 6.23452345, 2.53462, 6.432452345, 2.345345, 1.234563245, 7.3245345, 789.45978, 74.5789123, 74.23541, 4.7894113, 7894.5134};

  constexpr auto      num_rows2    = 5;
  constexpr auto      num_columns2 = 3;
  std::vector<double> val2         = {1.234234, 2.3462345, 789456.0, 345.324, 2.6345345, 74.48651, 634523.5, 2345345.3, 710.1846, 23453.345, 234534.6, 12.5487, 234523.5, 623452.1, 421.7456};

  const Matrix_Const_Wrapper cmw1   = {num_rows1, num_columns1, val1.data()};
  const Matrix_Const_Wrapper cmw2   = {num_rows2, num_columns2, val2.data()};
  const auto                 cmw2T  = cmw2.transposed_const_matrix_wrapper();
  const auto                 result = cmw1 + cmw2T;

  Matrix ref(3, 5, {2.468734, 347.670345, 634529.8262345, 23461.9124567, 234529.73452345, 4.8808545, 9.066986845, 2345347.645345, 234535.834563245, 623459.4245345, 790245.45978, 149.0654223, 784.42001, 17.3381113, 8316.259});
  matrix_test_API::compare_considering_4ULP(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_addition_8)
{
  constexpr auto             num_rows1          = 2;
  constexpr auto             num_columns1       = 3;
  constexpr auto             leading_dimension1 = 5;
  std::vector<double>        val1               = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  const Matrix_Const_Wrapper cmw1               = {num_rows1, num_columns1, val1.data(), leading_dimension1};
  // cmw1 =	¦¢ 1	2	3 ¦¢
  //				¦¢ 6	7	8 ¦¢

  constexpr auto             num_rows2          = 3;
  constexpr auto             num_columns2       = 2;
  constexpr auto             leading_dimension2 = 3;
  std::vector<double>        val2               = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  const Matrix_Const_Wrapper cmw2               = {num_rows2, num_columns2, val2.data(), leading_dimension2, Transpose_Type::transposed};
  // cmw2 =	¦¢ 1	2	¦¢T	=	¦¢ 1	4	7 ¦¢
  //				¦¢ 4	5	¦¢			¦¢ 2	5	8 ¦¢
  //				¦¢ 7	8	¦¢

  const auto   result = cmw1 + cmw2;
  const Matrix ref    = {2, 3, {2, 6, 10, 8, 12, 16}};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_addition_9)
{
  constexpr auto             num_rows1          = 2;
  constexpr auto             num_columns1       = 3;
  constexpr auto             leading_dimension1 = 5;
  std::vector<double>        val1               = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  const Matrix_Const_Wrapper cmw1               = {num_rows1, num_columns1, val1.data(), leading_dimension1};
  // cmw1 =	¦¢ 1	2	3 ¦¢
  //				¦¢ 6	7	8 ¦¢

  constexpr auto             num_rows2          = 3;
  constexpr auto             num_columns2       = 2;
  constexpr auto             leading_dimension2 = 3;
  std::vector<double>        val2               = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  const Matrix_Const_Wrapper cmw2               = {num_rows2, num_columns2, val2.data(), leading_dimension2, Transpose_Type::transposed};
  // cmw2 =	¦¢ 1	2	¦¢T	=	¦¢ 1	4	7 ¦¢
  //				¦¢ 4	5	¦¢			¦¢ 2	5	8 ¦¢
  //				¦¢ 7	8	¦¢

  const auto   result = cmw2 + cmw1;
  const Matrix ref    = {2, 3, {2, 6, 10, 8, 12, 16}};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_scalar_multiplication_1)
{
  constexpr auto num_rows    = 3;
  constexpr auto num_columns = 3;

  std::vector<double>  values(num_rows * num_columns, 1);
  Matrix_Const_Wrapper cmw(num_rows, num_columns, values.data());

  const auto result = cmw * 10;

  std::vector<double> ref_val(num_rows * num_columns, 10);

  const auto result_ptr = result.data();
  const auto ref_ptr    = ref_val.data();
  for (int i = 0; i < num_rows * num_columns; ++i)
  {
    EXPECT_EQ(result_ptr[i], ref_ptr[i]);
  }
}
TEST(Matrix_Const_Wrapper, operator_scalar_multiplication_2)
{
  constexpr auto num_rows    = 10;
  constexpr auto num_columns = 10;

  std::vector<double>  values(num_rows * num_columns, 1);
  Matrix_Const_Wrapper cmw(num_rows, num_columns, values.data());

  const auto result = cmw * 10;

  std::vector<double> ref_val(num_rows * num_columns, 10);

  const auto result_ptr = result.data();
  const auto ref_ptr    = ref_val.data();
  for (int i = 0; i < num_rows * num_columns; ++i)
  {
    EXPECT_EQ(result_ptr[i], ref_ptr[i]);
  }
}
TEST(Matrix_Const_Wrapper, operator_scalar_multiplication_3)
{
  constexpr auto num_rows          = 10;
  constexpr auto num_columns       = 10;
  constexpr auto leading_dimension = 20;

  std::vector<double>  values(num_rows * leading_dimension, 1);
  Matrix_Const_Wrapper cmw(num_rows, num_columns, values.data(), leading_dimension);

  const auto result = cmw * 10;

  std::vector<double> ref_val(num_rows * num_columns, 10);

  const auto result_ptr = result.data();
  const auto ref_ptr    = ref_val.data();
  for (int i = 0; i < num_rows * num_columns; ++i)
  {
    EXPECT_EQ(result_ptr[i], ref_ptr[i]);
  }
}
TEST(Matrix_Const_Wrapper, operator_scalar_multiplication_4)
{
  constexpr auto       num_rows    = 2;
  constexpr auto       num_columns = 3;
  std::vector<double>  values      = {1, 2, 3, 4, 5, 6};
  Matrix_Const_Wrapper cmw(num_rows, num_columns, values.data(), Transpose_Type::transposed);

  const auto result = cmw * 10;

  constexpr auto       ref_num_rows    = num_columns;
  constexpr auto       ref_num_columns = num_rows;
  std::vector<double>  ref_val         = {10, 40, 20, 50, 30, 60};
  Matrix_Const_Wrapper ref(ref_num_rows, ref_num_columns, ref_val.data());

  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mv_1)
{
  constexpr auto             num_row    = 3;
  constexpr auto             num_column = 3;
  std::vector<double>        val1       = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  const Matrix_Const_Wrapper m          = {num_row, num_column, val1.data()};

  const std::vector<double> vec = {-1, -1, 1};

  const auto result = m * vec;

  const Vector ref = {0, -3, -6};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mv_2)
{
  constexpr auto      num_row    = 2;
  constexpr auto      num_column = 2;
  std::vector<double> val1       = {1, 2, 3, 4};

  const Matrix_Const_Wrapper m = {num_row, num_column, val1.data()};

  const std::vector<double> vec = {1, 2};

  const auto result = m * vec;

  const Vector ref = {5, 11};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mv_3)
{
  constexpr auto             num_row           = 2;
  constexpr auto             num_column        = 2;
  std::vector<double>        val1              = {1, 2, 9, 3, 4, 2};
  constexpr auto             leading_dimension = 3;
  const Matrix_Const_Wrapper m                 = {num_row, num_column, val1.data(), leading_dimension};

  const std::vector<double> vec = {1, 2};

  const auto result = m * vec;

  const Vector ref = {5, 11};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mv_4)
{
  constexpr auto      num_row    = 3;
  constexpr auto      num_column = 3;
  std::vector<double> val1       = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  const Matrix_Const_Wrapper m = {num_row, num_column, val1.data(), Transpose_Type::transposed};

  const std::vector<double> vec = {-1, -1, 1};

  const auto result = m * vec;

  const Vector ref = {2, 1, 0};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mv_5)
{
  constexpr auto      num_row           = 2;
  constexpr auto      num_column        = 2;
  std::vector<double> val1              = {1, 2, 9, 3, 4, 2};
  constexpr auto      leading_dimension = 3;

  const Matrix_Const_Wrapper m      = {num_row, num_column, val1.data(), leading_dimension, Transpose_Type::transposed};
  const std::vector<double>  vec    = {1, 2};
  const auto                 result = m * vec;

  const Vector ref = {7, 10};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mv_6)
{
  constexpr auto      num_row           = 3;
  constexpr auto      num_column        = 2;
  std::vector<double> val1              = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  constexpr auto      leading_dimension = 4;

  const Matrix_Const_Wrapper m      = {num_row, num_column, val1.data(), leading_dimension, Transpose_Type::transposed};
  const std::vector<double>  v      = {-1, -1, 1};
  const auto                 result = m * v;

  const std::vector<double> ref = {3, 2};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mv_7)
{
  constexpr auto      num_row           = 3;
  constexpr auto      num_column        = 2;
  std::vector<double> val1              = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  constexpr auto      leading_dimension = 3;

  const Matrix_Const_Wrapper m      = {num_row, num_column, val1.data(), leading_dimension};
  const std::vector<double>  v      = {2, -1};
  const auto                 result = m * v;

  const std::vector<double> ref = {0, 3, 6};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mv_8)
{
  constexpr auto      num_row           = 3;
  constexpr auto      num_column        = 2;
  std::vector<double> val1              = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  constexpr auto      leading_dimension = 4;

  const Matrix_Const_Wrapper m      = {num_row, num_column, val1.data(), leading_dimension};
  const std::vector<double>  v      = {-1, 1};
  const auto                 result = m * v;

  const std::vector<double> ref = {1, 1, 1};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mv_9)
{
  constexpr auto      num_row           = 3;
  constexpr auto      num_column        = 2;
  std::vector<double> val1              = {1, 2, 3, 4, 5, 6};
  constexpr auto      leading_dimension = 2;

  const Matrix_Const_Wrapper m      = {num_row, num_column, val1.data(), leading_dimension, Transpose_Type::transposed};
  const std::vector<double>  v      = {-1, -1, 1};
  const auto                 result = m * v;

  const std::vector<double> ref = {1, 0};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mv_10)
{
  constexpr auto      num_row           = 3;
  constexpr auto      num_column        = 2;
  std::vector<double> val1              = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  constexpr auto      leading_dimension = 3;

  const Matrix_Const_Wrapper m = {num_row, num_column, val1.data(), leading_dimension, Transpose_Type::transposed};
  // m =	¦£ 1 2 ¦¤T	 =	¦¢	1 4 7	¦¢
  //			¦¢ 4 5 ¦¢				¦¢	2 5 8	¦¢
  //			¦¦ 7 8 ¦¥

  const std::vector<double> v      = {-1, -1, 1};
  const auto                result = m * v;

  const std::vector<double> ref = {2, 1};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mm_1)
{
  constexpr auto      A_num_rows          = 2;
  constexpr auto      A_num_columns       = 2;
  constexpr auto      A_leading_dimension = 2;
  std::vector<double> A_vals              = {1, 2, 3, 4};

  constexpr auto      B_num_rows          = 2;
  constexpr auto      B_num_columns       = 2;
  constexpr auto      B_leading_dimension = 2;
  std::vector<double> B_vals              = {1, 2, 3, 4};

  const Matrix_Const_Wrapper A      = {A_num_rows, A_num_columns, A_vals.data(), A_leading_dimension};
  const Matrix_Const_Wrapper B      = {B_num_rows, B_num_columns, B_vals.data(), B_leading_dimension};
  const auto                 result = A * B;

  Matrix ref = {2, 2, {7, 10, 15, 22}};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mm_2)
{
  constexpr auto      A_num_rows          = 2;
  constexpr auto      A_num_columns       = 8;
  constexpr auto      A_leading_dimension = 8;
  std::vector<double> A_vals              = {1, 2, 3, 4, 5, 6, 7, 8, 8, 7, 6, 5, 4, 3, 2, 1};

  constexpr auto      B_num_rows          = 8;
  constexpr auto      B_num_columns       = 2;
  constexpr auto      B_leading_dimension = 2;
  std::vector<double> B_vals              = {1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8};

  const Matrix_Const_Wrapper A      = {A_num_rows, A_num_columns, A_vals.data(), A_leading_dimension};
  const Matrix_Const_Wrapper B      = {B_num_rows, B_num_columns, B_vals.data(), B_leading_dimension};
  const auto                 result = A * B;

  Matrix ref = {2, 2, {204, 204, 120, 120}};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mm_3)
{
  constexpr auto      A_num_rows          = 2;
  constexpr auto      A_num_columns       = 8;
  constexpr auto      A_leading_dimension = 10;
  std::vector<double> A_vals              = {1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 8, 7, 6, 5, 4, 3, 2, 1, 0, 0};

  constexpr auto      B_num_rows          = 8;
  constexpr auto      B_num_columns       = 2;
  constexpr auto      B_leading_dimension = 3;
  std::vector<double> B_vals              = {1, 1, 0, 2, 2, 0, 3, 3, 0, 4, 4, 0, 5, 5, 0, 6, 6, 0, 7, 7, 0, 8, 8, 0};

  const Matrix_Const_Wrapper A      = {A_num_rows, A_num_columns, A_vals.data(), A_leading_dimension};
  const Matrix_Const_Wrapper B      = {B_num_rows, B_num_columns, B_vals.data(), B_leading_dimension};
  const auto                 result = A * B;

  Matrix ref = {2, 2, {204, 204, 120, 120}};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mm_4)
{
  constexpr auto      A_num_rows          = 2;
  constexpr auto      A_num_columns       = 2;
  constexpr auto      A_leading_dimension = 2;
  std::vector<double> A_vals              = {1, 2, 3, 4};

  constexpr auto      B_num_rows          = 2;
  constexpr auto      B_num_columns       = 2;
  constexpr auto      B_leading_dimension = 2;
  std::vector<double> B_vals              = {1, 3, 4, 5};

  const Matrix_Const_Wrapper A      = {A_num_rows, A_num_columns, A_vals.data(), A_leading_dimension};
  const Matrix_Const_Wrapper B      = {B_num_rows, B_num_columns, B_vals.data(), B_leading_dimension};
  const auto                 AT     = A.transposed_const_matrix_wrapper();
  const auto                 result = AT * B;

  Matrix ref = {2, 2, {13, 18, 18, 26}};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mm_5)
{
  constexpr auto      A_num_rows          = 2;
  constexpr auto      A_num_columns       = 2;
  constexpr auto      A_leading_dimension = 2;
  std::vector<double> A_vals              = {1, 3, 4, 5};

  constexpr auto      B_num_rows          = 2;
  constexpr auto      B_num_columns       = 2;
  constexpr auto      B_leading_dimension = 2;
  std::vector<double> B_vals              = {1, 2, 3, 4};

  const Matrix_Const_Wrapper A      = {A_num_rows, A_num_columns, A_vals.data(), A_leading_dimension};
  const Matrix_Const_Wrapper B      = {B_num_rows, B_num_columns, B_vals.data(), B_leading_dimension};
  const auto                 AT     = A.transposed_const_matrix_wrapper();
  const auto                 result = AT * B;

  Matrix ref = {2, 2, {13, 18, 18, 26}};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mm_6)
{
  constexpr auto      A_num_rows          = 2;
  constexpr auto      A_num_columns       = 3;
  constexpr auto      A_leading_dimension = 3;
  std::vector<double> A_vals              = {1, 2, 3, 4, 5, 6};

  constexpr auto      B_num_rows          = 2;
  constexpr auto      B_num_columns       = 3;
  constexpr auto      B_leading_dimension = 3;
  std::vector<double> B_vals              = {1, 2, 3, 3, 2, 1};

  const Matrix_Const_Wrapper A      = {A_num_rows, A_num_columns, A_vals.data(), A_leading_dimension};
  const Matrix_Const_Wrapper B      = {B_num_rows, B_num_columns, B_vals.data(), B_leading_dimension};
  const auto                 AT     = A.transposed_const_matrix_wrapper();
  const auto                 result = AT * B;

  Matrix ref = {3, 3, {13, 10, 7, 17, 14, 11, 21, 18, 15}};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mm_7)
{
  constexpr auto      A_num_rows          = 2;
  constexpr auto      A_num_columns       = 3;
  constexpr auto      A_leading_dimension = 3;
  std::vector<double> A_vals              = {1.5479, 2.4567123, 3.414878, 4.41487, 5.121, 6.15789};

  constexpr auto      B_num_rows          = 2;
  constexpr auto      B_num_columns       = 3;
  constexpr auto      B_leading_dimension = 3;
  std::vector<double> B_vals              = {1.1244, 2.48711, 3.12314, 3.789413, 2.9135491, 1.264863};

  const Matrix_Const_Wrapper A      = {A_num_rows, A_num_columns, A_vals.data(), A_leading_dimension};
  const Matrix_Const_Wrapper B      = {B_num_rows, B_num_columns, B_vals.data(), B_leading_dimension};
  const auto                 AT     = A.transposed_const_matrix_wrapper();
  const auto                 result = AT * B;

  Matrix ref = {3, 3, {18.470224531309999, 16.712738084116999, 10.418514118810000, 22.167911283120002, 21.030398669553001, 14.150019875622000, 27.174477241769999, 26.434492089979003, 18.454029295989997}};
  matrix_test_API::compare_considering_4ULP(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mm_8)
{
  constexpr auto      A_num_rows          = 2;
  constexpr auto      A_num_columns       = 3;
  constexpr auto      A_leading_dimension = 5;
  std::vector<double> A_vals              = {1.5479, 2.4567123, 3.414878, 0, 0, 4.41487, 5.121, 6.15789};

  constexpr auto      B_num_rows          = 2;
  constexpr auto      B_num_columns       = 3;
  constexpr auto      B_leading_dimension = 4;
  std::vector<double> B_vals              = {1.1244, 2.48711, 3.12314, 0, 3.789413, 2.9135491, 1.264863};

  const Matrix_Const_Wrapper A      = {A_num_rows, A_num_columns, A_vals.data(), A_leading_dimension};
  const Matrix_Const_Wrapper B      = {B_num_rows, B_num_columns, B_vals.data(), B_leading_dimension};
  const auto                 AT     = A.transposed_const_matrix_wrapper();
  const auto                 result = AT * B;

  Matrix ref = {3, 3, {18.470224531309999, 16.712738084116999, 10.418514118810000, 22.167911283120002, 21.030398669553001, 14.150019875622000, 27.174477241769999, 26.434492089979003, 18.454029295989997}};
  matrix_test_API::compare_considering_4ULP(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mm_9)
{
  constexpr auto      A_num_rows          = 2;
  constexpr auto      A_num_columns       = 2;
  constexpr auto      A_leading_dimension = 2;
  std::vector<double> A_vals              = {1, 2, 3, 4};

  constexpr auto      B_num_rows          = 2;
  constexpr auto      B_num_columns       = 2;
  constexpr auto      B_leading_dimension = 2;
  std::vector<double> B_vals              = {1, 3, 4, 5};

  const Matrix_Const_Wrapper A      = {A_num_rows, A_num_columns, A_vals.data(), A_leading_dimension};
  const Matrix_Const_Wrapper B      = {B_num_rows, B_num_columns, B_vals.data(), B_leading_dimension};
  const auto                 BT     = B.transposed_const_matrix_wrapper();
  const auto                 result = A * BT;

  Matrix ref = {2, 2, {7, 14, 15, 32}};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mm_10)
{
  constexpr auto      A_num_rows          = 2;
  constexpr auto      A_num_columns       = 3;
  constexpr auto      A_leading_dimension = 3;
  std::vector<double> A_vals              = {1, 2, 3, 4, 5, 6};

  constexpr auto      B_num_rows          = 2;
  constexpr auto      B_num_columns       = 3;
  constexpr auto      B_leading_dimension = 3;
  std::vector<double> B_vals              = {1, 2, 3, 3, 2, 1};

  const Matrix_Const_Wrapper A      = {A_num_rows, A_num_columns, A_vals.data(), A_leading_dimension};
  const Matrix_Const_Wrapper B      = {B_num_rows, B_num_columns, B_vals.data(), B_leading_dimension};
  const auto                 BT     = B.transposed_const_matrix_wrapper();
  const auto                 result = A * BT;

  Matrix ref = {2, 2, {14, 10, 32, 28}};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mm_11)
{
  constexpr auto      A_num_rows          = 2;
  constexpr auto      A_num_columns       = 3;
  constexpr auto      A_leading_dimension = 3;
  std::vector<double> A_vals              = {1.5479, 2.4567123, 3.414878, 4.41487, 5.121, 6.15789};

  constexpr auto      B_num_rows          = 2;
  constexpr auto      B_num_columns       = 3;
  constexpr auto      B_leading_dimension = 3;
  std::vector<double> B_vals              = {1.1244, 2.48711, 3.12314, 3.789413, 2.9135491, 1.264863};

  const Matrix_Const_Wrapper A      = {A_num_rows, A_num_columns, A_vals.data(), A_leading_dimension};
  const Matrix_Const_Wrapper B      = {B_num_rows, B_num_columns, B_vals.data(), B_leading_dimension};
  const auto                 BT     = B.transposed_const_matrix_wrapper();
  const auto                 result = A * BT;

  Matrix ref = {2, 2, {18.515714565372999, 17.342737125037932, 36.932522712599997, 39.438937931479998}};
  matrix_test_API::compare_considering_4ULP(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mm_12)
{
  constexpr auto      A_num_rows          = 2;
  constexpr auto      A_num_columns       = 3;
  constexpr auto      A_leading_dimension = 5;
  std::vector<double> A_vals              = {1.5479, 2.4567123, 3.414878, 0, 0, 4.41487, 5.121, 6.15789};

  constexpr auto      B_num_rows          = 2;
  constexpr auto      B_num_columns       = 3;
  constexpr auto      B_leading_dimension = 5;
  std::vector<double> B_vals              = {1.1244, 2.48711, 3.12314, 0, 0, 3.789413, 2.9135491, 1.264863};

  const Matrix_Const_Wrapper A      = {A_num_rows, A_num_columns, A_vals.data(), A_leading_dimension};
  const Matrix_Const_Wrapper B      = {B_num_rows, B_num_columns, B_vals.data(), B_leading_dimension};
  const auto                 BT     = B.transposed_const_matrix_wrapper();
  const auto                 result = A * BT;

  Matrix ref = {2, 2, {18.515714565372999, 17.342737125037932, 36.932522712599997, 39.438937931479998}};
  matrix_test_API::compare_considering_4ULP(result, ref);
}

TEST(Matrix_Const_Wrapper, operator_mm_13)
{
  constexpr auto      A_num_rows          = 2;
  constexpr auto      A_num_columns       = 2;
  constexpr auto      A_leading_dimension = 2;
  std::vector<double> A_vals              = {1, 2, 3, 4};

  constexpr auto      B_num_rows          = 2;
  constexpr auto      B_num_columns       = 2;
  constexpr auto      B_leading_dimension = 2;
  std::vector<double> B_vals              = {1, 3, 4, 5};

  const Matrix_Const_Wrapper A      = {A_num_rows, A_num_columns, A_vals.data(), A_leading_dimension};
  const Matrix_Const_Wrapper B      = {B_num_rows, B_num_columns, B_vals.data(), B_leading_dimension};
  const auto                 AT     = A.transposed_const_matrix_wrapper();
  const auto                 BT     = B.transposed_const_matrix_wrapper();
  const auto                 result = AT * BT;

  Matrix ref = {2, 2, {10, 19, 14, 28}};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mm_14)
{
  constexpr auto      A_num_rows          = 2;
  constexpr auto      A_num_columns       = 3;
  constexpr auto      A_leading_dimension = 3;
  std::vector<double> A_vals              = {1, 2, 3, 4, 5, 6};

  constexpr auto      B_num_rows          = 3;
  constexpr auto      B_num_columns       = 2;
  constexpr auto      B_leading_dimension = 2;
  std::vector<double> B_vals              = {1, 2, 3, 3, 2, 1};

  const Matrix_Const_Wrapper A      = {A_num_rows, A_num_columns, A_vals.data(), A_leading_dimension};
  const Matrix_Const_Wrapper B      = {B_num_rows, B_num_columns, B_vals.data(), B_leading_dimension};
  const auto                 AT     = A.transposed_const_matrix_wrapper();
  const auto                 BT     = B.transposed_const_matrix_wrapper();
  const auto                 result = AT * BT;

  Matrix ref = {3, 3, {9, 15, 6, 12, 21, 9, 15, 27, 12}};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mm_15)
{
  constexpr auto      A_num_rows          = 2;
  constexpr auto      A_num_columns       = 3;
  constexpr auto      A_leading_dimension = 4;
  std::vector<double> A_vals              = {1, 2, 3, 0, 4, 5, 6};

  constexpr auto      B_num_rows          = 3;
  constexpr auto      B_num_columns       = 2;
  constexpr auto      B_leading_dimension = 4;
  std::vector<double> B_vals              = {1, 2, 0, 0, 3, 3, 0, 0, 2, 1};

  const Matrix_Const_Wrapper A      = {A_num_rows, A_num_columns, A_vals.data(), A_leading_dimension};
  const Matrix_Const_Wrapper B      = {B_num_rows, B_num_columns, B_vals.data(), B_leading_dimension};
  const auto                 AT     = A.transposed_const_matrix_wrapper();
  const auto                 BT     = B.transposed_const_matrix_wrapper();
  const auto                 result = AT * BT;

  Matrix ref = {3, 3, {9, 15, 6, 12, 21, 9, 15, 27, 12}};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, operator_mm_16)
{
  constexpr auto             num_rows1          = 1;
  constexpr auto             num_columns1       = 2;
  constexpr auto             leading_dimension1 = 2;
  std::vector<double>        val1               = {1, 2};
  const Matrix_Const_Wrapper m1                 = {num_rows1, num_columns1, val1.data(), leading_dimension1};

  constexpr auto             num_rows2          = 2;
  constexpr auto             num_columns2       = 2;
  constexpr auto             leading_dimension2 = 2;
  std::vector<double>        val2               = {2, 2, 2, 2};
  const Matrix_Const_Wrapper m2                 = {num_rows2, num_columns2, val2.data(), leading_dimension2};

  const auto result = m1 * m2;

  const Matrix ref = {1, 2, {6, 6}};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, at_1)
{
  constexpr auto      num_row           = 3;
  constexpr auto      num_column        = 2;
  std::vector<double> val1              = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  constexpr auto      leading_dimension = 3;

  const Matrix_Const_Wrapper m = {num_row, num_column, val1.data(), leading_dimension, Transpose_Type::transposed};
  // m =	¦£ 1 2 ¦¤T	 =	¦¢	1 4 7	¦¢
  //			¦¢ 4 5 ¦¢				¦¢	2 5 8	¦¢
  //			¦¦ 7 8 ¦¥

  std::vector<double> ref = {1, 4, 7, 2, 5, 8};

  int index = 0;
  for (int i = 0; i < m.num_rows(); i++)
  {
    for (int j = 0; j < m.num_cols(); j++)
    {
      EXPECT_EQ(m.at(i, j), ref[index]);
      index++;
    }
  }
}
TEST(Matrix_Const_Wrapper, const_column_vector_wrapper_1)
{
  constexpr auto             num_row           = 2;
  constexpr auto             num_column        = 2;
  std::vector<double>        val1              = {1, 2, 9, 3, 4, 2};
  constexpr auto             leading_dimension = 3;
  const Matrix_Const_Wrapper m                 = {num_row, num_column, val1.data(), leading_dimension};

  constexpr auto column_index = 0;
  const auto     result       = m.const_column_vector_wrapper(column_index);

  std::vector<double> ref = {1, 3};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, const_column_vector_wrapper_2)
{
  constexpr auto      num_row           = 3;
  constexpr auto      num_column        = 2;
  std::vector<double> val1              = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  constexpr auto      leading_dimension = 3;

  const Matrix_Const_Wrapper m = {num_row, num_column, val1.data(), leading_dimension, Transpose_Type::transposed};
  // m =	¦£ 1 2 ¦¤T	 =	¦¢	1 4 7	¦¢
  //			¦¢ 4 5 ¦¢				¦¢	2 5 8	¦¢
  //			¦¦ 7 8 ¦¥

  constexpr auto column_index = 0;
  const auto     result       = m.const_column_vector_wrapper(column_index);

  std::vector<double> ref = {1, 2};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, const_row_vector_wrapper_1)
{
  constexpr auto             num_row           = 2;
  constexpr auto             num_column        = 2;
  std::vector<double>        val1              = {1, 2, 9, 3, 4, 2};
  constexpr auto             leading_dimension = 3;
  const Matrix_Const_Wrapper m                 = {num_row, num_column, val1.data(), leading_dimension};

  constexpr auto row_index = 0;
  const auto     result    = m.const_row_vector_wrapper(row_index);

  std::vector<double> ref = {1, 2};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, const_row_vector_wrapper_2)
{
  constexpr auto      num_row           = 3;
  constexpr auto      num_column        = 2;
  std::vector<double> val1              = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  constexpr auto      leading_dimension = 3;

  const Matrix_Const_Wrapper m = {num_row, num_column, val1.data(), leading_dimension, Transpose_Type::transposed};
  // m =	¦£ 1 2 ¦¤T	 =	¦¢	1 4 7	¦¢
  //			¦¢ 4 5 ¦¢				¦¢	2 5 8	¦¢
  //			¦¦ 7 8 ¦¥

  constexpr auto row_index = 0;
  const auto     result    = m.const_row_vector_wrapper(row_index);

  std::vector<double> ref = {1, 4, 7};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, const_row_vector_wrapper_3)
{
  constexpr auto      num_row           = 3;
  constexpr auto      num_column        = 2;
  std::vector<double> val1              = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  constexpr auto      leading_dimension = 3;

  const Matrix_Const_Wrapper m = {num_row, num_column, val1.data(), leading_dimension, Transpose_Type::transposed};
  // m =	¦£ 1 2 ¦¤T	 =	¦¢	1 4 7	¦¢
  //			¦¢ 4 5 ¦¢				¦¢	2 5 8	¦¢
  //			¦¦ 7 8 ¦¥

  constexpr auto row_index = 1;
  const auto     result    = m.const_row_vector_wrapper(row_index);

  std::vector<double> ref = {2, 5, 8};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, const_part_1)
{
  constexpr auto       num_row    = 2;
  constexpr auto       num_column = 4;
  std::vector<double>  val        = {1, 2, 3, 4, 5, 6, 7, 8};
  Matrix_Const_Wrapper m          = {num_row, num_column, val.data()};
  // m =	¦£ 1 2 3 4 ¦¤
  //			¦¦ 5 6 7 8 ¦¥

  constexpr auto start_row_index    = 0;
  constexpr auto end_row_index      = 2;
  constexpr auto start_column_index = 1;
  constexpr auto end_column_index   = 3;
  auto           result             = m.const_part(start_row_index, end_row_index, start_column_index, end_column_index);

  constexpr auto       ref_num_row    = 2;
  constexpr auto       ref_num_column = 2;
  std::vector<double>  ref_val        = {2, 3, 6, 7};
  Matrix_Const_Wrapper ref            = {ref_num_row, ref_num_column, ref_val.data()};

  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, const_part_2)
{
  constexpr auto       num_row    = 2;
  constexpr auto       num_column = 4;
  std::vector<double>  val        = {1, 2, 3, 4, 5, 6, 7, 8};
  Matrix_Const_Wrapper m          = {num_row, num_column, val.data()};
  // m =	¦£ 1 2 3 4 ¦¤
  //			¦¦ 5 6 7 8 ¦¥

  constexpr auto start_row_index    = 0;
  constexpr auto end_row_index      = 2;
  constexpr auto start_column_index = 0;
  constexpr auto end_column_index   = 3;
  auto           result             = m.const_part(start_row_index, end_row_index, start_column_index, end_column_index);

  constexpr auto       ref_num_row    = 2;
  constexpr auto       ref_num_column = 3;
  std::vector<double>  ref_val        = {1, 2, 3, 5, 6, 7};
  Matrix_Const_Wrapper ref            = {ref_num_row, ref_num_column, ref_val.data()};

  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, const_part_3)
{
  constexpr auto      num_row    = 2;
  constexpr auto      num_column = 4;
  std::vector<double> val        = {1, 2, 3, 4, 5, 6, 7, 8};

  Matrix_Const_Wrapper m = {num_row, num_column, val.data(), Transpose_Type::transposed};
  // m =	¦£ 1 2 3 4 ¦¤T	 =	¦¢	1 5 ¦¢
  //			¦¦ 5 6 7 8 ¦¥				¦¢	2 6 ¦¢
  //												¦¢ 3 7 ¦¢
  //												¦¢ 4 8 ¦¢

  constexpr auto start_row_index    = 1;
  constexpr auto end_row_index      = 3;
  constexpr auto start_column_index = 0;
  constexpr auto end_column_index   = 2;
  auto           result             = m.const_part(start_row_index, end_row_index, start_column_index, end_column_index);

  constexpr auto       ref_num_row    = 2;
  constexpr auto       ref_num_column = 2;
  std::vector<double>  ref_val        = {2, 6, 3, 7};
  Matrix_Const_Wrapper ref            = {ref_num_row, ref_num_column, ref_val.data()};

  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, const_part_4)
{
  constexpr auto      num_row    = 2;
  constexpr auto      num_column = 4;
  std::vector<double> val        = {1, 2, 3, 4, 5, 6, 7, 8};

  Matrix_Const_Wrapper m = {num_row, num_column, val.data(), Transpose_Type::transposed};
  // m =	¦£ 1 2 3 4 ¦¤T	 =	¦¢	1 5 ¦¢
  //			¦¦ 5 6 7 8 ¦¥				¦¢	2 6 ¦¢
  //												¦¢ 3 7 ¦¢
  //												¦¢ 4 8 ¦¢

  constexpr auto start_row_index    = 1;
  constexpr auto end_row_index      = 4;
  constexpr auto start_column_index = 0;
  constexpr auto end_column_index   = 2;
  auto           result             = m.const_part(start_row_index, end_row_index, start_column_index, end_column_index);

  constexpr auto       ref_num_row    = 3;
  constexpr auto       ref_num_column = 2;
  std::vector<double>  ref_val        = {2, 6, 3, 7, 4, 8};
  Matrix_Const_Wrapper ref            = {ref_num_row, ref_num_column, ref_val.data()};

  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, const_part_5)
{
  constexpr auto      num_row           = 3;
  constexpr auto      num_column        = 4;
  std::vector<double> val               = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  constexpr auto      leading_dimension = 5;

  Matrix_Const_Wrapper m = {num_row, num_column, val.data(), leading_dimension, Transpose_Type::transposed};
  // m =	¦¢ 1  2  3  4  ¦¢T	 =	¦¢	1 6 11 ¦¢
  //			¦¢ 6  7  8  9  ¦¢				¦¢	2 7 12 ¦¢
  //			¦¢ 11 12 13 14	¦¢				¦¢ 3	8 13 ¦¢
  //														¦¢ 4 9 14 ¦¢

  constexpr auto start_row_index    = 0;
  constexpr auto end_row_index      = 3;
  constexpr auto start_column_index = 1;
  constexpr auto end_column_index   = 3;
  auto           result             = m.const_part(start_row_index, end_row_index, start_column_index, end_column_index);

  constexpr auto       ref_num_row    = 3;
  constexpr auto       ref_num_column = 2;
  std::vector<double>  ref_val        = {6, 11, 7, 12, 8, 13};
  Matrix_Const_Wrapper ref            = {ref_num_row, ref_num_column, ref_val.data()};

  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, const_part_6)
{
  constexpr auto      num_row           = 3;
  constexpr auto      num_column        = 4;
  std::vector<double> val               = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  constexpr auto      leading_dimension = 5;

  Matrix_Const_Wrapper m = {num_row, num_column, val.data(), leading_dimension, Transpose_Type::transposed};
  // m =	¦¢ 1  2  3  4  ¦¢T	 =	¦¢	1 6 11 ¦¢
  //			¦¢ 6  7  8  9  ¦¢				¦¢	2 7 12 ¦¢
  //			¦¢ 11 12 13 14	¦¢				¦¢ 3	8 13 ¦¢
  //														¦¢ 4 9 14 ¦¢

  constexpr auto start_row_index    = 2;
  constexpr auto end_row_index      = 4;
  constexpr auto start_column_index = 0;
  constexpr auto end_column_index   = 3;
  auto           result             = m.const_part(start_row_index, end_row_index, start_column_index, end_column_index);

  constexpr auto       ref_num_row    = 2;
  constexpr auto       ref_num_column = 3;
  std::vector<double>  ref_val        = {3, 8, 13, 4, 9, 14};
  Matrix_Const_Wrapper ref            = {ref_num_row, ref_num_column, ref_val.data()};

  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, inverse_matrix_1)
{
  constexpr auto       num_rows          = 5;
  constexpr auto       num_columns       = 5;
  constexpr auto       leading_dimension = 7;
  std::vector<double>  val               = {1, 2, 3, 4, 5, 0, 0, 2, 3, 4, 1, 6, 0, 0, 1, 6, 5, 4, 3, 0, 0, 4, 1, 2, 7, 8, 0, 0, 9, 8, 7, 6, 5, 0, 0};
  Matrix_Const_Wrapper cmw               = {num_rows, num_columns, val.data(), leading_dimension};

  const auto result = cmw.inverse_matrix();

  std::vector<double> ref_val = {-0.033333333333333, -0.000000000000000, -0.166666666666667, -0.000000000000000, 0.133333333333333, -1.291666666666667, 0.250000000000000, 0.666666666666667, 0.500000000000000, -0.208333333333333, 1.616666666666666, -0.250000000000000, -0.666666666666667, -0.750000000000000, 0.283333333333333, 0.275000000000000, -0.250000000000000, 0.000000000000000, 0.000000000000000, 0.025000000000000, -0.466666666666667, 0.250000000000000, 0.166666666666667, 0.250000000000000, -0.133333333333333};
  Matrix              ref     = {num_rows, num_columns, std::move(ref_val)};

  const auto epsilon = 1.0E-15;
  matrix_test_API::compare_considering_epsilon(result, ref, epsilon);
}
TEST(Matrix_Const_Wrapper, inverse_matrix_2)
{
  constexpr auto       num_rows          = 2;
  constexpr auto       num_columns       = 2;
  constexpr auto       leading_dimension = 3;
  std::vector<double>  val               = {2, 1, 3, 1, 1, 3};
  Matrix_Const_Wrapper cmw               = {num_rows, num_columns, val.data(), leading_dimension};

  const auto result = cmw.inverse_matrix();

  std::vector<double> ref_val = {1, -1, -1, 2};
  Matrix              ref     = {num_rows, num_columns, std::move(ref_val)};

  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, transposed_matrix_1)
{
  constexpr auto      num_row           = 3;
  constexpr auto      num_column        = 4;
  std::vector<double> val               = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  constexpr auto      leading_dimension = 5;

  Matrix_Const_Wrapper m = {num_row, num_column, val.data(), leading_dimension};
  // m =	¦¢ 1  2  3  4  ¦¢
  //			¦¢ 6  7  8  9  ¦¢
  //			¦¢ 11 12 13 14	¦¢

  const auto result = m.transposed_matrix();
  //	mT =	¦¢	1 6 11 ¦¢
  //				¦¢	2 7 12 ¦¢
  //				¦¢ 3	8 13 ¦¢
  //				¦¢ 4 9 14 ¦¢

  Matrix ref = {4, 3, {1, 6, 11, 2, 7, 12, 3, 8, 13, 4, 9, 14}};
  EXPECT_EQ(result, ref);
}
TEST(Matrix_Const_Wrapper, transposed_matrix_2)
{
  constexpr auto      num_row           = 3;
  constexpr auto      num_column        = 4;
  std::vector<double> val               = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  constexpr auto      leading_dimension = 5;

  Matrix_Const_Wrapper m = {num_row, num_column, val.data(), leading_dimension, Transpose_Type::transposed};
  // m =	¦¢ 1  2  3  4  ¦¢T =	¦¢	1 6 11 ¦¢
  //			¦¢ 6  7  8  9  ¦¢			¦¢	2 7 12 ¦¢
  //			¦¢ 11 12 13 14	¦¢			¦¢ 3	8 13 ¦¢
  //													¦¢ 4 9 14 ¦¢

  const auto result = m.transposed_matrix();
  //	mT =	¦¢ 1  2  3  4  ¦¢
  //				¦¢ 6  7  8  9  ¦¢
  //				¦¢ 11 12 13 14	¦¢

  Matrix ref = {3, 4, {1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14}};
  EXPECT_EQ(result, ref);
}

// Caution
#ifdef _DEBUG

TEST(Matrix_Const_Wrapper, operator_addition_exception_1)
{
  Matrix m1(2, 3, {1, 2, 3, 4, 5, 6});
  Matrix m2(3, 2, {1, 3, 4, 5, 5, 6});
  EXPECT_ANY_THROW(m1 + m2);
}

#endif
//
}



//TEST(Const_Matrix_Wrapper, column_1)
//{
//	const Matrix m = { 2,3,{ 1,2,3,4,5,6 } };
//	const auto result = m.column(0);
//
//	const std::vector<double> ref = { 1,4 };
//	EXPECT_EQ(result, ref);
//}
//TEST(Const_Matrix_Wrapper, column_2)
//{
//	const Matrix m = { 2,3,{ 1,2,3,4,5,6 } };
//	const auto result = m.column(2);
//
//	const std::vector<double> ref = { 3,6 };
//	EXPECT_EQ(result, ref);
//}














//TEST(Const_Matrix_Wrapper, operator_mv_3)
//{
//	constexpr auto num_row = 1;
//	constexpr auto num_column = 2;
//	std::vector<double> val1 = { 1,1 };
//
//	const Const_Matrix_Wrapper m = { num_row, num_column, val1.data() };
//
//	const std::vector<double> vec = { 1,2 };
//	const Const_Real_Vector_Wrapper cvw = vec;
//
//	const auto result = (m + m) * vec;
//
//	const Real_Vector ref = { 6 };
//	EXPECT_EQ(result, ref);
//}
//TEST(Const_Matrix_Wrapper, operator_substitution_1)
//{
//	Matrix m = { 1,2,{ 1,1 } };
//	Matrix result;
//
//	result = std::move(m);
//	result.at(0, 1) = 2;
//
//	const Matrix ref = { 1,2,{1,2} };
//	EXPECT_EQ(result, ref);
//}
//TEST(Const_Matrix_Wrapper, row1)
//{
//	Matrix dm(2, 3, { 1,2,3,4,5,6 });
//	const auto result = dm.row(0);
//
//	const std::vector<double> ref = { 1,2,3 };
//	EXPECT_EQ(result, ref);
//}

//TEST(Const_Matrix_Wrapper, scalar_multiplication_at_columns_1)
//{
//	Matrix m = { 3,3,{ 1,2,3,4,5,6,7,8,9 } };
//	m.scalar_multiplcation_at_columns(0, 1, 3);
//
//	const Matrix ref = { 3,3,{3,2,3,12,5,6,21,8,9} };
//	EXPECT_EQ(m, ref);
//}
//TEST(Const_Matrix_Wrapper, scalar_multiplication_at_columns_2)
//{
//	Matrix m = { 3,3,{ 1,2,3,4,5,6,7,8,9 } };
//	m.scalar_multiplcation_at_columns(1, 3, 3);
//
//	const Matrix ref = { 3,3,{1,6,9,4,15,18,7,24,27} };
//	EXPECT_EQ(m, ref);
//}











//============================================================================================================

//TEST(Const_Matrix_Wrapper, construct_diagonal_matrix1)
//{
//	Matrix m(2, { 3, 3 });
//	Matrix ref(2, 2, { 3,0,0,3 });
//
//	EXPECT_EQ(m, ref);
//}

//TEST(Const_Matrix_Wrapper, change_row1)
//{
//	Matrix m(2, 2, { 1,2,3,4 });
//
//	std::vector<double> v = { 5,6 };
//	m.change_row(0, v);
//
//	Matrix ref(2, 2, { 5, 6, 3, 4 });
//	EXPECT_EQ(m, ref);
//}
//TEST(Const_Matrix_Wrapper, change_row2)
//{
//	Matrix m(2, 2, { 1,2,3,4 });
//
//	std::vector<double> v = { 5,6 };
//	m.change_row(0, v);
//
//	Matrix ref(2, 2, { 5, 6, 3, 4 });
//	EXPECT_EQ(m, ref);
//}
//TEST(Const_Matrix_Wrapper, change_rows1)
//{
//	const Matrix m = { 2, 2,  {1,2,3,4} };
//
//	Matrix dm(4, 2);
//	dm.change_rows(0, m);
//
//	const Matrix ref(4, 2, { 1,2,3,4,0,0,0,0 });
//	EXPECT_EQ(dm, ref);
//}
//TEST(Const_Matrix_Wrapper, change_columns_1)
//{
//	const Matrix m = { 2,2,{ 1,2,3,4 } };
//
//	Matrix dm(2, 4);
//	dm.change_columns(0, m);
//
//	const Matrix ref(2, 4, { 1,2,0,0,3,4,0,0 });
//	EXPECT_EQ(dm, ref);
//}
//TEST(Const_Matrix_Wrapper, change_columns_2)
//{
//	Matrix m = { 2,3,{ 1,2,3,4,5,6 } };
//	const Matrix m1 = { 2,2,{ 1,1,4,4 } };
//	m.change_columns(1, m1);
//
//	const Matrix ref = { 2,3,{ 1,1,1,4,4,4 } };
//	EXPECT_EQ(m, ref);
//}
//TEST(Const_Matrix_Wrapper, change_columns_3)
//{
//	Matrix m = { 2,3,{ 1,2,3,4,5,6 } };
//	m.change_columns(0, 1, 4);
//
//	const Matrix ref = { 2,3,{ 4,2,3,4,5,6 } };
//	EXPECT_EQ(m, ref);
//}


//TEST(Const_Matrix_Wrapper_Wrapper, change_columns_1)
//{
//	std::vector<double> value = { 1,2,3,4,5,6 };
//	Matrix_Wrapper mw(2, 2, value.data() + 2);
//	mw.change_columns(0, 2, 0.0);
//	std::vector<double> ref = { 1,2,0,0,0,0 };
//	EXPECT_EQ(value, ref);
//}
//TEST(Const_Matrix_Wrapper_Wrapper, change_columns_2)
//{
//	std::vector<double> value = { 1,2,3,4,5,6 };
//	Matrix_Wrapper mw(2, 3, value.data());
//	mw.change_columns(1, 3, 0.0);
//	std::vector<double> ref = { 1,0,0,4,0,0 };
//	EXPECT_EQ(value, ref);
//}


//TEST(Const_Matrix_Wrapper_Function, change_column_1)
//{
//	Polynomial x("x0");
//	Polynomial y("x1");
//
//	Vector_Function<Polynomial> vf = { x + y , 2 * x + y };
//
//	Matrix_Function<Polynomial> result(2, 2);
//	result.change_column(0, vf[0].gradient());
//	result.change_column(1, vf[1].gradient());
//
//	Matrix_Function<Polynomial> ref(2, 2, { 1,2,1,1 });
//	EXPECT_EQ(result, ref);
//}
//TEST(ms, Jacobian_1) 
//{
//	constexpr ushort domain_dimension = 2;
//
//	Polynomial x("x0");
//	Polynomial y("x1");
//
//	Vector_Function<Polynomial, domain_dimension> vf = { x * y , 2 * x + y };
//
//	const auto result = ms::Jacobian(vf);
//
//	Matrix_Function<Polynomial, domain_dimension, domain_dimension> ref = { y,x,2,1 };
//	EXPECT_EQ(result, ref);
//}
//TEST(Vector_Function, mv_1) 
//{
//	constexpr ushort domain_dimension = 2;
//
//	Polynomial x("x0");
//	Polynomial y("x1");
//
//	Static_Matrix<2, 2> m = { 1,2,3,4 };
//	Vector_Function<Polynomial<2>, 2> vf = { x , y };
//	const auto result = m * vf;
//
//	Vector_Function<Polynomial<2>, 2> ref = { x + 2 * y ,3 * x + 4 * y };
//	EXPECT_EQ(result, ref);
//}