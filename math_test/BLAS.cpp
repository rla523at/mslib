#include "API_BLAS.h"


TEST(MSBLAS, abs_x_1)
{
	std::vector<double> x = { -1,-2,3 };

	const auto n = static_cast<int>(x.size());
	ms::math::blas::abs_x(x.data(), n); 

	const std::vector<double> ref = { 1,2,3 };
	EXPECT_EQ(x, ref);
}
TEST(MSBLAS, abs_sum_x_1)
{
	std::vector<double> x = { -1,-2,3 };

	const auto n = static_cast<int>(x.size());
	const auto result = ms::math::blas::abs_sum_x(x.data(), n);

	const double ref = 6;
	EXPECT_EQ(result, ref);
}
TEST(MSBLAS, cx_1)
{
	std::vector<double> x = { 1,2,3 };

	constexpr auto c = 2;
	const auto n = static_cast<int>(x.size());
	ms::math::blas::cx(c, x.data(), n);

	const std::vector<double> ref = { 2,4,6 };
	EXPECT_EQ(x, ref);
}
TEST(MSBLAS, cx_2)
{
	std::vector<double> x = { 1,2,3 };

	constexpr auto c = 2;
	const auto n = static_cast<int>(x.size());

	std::vector<double> result(n);
	ms::math::blas::cx(result.data(), c, x.data(), n);

	const std::vector<double> ref = { 2,4,6 };
	EXPECT_EQ(result, ref);
}
TEST(MSBLAS, x_dot_y_1)
{
	const std::vector<double> x = { 1,2,3 };
	const std::vector<double> y = { 4,5,6 };

	const auto n = static_cast<int>(x.size());
	const auto result = ms::math::blas::x_dot_y(x.data(), y.data(), n);

	const double ref = 32;
	EXPECT_EQ(result, ref);
}
TEST(MSBLAS, x_dot_y_2)
{
	const std::vector<double> x = { 3,4 };

	const auto n = static_cast<int>(x.size());
	const auto result = ms::math::blas::x_dot_y(x.data(), x.data(), n);

	const auto ref = 25;
	EXPECT_EQ(result, ref);
}
TEST(MSBLAS, x_plus_assign_y_1)
{
	std::vector<double> x = { 1,2,3 };
	const std::vector<double> y = { 4,5,7 };

	const auto n = static_cast<int>(x.size());
	ms::math::blas::x_plus_assign_y(x.data(), y.data(), n);

	const std::vector<double> ref = { 5,7,10 };
	EXPECT_EQ(x, ref);
}
TEST(MSBLAS, x_minus_assign_y_1)
{
	std::vector<double> x = { 1,2,3 };
	const std::vector<double> y = { 4,5,7 };

	const auto n = static_cast<int>(x.size());
	ms::math::blas::x_minus_assign_y(x.data(), y.data(), n);

	const std::vector<double> ref = { -3,-3,-4 };
	EXPECT_EQ(x, ref);
}
TEST(MSBLAS, find_maximum_element_pos_1)
{
	std::vector<double> x = { 1,2,3,4,5,6 };

	constexpr auto n = 3;
	constexpr auto inc = 2;
	const auto result = ms::math::blas::find_maximum_element_pos(x.data(), n, inc);

	constexpr auto ref = 2;
	EXPECT_EQ(result, ref);
}
TEST(MSBLAS, mv_1)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 3;
	constexpr auto leading_dimension = 3;
	std::vector<double> m = { 1,2,3,4,5,6,7,8,9 };

	const std::vector<double> v = { -1,-1,1 };

	std::vector<double> result(num_row);
	ms::math::blas::Ax(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 0,-3,-6 };
	EXPECT_EQ(result, ref);
}
TEST(MSBLAS, mv_2)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	constexpr auto leading_dimension = 3;
	std::vector<double> m = { 1,2,3,4,5,6,7,8,9 };

	const std::vector<double> v = { 2,-1 };

	std::vector<double> result(num_row);
	ms::math::blas::Ax(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 0,3,6 };
	EXPECT_EQ(result, ref);
}
TEST(MSBLAS, mv_3)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	constexpr auto leading_dimension = 4;
	std::vector<double> m = { 1,2,3,4,5,6,7,8,9,10,11,12 };

	const std::vector<double> v = { -1,1 };

	std::vector<double> result(num_row);
	ms::math::blas::Ax(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 1,1,1 };
	EXPECT_EQ(result, ref);
}
TEST(MSBLAS, mtv_1)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 3;
	constexpr auto leading_dimension = 3;
	std::vector<double> m = { 1,2,3,4,5,6,7,8,9 };

	const std::vector<double> v = { -1,-1,1 };

	std::vector<double> result(num_row);
	ms::math::blas::ATx(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 2,1,0 };
	EXPECT_EQ(result, ref);
}
TEST(MSBLAS, mtv_2)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	constexpr auto leading_dimension = 2;
	std::vector<double> m = { 1,2,3,4,5,6 };

	const std::vector<double> v = { -1,-1,1 };

	std::vector<double> result(num_column);
	ms::math::blas::ATx(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 1,0 };
	EXPECT_EQ(result, ref);
}
TEST(MSBLAS, mtv_3)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	constexpr auto leading_dimension = 3;
	std::vector<double> m = { 1,2,3,4,5,6,7,8,9 };

	const std::vector<double> v = { -1,-1,1 };

	std::vector<double> result(num_column);
	ms::math::blas::ATx(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 2,1 };
	EXPECT_EQ(result, ref);
}
TEST(MSBLAS, mtv_4)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	constexpr auto leading_dimension = 4;
	std::vector<double> m = { 1,2,3,4,5,6,7,8,9,10,11,12 };

	const std::vector<double> v = { -1,-1,1 };

	std::vector<double> result(num_column);
	ms::math::blas::ATx(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 3,2 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas, A_plus_cBT_1)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 3;
	constexpr auto A_leading_dimension = 5;
	std::vector<double> A = { 1,2,3,9,9,4,5,6,9,9 };

	//constexpr auto B_num_rows = 3;
	//constexpr auto B_num_columns = 2;
	constexpr auto B_leading_dimension = 3;
	std::vector<double> B = { 1,2,0,3,4,0,5,6,0 };

	constexpr auto c = 1.0;
	ms::math::blas::manual::A_plus_cBT(A.data(), A.data(), c, B.data(), A_num_rows, A_num_columns, A_leading_dimension, A_leading_dimension, B_leading_dimension);

	const std::vector<double> ref = { 2,5,8,9,9,6,9,12,9,9 };
	EXPECT_EQ(A, ref);
}
TEST(ms_blas, cAx)
{
  constexpr auto      A_num_rows          = 2;
  constexpr auto      A_num_columns       = 3;
  constexpr auto      A_leading_dimension = 5;
  std::vector<double> A                   = {1, 2, 3, 9, 9, 4, 5, 6, 9, 9};
  std::vector<double> x                   = {1, 2, 1};

	std::vector<double> result(2);

  constexpr auto c = 0.5;
  ms::math::blas::cAx(result.data(), c, A.data(), x.data(), A_num_rows, A_num_columns, A_leading_dimension);

  const std::vector<double> ref = {4,10};
  EXPECT_EQ(result, ref);
}