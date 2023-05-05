#include "API_BLAS.h"

TEST(ms_blas_manual, abs_x_1)
{
	std::vector<double> x = { -1,-2,3 };

	const auto n = static_cast<int>(x.size());
	ms::math::blas::manual::abs_x(x.data(), n);

	const std::vector<double> ref = { 1,2,3 };
	EXPECT_EQ(x, ref);
}
TEST(ms_blas_manual, abs_x_2)
{
	std::vector<double> x = { -1,-2,3,-4,-6 };

	constexpr auto n = 3;
	constexpr auto incx = 2;
	ms::math::blas::manual::abs_x(x.data(), n, incx);

	const std::vector<double> ref = { 1,-2,3,-4,6 };
	EXPECT_EQ(x, ref);
}
TEST(ms_blas_manual, abs_x_3)
{
	std::vector<double> x = { -1,-2,3,-4,-6 };

	constexpr auto n = 3;
	constexpr auto incx = 2;
	std::vector<double> result(n);
	constexpr auto incr = 1;
	ms::math::blas::manual::abs_x(result.data(),x.data(), n, incr, incx);

	const std::vector<double> ref = { 1,3,6 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, abs_sum_x_1)
{
	std::vector<double> x = { -1,-2,3 };

	const auto n = static_cast<int>(x.size());
	const auto result = ms::math::blas::manual::abs_sum_x(x.data(), n);

	const double ref = 6;
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, abs_sum_x_2)
{
	std::vector<double> x = { -1,-2,3 };

	constexpr auto n = 2;
	constexpr auto incx = 2;
	const auto result = ms::math::blas::manual::abs_sum_x(x.data(), n, incx);

	const double ref = 4;
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, cx_1)
{
	std::vector<double> x = { 1,2,3 };

	constexpr auto c = 2;
	const auto n = static_cast<int>(x.size());
	ms::math::blas::manual::cx(c, x.data(), n);

	const std::vector<double> ref = { 2,4,6 };
	EXPECT_EQ(x, ref);
}
TEST(ms_blas_manual, cx_2)
{
	std::vector<double> x = { 1,2,3 };

	constexpr auto c = 2;
	const auto n = static_cast<int>(x.size());
	
	std::vector<double> result(n);	
	ms::math::blas::manual::cx(result.data(), c, x.data(), n);

	const std::vector<double> ref = { 2,4,6 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, x_dot_y_1)
{
	const std::vector<double> x = { 1,2,3 };
	const std::vector<double> y = { 4,5,6 };

	const auto n = static_cast<int>(x.size());
	const auto result = ms::math::blas::manual::x_dot_y(x.data(), y.data(), n);

	const double ref = 32;
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, x_dot_y_2)
{
	const std::vector<double> x = { 3,4 };

	const auto n = static_cast<int>(x.size());
	const auto result = ms::math::blas::manual::x_dot_y(x.data(), x.data(), n);

	const auto ref = 25;
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, x_plus_assign_y_1)
{
	std::vector<double> x = { 1,2,3 };
	const std::vector<double> y = { 4,5,7 };

	const auto n = static_cast<int>(x.size());
	ms::math::blas::manual::x_plus_assign_y(x.data(), y.data(), n);

	const std::vector<double> ref = { 5,7,10 };
	EXPECT_EQ(x, ref);
}
TEST(ms_blas_manual, x_plus_assign_y_2)
{
	std::vector<double> x = { 1,2,3,4,5,6,7,8,9 };

	constexpr auto n = 3;
	constexpr auto inc = 3;

	const auto start_ptr1 = x.data();
	const auto start_ptr2 = x.data() + 1;

	ms::math::blas::manual::x_plus_assign_y(start_ptr1, start_ptr2, n, inc, inc);

	const std::vector<double> ref = { 3,2,3,9,5,6,15,8,9 };
	EXPECT_EQ(x, ref);
}
TEST(ms_blas_manual, x_minus_assign_y_1)
{
	std::vector<double> x = { 1,2,3 };
	const std::vector<double> y = { 4,5,7 };

	const auto n = static_cast<int>(x.size());
	ms::math::blas::manual::x_minus_assign_y(x.data(), y.data(), n);

	const std::vector<double> ref = { -3,-3,-4 };
	EXPECT_EQ(x, ref);
}
TEST(ms_blas_manual, cA_1)
{
	std::vector<double> A = { 1,2,3,4,5,6 };
	constexpr auto num_rows = 2;
	constexpr auto num_columns = 3;
	constexpr auto A_leading_dimension = 3;

	constexpr auto c = 2.0;
	ms::math::blas::manual::cA(c, A.data(), num_rows, num_columns, A_leading_dimension);

	std::vector<double> ref = { 2,4,6,8,10,12 };
	EXPECT_EQ(A, ref);
}
TEST(ms_blas_manual, cA_2)
{
	std::vector<double> A = { 1,2,3,8,9,4,5,6,3,4 };
	constexpr auto num_rows = 2;
	constexpr auto num_columns = 3;
	constexpr auto A_leading_dimension = 5;

	constexpr auto c = -1.0;
	ms::math::blas::manual::cA(c, A.data(), num_rows, num_columns, A_leading_dimension);

	std::vector<double> ref = { -1,-2,-3,8,9,-4,-5,-6,3,4 };
	EXPECT_EQ(A, ref);
}
TEST(ms_blas_manual, cA_3)
{
	std::vector<double> A = { 1,2,3,4,5,6 };
	constexpr auto num_rows = 2;
	constexpr auto num_columns = 3;
	constexpr auto A_leading_dimension = 3;

	constexpr auto c = 2.0;
	
	std::vector<double> result(num_rows * num_columns);
	ms::math::blas::manual::cA(result.data(), c, A.data(), num_rows, num_columns, A_leading_dimension);

	std::vector<double> ref = { 2,4,6,8,10,12 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, cA_4)
{
	std::vector<double> A = { 1,2,3,8,9,4,5,6,3,4 };
	constexpr auto num_rows = 2;
	constexpr auto num_columns = 3;
	constexpr auto A_leading_dimension = 5;

	constexpr auto c = -1.0;

	std::vector<double> result(num_rows * num_columns);
	ms::math::blas::manual::cA(result.data(), c, A.data(), num_rows, num_columns, A_leading_dimension);

	std::vector<double> ref = { -1,-2,-3,-4,-5,-6 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, A_plus_cB_1)
{
	std::vector<double> A = { 1,2,3,4,5,6 };
	std::vector<double> B = { 1,2,3,4,5,6 };

	constexpr auto num_rows = 2;
	constexpr auto num_columns = 3;
	constexpr auto A_leading_dimension = 3;
	constexpr auto B_leading_dimension = 3;

	constexpr auto c = 2.0;
	std::vector<double> result(num_rows * num_columns);
	ms::math::blas::manual::A_plus_cB(result.data(), A.data(), c, B.data(), num_rows, num_columns, num_columns, A_leading_dimension, B_leading_dimension);

	const std::vector<double> ref = { 3,6,9,12,15,18 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, A_plus_cB_2)
{
	std::vector<double> A = { 1,2,3,7,4,5,6,9 };
	std::vector<double> B = { 1,2,3,8,9,4,5,6,3,4 };

	constexpr auto num_rows = 2;
	constexpr auto num_columns = 3;
	constexpr auto A_leading_dimension = 4;
	constexpr auto B_leading_dimension = 5;

	constexpr auto c = 2.0;
	std::vector<double> result(num_rows * num_columns);
	ms::math::blas::manual::A_plus_cB(result.data(), A.data(), c, B.data(), num_rows, num_columns, num_columns, A_leading_dimension, B_leading_dimension);

	const std::vector<double> ref = { 3,6,9,12,15,18 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, AT_plus_cB_1)
{
	std::vector<double> A = { 1,2,3,4,5,6 };
	std::vector<double> B = { 1,2,3,4,5,6 };

	//constexpr auto A_num_rows = 2;
	//constexpr auto A_num_columns = 3;
	constexpr auto A_leading_dimension = 3;
	constexpr auto B_num_rows = 3;
	constexpr auto B_num_columns = 2;
	constexpr auto B_leading_dimension = 2;

	constexpr auto c = 1.0;
	std::vector<double> result(B_num_rows * B_num_columns);
	ms::math::blas::manual::AT_plus_cB(result.data(), A.data(), c, B.data(), B_num_rows, B_num_columns, B_num_columns, A_leading_dimension, B_leading_dimension);

	const std::vector<double> ref = { 2,6,5,9,8,12 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, AT_plus_cB_2)
{
	std::vector<double> A = { 1,2,3,9,4,5,6,9 };
	std::vector<double> B = { 1,2,9,3,4,9,5,6,9 };

	//constexpr auto A_num_rows = 2;
	//constexpr auto A_num_columns = 3;
	constexpr auto A_leading_dimension = 4;
	constexpr auto B_num_rows = 3;
	constexpr auto B_num_columns = 2;
	constexpr auto B_leading_dimension = 3;

	constexpr auto c = 1.0;
	std::vector<double> result(B_num_rows * B_num_columns);
	ms::math::blas::manual::AT_plus_cB(result.data(), A.data(), c, B.data(), B_num_rows, B_num_columns, B_num_columns, A_leading_dimension, B_leading_dimension);

	const std::vector<double> ref = { 2,6,5,9,8,12 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, AT_plus_cB_3)
{
	std::vector<double> A = { 1,2,3,4,5,6,7,8 };
	std::vector<double> B = { 1,2,3,4,9,5,6,7,8,9 };

	//constexpr auto A_num_rows = 4;
	//constexpr auto A_num_columns = 2;
	constexpr auto A_leading_dimension = 2;
	constexpr auto B_num_rows = 2;
	constexpr auto B_num_columns = 4;
	constexpr auto B_leading_dimension = 5;

	constexpr auto c = 1.0;
	std::vector<double> result(B_num_rows * B_num_columns);
	ms::math::blas::manual::AT_plus_cB(result.data(), A.data(), c, B.data(), B_num_rows, B_num_columns, B_num_columns, A_leading_dimension, B_leading_dimension);

	const std::vector<double> ref = { 2,5,8,11,7,10,13,16 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, A_plus_cBT_1)
{
	std::vector<double> A = { 1,2,3,4,5,6 };
	std::vector<double> B = { 1,2,3,4,5,6 };

	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 3;
	constexpr auto A_leading_dimension = 3;
	//constexpr auto B_num_rows = 3;
	//constexpr auto B_num_columns = 2;
	constexpr auto B_leading_dimension = 2;

	constexpr auto c = 1.0;
	std::vector<double> result(A_num_rows * A_num_columns);
	ms::math::blas::manual::A_plus_cBT(result.data(), A.data(), c, B.data(), A_num_rows, A_num_columns, A_num_columns, A_leading_dimension, B_leading_dimension);

	const std::vector<double> ref = { 2,5,8,6,9,12 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, A_plus_cBT_2)
{
	std::vector<double> A = { 1,2,3,9,9,4,5,6,9,9 };
	std::vector<double> B = { 1,2,0,3,4,0,5,6,0 };

	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 3;
	constexpr auto A_leading_dimension = 5;
	//constexpr auto B_num_rows = 3;
	//constexpr auto B_num_columns = 2;
	constexpr auto B_leading_dimension = 3;

	constexpr auto c = 1.0;
	std::vector<double> result(A_num_rows * A_num_columns);
	ms::math::blas::manual::A_plus_cBT(result.data(), A.data(), c, B.data(), A_num_rows, A_num_columns, A_num_columns, A_leading_dimension, B_leading_dimension);

	const std::vector<double> ref = { 2,5,8,6,9,12 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, AT_plus_cBT_1)
{
	constexpr auto result_num_rows = 3;
	constexpr auto result_num_columns = 2;
	constexpr auto result_leading_dimension = 2;

	//constexpr auto A_num_rows = result_num_columns;
	//constexpr auto A_num_columns = result_num_rows;
	constexpr auto A_leading_dimension = 3;
	std::vector<double> A = { 1,2,3,4,5,6 };

	//constexpr auto B_num_rows = result_num_columns;
	//constexpr auto B_num_columns = result_num_rows;
	constexpr auto B_leading_dimension = 3;
	std::vector<double> B = { 1,2,3,4,5,6 };

	constexpr auto c = 1.0;	
	std::vector<double> result(result_num_rows * result_num_columns);
	ms::math::blas::manual::AT_plus_cBT(result.data(), A.data(), c, B.data(), result_num_rows, result_num_columns, result_leading_dimension, A_leading_dimension, B_leading_dimension);

	const std::vector<double> ref = { 2,8,4,10,6,12 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, AT_plus_cBT_2)
{
	constexpr auto result_num_rows = 3;
	constexpr auto result_num_columns = 2;
	constexpr auto result_leading_dimension = 2;

	constexpr auto A_num_rows = result_num_columns;
	constexpr auto A_num_columns = result_num_rows;
	constexpr auto A_leading_dimension = 5;
	std::vector<double> A = { 1,2,3,0,0,4,5,6,0,0 };

	//constexpr auto B_num_rows = result_num_columns;
	//constexpr auto B_num_columns = result_num_rows;
	constexpr auto B_leading_dimension = 4;
	std::vector<double> B = { 1,2,3,0,4,5,6,0 };

	constexpr auto c = 2.0;
	std::vector<double> result(A_num_rows * A_num_columns);
	ms::math::blas::manual::AT_plus_cBT(result.data(), A.data(), c, B.data(), result_num_rows, result_num_columns, result_leading_dimension, A_leading_dimension, B_leading_dimension);

	const std::vector<double> ref = { 3,12,6,15,9,18 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, mv_1)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 3;
	constexpr auto leading_dimension = 3;
	std::vector<double> m = { 1,2,3,4,5,6,7,8,9 };

	const std::vector<double> v = { -1,-1,1 };

	std::vector<double> result(num_row);
	ms::math::blas::manual::Ax(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 0,-3,-6 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, mv_2)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	constexpr auto leading_dimension = 3;
	std::vector<double> m = { 1,2,3,4,5,6,7,8,9 };

	const std::vector<double> v = { 2,-1 };

	std::vector<double> result(num_row);
	ms::math::blas::manual::Ax(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 0,3,6 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, mv_3)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	constexpr auto leading_dimension = 4;
	std::vector<double> m = { 1,2,3,4,5,6,7,8,9,10,11,12 };

	const std::vector<double> v = { -1,1 };

	std::vector<double> result(num_row);
	ms::math::blas::manual::Ax(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 1,1,1 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, mtv_1)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 3;
	constexpr auto leading_dimension = 3;
	std::vector<double> m = { 1,2,3,4,5,6,7,8,9 };
	const std::vector<double> v = { -1,-1,1 };

	std::vector<double> result(num_row);
	ms::math::blas::manual::ATx(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 2,1,0 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, mtv_2)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	constexpr auto leading_dimension = 2;
	std::vector<double> m = { 1,2,3,4,5,6 };
	const std::vector<double> v = { -1,-1,1 };

	std::vector<double> result(num_column);
	ms::math::blas::manual::ATx(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 1,0 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, mtv_3)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	constexpr auto leading_dimension = 3;
	std::vector<double> m = { 1,2,3,4,5,6,7,8,9 };
	const std::vector<double> v = { -1,-1,1 };

	std::vector<double> result(num_column);
	ms::math::blas::manual::ATx(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 2,1 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, mtv_4)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	constexpr auto leading_dimension = 4;
	std::vector<double> m = { 1,2,3,4,5,6,7,8,9,10,11,12 };
	const std::vector<double> v = { -1,-1,1 };

	std::vector<double> result(num_column);
	ms::math::blas::manual::ATx(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 3,2 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, cAB_1)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 2;
	constexpr auto A_leading_dimeension = 2;
	std::vector<double> A_vals = { 1,2,3,4 };

	//constexpr auto B_num_rows = 2;
	constexpr auto B_num_columns = 2;
	constexpr auto B_leading_dimeension = 2;
	std::vector<double> B_vals = { 1,2,3,4 };

	constexpr double c = 1.0;
	std::vector<double> result(A_num_rows * B_num_columns);	
	ms::math::blas::manual::cAB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 7,10,15,22 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, cAB_2)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 8;
	constexpr auto A_leading_dimeension = 8;
	std::vector<double> A_vals = { 1,2,3,4,5,6,7,8,8,7,6,5,4,3,2,1 };

	//constexpr auto B_num_rows = 8;
	constexpr auto B_num_columns = 2;
	constexpr auto B_leading_dimeension = 2;
	std::vector<double> B_vals = { 1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8 };

	constexpr double c = 1.0;
	std::vector<double> result(A_num_rows * B_num_columns);
	ms::math::blas::manual::cAB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 204,204,120,120 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, cAB_3)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 8;
	constexpr auto A_leading_dimeension = 10;
	std::vector<double> A_vals = { 1,2,3,4,5,6,7,8,0,0,8,7,6,5,4,3,2,1,0,0 };

	//constexpr auto B_num_rows = 8;
	constexpr auto B_num_columns = 2;
	constexpr auto B_leading_dimeension = 3;
	std::vector<double> B_vals = { 1,1,0,2,2,0,3,3,0,4,4,0,5,5,0,6,6,0,7,7,0,8,8,0 };

	constexpr double c = 1.0;
	std::vector<double> result(A_num_rows * B_num_columns);
	ms::math::blas::manual::cAB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 204,204,120,120 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, cATB_1)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 2;
	constexpr auto A_leading_dimeension = 2;
	std::vector<double> A_vals = { 1,2,3,4 };

	//constexpr auto B_num_rows = 2;
	constexpr auto B_num_columns = 2;
	constexpr auto B_leading_dimeension = 2;
	std::vector<double> B_vals = { 1,3,4,5 };

	constexpr double c = 1.0;
	std::vector<double> result(A_num_columns * B_num_columns);
	ms::math::blas::manual::cATB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 13,18,18,26 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, cATB_2)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 2;
	constexpr auto A_leading_dimeension = 2;
	std::vector<double> A_vals = { 1,3,4,5 };

	//constexpr auto B_num_rows = 2;
	constexpr auto B_num_columns = 2;
	constexpr auto B_leading_dimeension = 2;
	std::vector<double> B_vals = { 1,2,3,4 };

	constexpr double c = 1.0;
	std::vector<double> result(A_num_columns * B_num_columns);
	ms::math::blas::manual::cATB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 13,18,18,26 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, cATB_3)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 3;
	constexpr auto A_leading_dimeension = 3;
	std::vector<double> A_vals = { 1,2,3,4,5,6 };

	//constexpr auto B_num_rows = 2;
	constexpr auto B_num_columns = 3;
	constexpr auto B_leading_dimeension = 3;
	std::vector<double> B_vals = { 1,2,3,3,2,1 };

	constexpr double c = 1.0;
	std::vector<double> result(A_num_columns * B_num_columns);
	ms::math::blas::manual::cATB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 13,10,7,17,14,11,21,18,15 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, cATB_4)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 3;
	constexpr auto A_leading_dimeension = 3;
	std::vector<double> A_vals = { 1.5479,2.4567123,3.414878,4.41487,5.121,6.15789 };

	//constexpr auto B_num_rows = 2;
	constexpr auto B_num_columns = 3;
	constexpr auto B_leading_dimeension = 3;
	std::vector<double> B_vals = { 1.1244,2.48711,3.12314,3.789413,2.9135491,1.264863 };

	constexpr double c = 1.0;
	std::vector<double> result(A_num_columns * B_num_columns);
	ms::math::blas::manual::cATB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 18.470224531309999,  16.712738084116999,  10.418514118810000,  22.167911283120002,  21.030398669553001,  14.150019875622000,  27.174477241769999,  26.434492089979003,  18.454029295989997 };
	blas_test_API::compare_considering_4ULP(result, ref);
}
TEST(ms_blas_manual, cATB_8)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 3;
	constexpr auto A_leading_dimeension = 5;
	std::vector<double> A_vals = { 1.5479,2.4567123,3.414878,0,0,4.41487,5.121,6.15789 };

	//constexpr auto B_num_rows = 2;
	constexpr auto B_num_columns = 3;
	constexpr auto B_leading_dimeension = 4;
	std::vector<double> B_vals = { 1.1244,2.48711,3.12314,0,3.789413,2.9135491,1.264863 };

	constexpr double c = 1.0;
	std::vector<double> result(A_num_columns * B_num_columns);
	ms::math::blas::manual::cATB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 18.470224531309999,  16.712738084116999,  10.418514118810000,  22.167911283120002,  21.030398669553001,  14.150019875622000,  27.174477241769999,  26.434492089979003,  18.454029295989997 };
	blas_test_API::compare_considering_4ULP(result, ref);
}
TEST(ms_blas_manual, cABT_1)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 2;
	constexpr auto A_leading_dimeension = 2;
	std::vector<double> A_vals = { 1,2,3,4 };

	constexpr auto B_num_rows = 2;
	//constexpr auto B_num_columns = 2;
	constexpr auto B_leading_dimeension = 2;
	std::vector<double> B_vals = { 1,3,4,5 };

	constexpr double c = 1.0;
	std::vector<double> result(A_num_rows * B_num_rows);
	ms::math::blas::manual::cABT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 7,14,15,32 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, cABT_2)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 3;
	constexpr auto A_leading_dimeension = 3;
	std::vector<double> A_vals = { 1,2,3,4,5,6 };

	constexpr auto B_num_rows = 2;
	//constexpr auto B_num_columns = 3;
	constexpr auto B_leading_dimeension = 3;
	std::vector<double> B_vals = { 1,2,3,3,2,1 };

	constexpr double c = 1.0;
	std::vector<double> result(A_num_rows * B_num_rows);
	ms::math::blas::manual::cABT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 14,10,32,28 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, cABT_3)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 3;
	constexpr auto A_leading_dimeension = 3;
	std::vector<double> A_vals = { 1.5479,2.4567123,3.414878,4.41487,5.121,6.15789 };

	constexpr auto B_num_rows = 2;
	//constexpr auto B_num_columns = 3;
	constexpr auto B_leading_dimeension = 3;
	std::vector<double> B_vals = { 1.1244,2.48711,3.12314,3.789413,2.9135491,1.264863 };

	constexpr double c = 1.0;
	std::vector<double> result(A_num_rows * B_num_rows);
	ms::math::blas::manual::cABT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 18.515714565372999,  17.342737125037932, 36.932522712599997 , 39.438937931479998 };
	blas_test_API::compare_considering_4ULP(result, ref);
}
TEST(ms_blas_manual, cABT_4)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 3;
	constexpr auto A_leading_dimeension = 3;
	std::vector<double> A_vals = { 1,2,3,4,5,6 };

	constexpr auto B_num_rows = 2;
	//constexpr auto B_num_columns = 3;
	constexpr auto B_leading_dimeension = 3;
	std::vector<double> B_vals = { 1,2,3,3,2,1 };

	constexpr double c = 1.0;
	std::vector<double> result(A_num_rows * B_num_rows);
	ms::math::blas::manual::cABT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 14,10,32,28 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, cABT_5)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 3;
	constexpr auto A_leading_dimeension = 3;
	std::vector<double> A_vals = { 1.5479,2.4567123,3.414878,4.41487,5.121,6.15789 };

	constexpr auto B_num_rows = 2;
	//constexpr auto B_num_columns = 3;
	constexpr auto B_leading_dimeension = 3;
	std::vector<double> B_vals = { 1.1244,2.48711,3.12314,3.789413,2.9135491,1.264863 };

	constexpr double c = 1.0;
	std::vector<double> result(A_num_rows * B_num_rows);
	ms::math::blas::manual::cABT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 18.515714565372999,  17.342737125037932, 36.932522712599997 , 39.438937931479998 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, cABT_6)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 3;
	constexpr auto A_leading_dimeension = 5;
	std::vector<double> A_vals = { 1.5479,2.4567123,3.414878,0,0,4.41487,5.121,6.15789 };

	constexpr auto B_num_rows = 2;
	//constexpr auto B_num_columns = 3;
	constexpr auto B_leading_dimeension = 5;
	std::vector<double> B_vals = { 1.1244,2.48711,3.12314,0,0,3.789413,2.9135491,1.264863 };

	constexpr double c = 1.0;
	std::vector<double> result(A_num_rows * B_num_rows);
	ms::math::blas::manual::cABT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 18.515714565372999,  17.342737125037932, 36.932522712599997 , 39.438937931479998 };
	EXPECT_EQ(result, ref);
}

TEST(ms_blas_manual, cATBT_1)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 2;
	constexpr auto A_leading_dimeension = 2;
	std::vector<double> A_vals = { 1,2,3,4 };

	constexpr auto B_num_rows = 2;
	//constexpr auto B_num_columns = 2;
	constexpr auto B_leading_dimeension = 2;
	std::vector<double> B_vals = { 1,3,4,5 };

	constexpr double c = 1.0;
	std::vector<double> result(A_num_columns * B_num_rows);
	ms::math::blas::manual::cATBT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 10,19,14,28 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, cATBT_2)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 3;
	constexpr auto A_leading_dimeension = 3;
	std::vector<double> A_vals = { 1,2,3,4,5,6 };

	constexpr auto B_num_rows = 3;
	//constexpr auto B_num_columns = 2;
	constexpr auto B_leading_dimeension = 2;
	std::vector<double> B_vals = { 1,2,3,3,2,1 };

	constexpr double c = 1.0;
	std::vector<double> result(A_num_columns * B_num_rows);
	ms::math::blas::manual::cATBT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 9,15,6,12,21,9,15,27,12 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_manual, cATBT_3)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 3;
	constexpr auto A_leading_dimeension = 4;
	std::vector<double> A_vals = { 1,2,3,0,4,5,6 };

	constexpr auto B_num_rows = 3;
	//constexpr auto B_num_columns = 2;
	constexpr auto B_leading_dimeension = 4;
	std::vector<double> B_vals = { 1,2,0,0,3,3,0,0,2,1 };

	constexpr double c = 1.0;
	std::vector<double> result(A_num_columns * B_num_rows);
	ms::math::blas::manual::cATBT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 9,15,6,12,21,9,15,27,12 };
	EXPECT_EQ(result, ref);
}


//TEST(Const_Matrix_Wrapper, operator_mm_5)
//{
//	Matrix m1(2, 2, { 1,2,3,4 });
//	Matrix m2(2, 2, { 1,3,4,5 });
//	m1.transpose();
//	m2.transpose();
//	const auto result = m1 * m2;
//
//	Matrix ref(2, 2, { 10,19,14,28 });
//	const auto [rows, cols] = ref.size();
//	for (size_t i = 0; i < rows; ++i)
//		for (size_t j = 0; j < cols; ++j)
//			EXPECT_DOUBLE_EQ(result.at(i, j), ref.at(i, j));
//}


//TEST(Const_Matrix_Wrapper, operator_mm_19)
//{
//	Matrix m1(2, 3, { 1.5479,2.4567123,3.414878,4.41487,5.121,6.15789 });
//	Matrix m2(2, 3, { 1.1244,2.48711,3.12314,3.789413,2.9135491,1.264863 });
//	m2.transpose();
//	const auto result = m1 * m2;
//
//	Matrix ref(2, 2, { 18.515714565372999,  17.342737125037932, 36.932522712599997 , 39.438937931479998 });
//	EXPECT_EQ(result, ref);
//}