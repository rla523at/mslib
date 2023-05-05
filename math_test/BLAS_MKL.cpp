#include "API_BLAS.h"

TEST(ms_blas_mkl, abs_sum_x_1)
{
	std::vector<double> x = { -1,-2,3 };

	const auto n = static_cast<int>(x.size());
	const auto result = ms::math::blas::mkl::abs_sum_x(x.data(), n);

	const double ref = 6;
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cx_1)
{
	std::vector<double> x = { 1,2,3 };

	constexpr auto c = 2;
	const auto n = static_cast<int>(x.size());
	ms::math::blas::mkl::cx(c, x.data(), n);

	const std::vector<double> ref = { 2,4,6 };
	EXPECT_EQ(x, ref);
}
TEST(ms_blas_mkl, cx_2)
{
	std::vector<double> x = { 1,2,3 };

	constexpr auto c = 2;
	const auto n = static_cast<int>(x.size());

	std::vector<double> result(n);
	ms::math::blas::mkl::cx(result.data(), c, x.data(), n);

	const std::vector<double> ref = { 2,4,6 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cx_3)
{
	std::vector<double> x = { 1,2,3,4,5,6,7,8,9 };

	constexpr auto n = 3;
	constexpr auto inc = 3;

	const auto start_ptr = x.data() + 1;
	
	constexpr auto c = 2;
	ms::math::blas::mkl::cx(c, start_ptr, n, inc);

	std::vector<double> ref = { 1,4,3,4,10,6,7,16,9 };
	EXPECT_EQ(x, ref);
}
TEST(ms_blas_mkl, x_dot_y_1)
{
	const std::vector<double> x = { 1,2,3 };
	const std::vector<double> y = { 4,5,6 };

	const auto n = static_cast<int>(x.size());
	const auto result = ms::math::blas::mkl::x_dot_y(x.data(), y.data(), n);

	const double ref = 32;
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, x_dot_y_2)
{
	const std::vector<double> x = { 3,4 };

	const auto n = static_cast<int>(x.size());
	const auto result = ms::math::blas::mkl::x_dot_y(x.data(), x.data(), n);

	const auto ref = 25;
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, x_plus_assign_y_1)
{
	std::vector<double> x = { 1,2,3 };
	const std::vector<double> y = { 4,5,7 };

	const auto n = static_cast<int>(x.size());
	ms::math::blas::mkl::x_plus_assign_y(x.data(), y.data(), n);

	const std::vector<double> ref = { 5,7,10 };
	EXPECT_EQ(x, ref);
}
TEST(ms_blas_mkl, x_plus_assign_y_2)
{
	std::vector<double> x = { 1,2,3,4,5,6,7,8,9 };
	
	constexpr auto n = 3;
	constexpr auto inc = 3;

	const auto start_ptr1 = x.data();
	const auto start_ptr2 = x.data() + 1;

	ms::math::blas::mkl::x_plus_assign_y(start_ptr1, start_ptr2, n, inc, inc);

	const std::vector<double> ref = { 3,2,3,9,5,6,15,8,9 };
	EXPECT_EQ(x, ref);
}
TEST(ms_blas_mkl, cA_1)
{
	std::vector<double> A = { 1,2,3,4,5,6 };
	constexpr auto num_rows = 2;
	constexpr auto num_columns = 3;
	constexpr auto A_leading_dimension = 3;

	constexpr auto c = 2.0;
	ms::math::blas::mkl::cA(c, A.data(), num_rows, num_columns, A_leading_dimension);

	std::vector<double> ref = { 2,4,6,8,10,12 };
	EXPECT_EQ(A, ref);
}
TEST(ms_blas_mkl, cA_2)
{
	std::vector<double> A = { 1,2,3,8,9,4,5,6,3,4 };
	constexpr auto num_rows = 2;
	constexpr auto num_columns = 3;
	constexpr auto A_leading_dimension = 5;

	constexpr auto c = -1.0;
	ms::math::blas::mkl::cA(c, A.data(), num_rows, num_columns, A_leading_dimension);

	std::vector<double> ref = { -1,-2,-3,8,9,-4,-5,-6,3,4 };
	EXPECT_EQ(A, ref);
}
TEST(ms_blas_mkl, cA_3)
{
	std::vector<double> A = { 1,2,3,4,5,6 };
	constexpr auto num_rows = 2;
	constexpr auto num_columns = 3;
	constexpr auto A_leading_dimension = 3;

	constexpr auto c = 2.0;

	std::vector<double> result(num_rows * num_columns);
	ms::math::blas::mkl::cA(result.data(), c, A.data(), num_rows, num_columns, A_leading_dimension);

	std::vector<double> ref = { 2,4,6,8,10,12 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cA_4)
{
	std::vector<double> A = { 1,2,3,8,9,4,5,6,3,4 };
	constexpr auto num_rows = 2;
	constexpr auto num_columns = 3;
	constexpr auto A_leading_dimension = 5;

	constexpr auto c = -1.0;

	std::vector<double> result(num_rows * num_columns);
	ms::math::blas::mkl::cA(result.data(), c, A.data(), num_rows, num_columns, A_leading_dimension);

	std::vector<double> ref = { -1,-2,-3,-4,-5,-6 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, A_plus_cB_1)
{
	std::vector<double> A = { 1,2,3,4,5,6 };
	std::vector<double> B = { 1,2,3,4,5,6 };

	constexpr auto num_rows = 2;
	constexpr auto num_columns = 3;
	constexpr auto A_leading_dimension = 3;
	constexpr auto B_leading_dimension = 3;

	constexpr auto c = 2.0;
	std::vector<double> result(num_rows * num_columns);
	ms::math::blas::mkl::A_plus_cB(result.data(), A.data(), c, B.data(), num_rows, num_columns, num_columns, A_leading_dimension, B_leading_dimension);

	const std::vector<double> ref = { 3,6,9,12,15,18 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, A_plus_cB_2)
{
	std::vector<double> A = { 1,2,3,7,4,5,6,9 };
	std::vector<double> B = { 1,2,3,8,9,4,5,6,3,4 };

	constexpr auto num_rows = 2;
	constexpr auto num_columns = 3;
	constexpr auto A_leading_dimension = 4;
	constexpr auto B_leading_dimension = 5;

	constexpr auto c = 2.0;
	std::vector<double> result(num_rows * num_columns);
	ms::math::blas::mkl::A_plus_cB(result.data(), A.data(), c, B.data(), num_rows, num_columns, num_columns, A_leading_dimension, B_leading_dimension);

	const std::vector<double> ref = { 3,6,9,12,15,18 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, AT_plus_cB_1)
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
	ms::math::blas::mkl::AT_plus_cB(result.data(), A.data(), c, B.data(), B_num_rows, B_num_columns, B_num_columns, A_leading_dimension, B_leading_dimension);

	const std::vector<double> ref = { 2,6,5,9,8,12 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, AT_plus_cB_2)
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
	ms::math::blas::mkl::AT_plus_cB(result.data(), A.data(), c, B.data(), B_num_rows, B_num_columns, B_num_columns, A_leading_dimension, B_leading_dimension);

	const std::vector<double> ref = { 2,6,5,9,8,12 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, AT_plus_cB_3)
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
	ms::math::blas::mkl::AT_plus_cB(result.data(), A.data(), c, B.data(), B_num_rows, B_num_columns, B_num_columns, A_leading_dimension, B_leading_dimension);

	const std::vector<double> ref = { 2,5,8,11,7,10,13,16 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, A_plus_cBT_1)
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
	ms::math::blas::mkl::A_plus_cBT(result.data(), A.data(), c, B.data(), A_num_rows, A_num_columns, A_num_columns, A_leading_dimension, B_leading_dimension);

	const std::vector<double> ref = { 2,5,8,6,9,12 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, A_plus_cBT_2)
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
	ms::math::blas::mkl::A_plus_cBT(result.data(), A.data(), c, B.data(), A_num_rows, A_num_columns, A_num_columns, A_leading_dimension, B_leading_dimension);

	const std::vector<double> ref = { 2,5,8,6,9,12 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, AT_plus_cBT_1)
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
	ms::math::blas::mkl::AT_plus_cBT(result.data(), A.data(), c, B.data(), result_num_rows, result_num_columns, result_leading_dimension, A_leading_dimension, B_leading_dimension);

	const std::vector<double> ref = { 2,8,4,10,6,12 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, AT_plus_cBT_2)
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
	ms::math::blas::mkl::AT_plus_cBT(result.data(), A.data(), c, B.data(), result_num_rows, result_num_columns, result_leading_dimension, A_leading_dimension, B_leading_dimension);

	const std::vector<double> ref = { 3,12,6,15,9,18 };
	EXPECT_EQ(result, ref);
}

//TEST(Matrix, operator_addition1)
//{
//	const Matrix m1 = { 1,2, { 1,2 } };
//	const Matrix m2 = { 1,2, { 2,2 } };
//	const auto result = m1 + m2;
//
//	const Matrix ref = { 1,2,{ 3,4 } };
//	EXPECT_EQ(result, ref);
//}
//TEST(Matrix, operator_addition2)
//{
//	const Matrix m1 = { 2,1, { 1,2 } };
//	const Matrix m2 = { 2,1, { 2,2 } };
//	const auto result = m1 + m2;
//
//	const Matrix ref = { 2,1,{ 3,4 } };
//	EXPECT_EQ(result, ref);
//}
//TEST(Matrix, operator_addition3)
//{
//	Matrix m1(2, 2, { 1,2,3,4 });
//	Matrix m2(2, 2, { 1,3,4,5 });
//	const auto result = m1 + m2;
//
//	Matrix ref(2, 2, { 2,5,7,9 });
//	EXPECT_EQ(result, ref);
//}
//TEST(Matrix, operator_addition4) 
//{
//	Matrix m1(2, 3, { 1,2,3,4,5,6 });
//	Matrix m2(3, 2, { 1,3,4,5,5,6 });
//	EXPECT_ANY_THROW(m1 + m2);
//}
//TEST(Matrix, operator_addition5) 
//{
//	Matrix m1(2, 3, { 1,2,3,4,5,6 });
//	Matrix m2(3, 2, { 1,3,4,5,5,6 });
//	m2.transpose();
//	EXPECT_ANY_THROW(m1 + m2);
//
//	//const auto result = m1 + m2;
//
//	//Matrix ref(2, 3, { 2,6,8,7,10,12 });
//	//EXPECT_EQ(result, ref);
//}
//TEST(Matrix, operator_addition6) 
//{
//	Matrix m1(2, 3, { 1,2,3,4,5,6 });
//	Matrix m2(3, 2, { 1,3,4,5,5,6 });
//	m2.transpose();
//	EXPECT_ANY_THROW(m2 + m1);
//
//	//const auto result = m2 + m1;
//
//	//Matrix ref(2, 3, { 2,6,8,7,10,12 });
//	//EXPECT_EQ(result, ref);
//}
//TEST(Matrix, operator_addition7) 
//{
//	Matrix m1(2, 2, { 1,2,3,4 });
//	Matrix m2(2, 2, { 1,3,4,5 });
//	m1.transpose();
//	m2.transpose();
//	EXPECT_ANY_THROW(m1 + m2);
//
//	//const auto result = m1 + m2;
//
//	//Matrix ref(2, 2, { 2,7,5,9 });
//	//EXPECT_EQ(result, ref);
//}
//TEST(Matrix, operator_addition8) 
//{
//	Matrix m1(2, 5, { 1.2345, 2.346345, 6.3262345, 8.5674567, 6.23452345, 2.53462, 6.432452345, 2.345345, 1.234563245, 7.3245345 });
//	Matrix m2(5, 2, { 1.234234, 2.3462345, 345.324, 2.6345345, 634523.5, 2345345.3,	 23453.345, 234534.6,	 234523.5, 623452.1 });
//	m2.transpose();
//	EXPECT_ANY_THROW(m1 + m2);
//	//const auto result = m1 + m2;
//
//	//Matrix ref(2, 5, { 2.468734,  347.670345, 634529.8262345,    23461.9124567, 234529.73452345,4.8808545, 9.066986845, 2345347.645345, 234535.834563245,  623459.4245345 });
//	//EXPECT_EQ(result, ref);
//}
//TEST(Matrix, operator_addition9) 
//{
//	Matrix m1(3, 5, { 1.2345, 2.346345, 6.3262345, 8.5674567, 6.23452345, 2.53462, 6.432452345, 2.345345, 1.234563245, 7.3245345,  789.45978, 74.5789123, 74.23541, 4.7894113, 7894.5134 });
//	Matrix m2(5, 3, { 1.234234, 2.3462345, 789456.0,   345.324, 2.6345345, 74.48651,  634523.5, 2345345.3, 710.1846, 23453.345,  234534.6,  12.5487,  234523.5,  623452.1, 421.7456 });
//	m2.transpose();
//	EXPECT_ANY_THROW(m1 + m2);
//	//const auto result = m1 + m2;
//
//	//Matrix ref(3, 5, { 2.468734,  347.670345, 634529.8262345,    23461.9124567, 234529.73452345,    4.8808545, 9.066986845, 2345347.645345, 234535.834563245,  623459.4245345, 790245.45978, 149.0654223,      784.42001,       17.3381113,        8316.259 });
//	//EXPECT_EQ(result, ref);
//}

TEST(ms_blas_mkl, mv_1)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 3;
	constexpr auto leading_dimension = 3;
	std::vector<double> m = { 1,2,3,4,5,6,7,8,9 };

	const std::vector<double> v = { -1,-1,1 };

	std::vector<double> result(num_row);
	ms::math::blas::mkl::Ax(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 0,-3,-6 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, mv_2)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	constexpr auto leading_dimension = 3;
	std::vector<double> m = { 1,2,3,4,5,6,7,8,9 };

	const std::vector<double> v = { 2,-1 };

	std::vector<double> result(num_row);
	ms::math::blas::mkl::Ax(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 0,3,6 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, mv_3)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	constexpr auto leading_dimension = 4;
	std::vector<double> m = { 1,2,3,4,5,6,7,8,9,10,11,12 };

	const std::vector<double> v = { -1,1 };

	std::vector<double> result(num_row);
	ms::math::blas::mkl::Ax(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 1,1,1 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, mtv_1)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 3;
	constexpr auto leading_dimension = 3;
	std::vector<double> m = { 1,2,3,4,5,6,7,8,9 };

	const std::vector<double> v = { -1,-1,1 };

	std::vector<double> result(num_row);
	ms::math::blas::mkl::ATx(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 2,1,0 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, mtv_2)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	constexpr auto leading_dimension = 2;
	std::vector<double> m = { 1,2,3,4,5,6 };

	const std::vector<double> v = { -1,-1,1 };

	std::vector<double> result(num_column);
	ms::math::blas::mkl::ATx(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 1,0 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, mtv_3)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	constexpr auto leading_dimension = 3;
	std::vector<double> m = { 1,2,3,4,5,6,7,8,9 };

	const std::vector<double> v = { -1,-1,1 };

	std::vector<double> result(num_column);
	ms::math::blas::mkl::ATx(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 2,1 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, mtv_4)
{
	constexpr auto num_row = 3;
	constexpr auto num_column = 2;
	constexpr auto leading_dimension = 4;
	std::vector<double> m = { 1,2,3,4,5,6,7,8,9,10,11,12 };

	const std::vector<double> v = { -1,-1,1 };

	std::vector<double> result(num_column);
	ms::math::blas::mkl::ATx(result.data(), m.data(), v.data(), num_row, num_column, leading_dimension);

	const std::vector<double> ref = { 3,2 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cAB_1)
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
	ms::math::blas::mkl::cAB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 7,10,15,22 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cAB_2)
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
	ms::math::blas::mkl::cAB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 204,204,120,120 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cAB_3)
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
	ms::math::blas::mkl::cAB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 204,204,120,120 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cAB_4)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 2;
	constexpr auto A_leading_dimeension = 2;
	std::vector<double> A_vals = { 1,2,3,4 };

	//constexpr auto B_num_rows = 2;
	constexpr auto B_num_columns = 2;
	constexpr auto B_leading_dimeension = 2;
	std::vector<double> B_vals = { 1,2,3,4 };

	constexpr double c = -1.0;
	std::vector<double> result(A_num_rows * B_num_columns);
	ms::math::blas::mkl::cAB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { -7,-10,-15,-22 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl , cATB_1)
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
	ms::math::blas::mkl::cATB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 13,18,18,26 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cATB_2)
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
	ms::math::blas::mkl::cATB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 13,18,18,26 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cATB_3)
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
	ms::math::blas::mkl::cATB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 13,10,7,17,14,11,21,18,15 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cATB_4)
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
	ms::math::blas::mkl::cATB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 18.470224531309999,  16.712738084116999,  10.418514118810000,  22.167911283120002,  21.030398669553001,  14.150019875622000,  27.174477241769999,  26.434492089979003,  18.454029295989997 };

	for (size_t i = 0; i < result.size(); ++i)
	{
		EXPECT_DOUBLE_EQ(result[i], ref[i]);
	}
}
TEST(ms_blas_mkl, cATB_5)
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
	ms::math::blas::mkl::cATB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 13,18,18,26 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cATB_6)
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
	ms::math::blas::mkl::cATB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 13,10,7,17,14,11,21,18,15 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cATB_7)
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
	ms::math::blas::mkl::cATB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 18.470224531309999,  16.712738084116999,  10.418514118810000,  22.167911283120002,  21.030398669553001,  14.150019875622000,  27.174477241769999,  26.434492089979003,  18.454029295989997 };

	for (size_t i = 0; i < result.size(); ++i)
	{
		EXPECT_DOUBLE_EQ(result[i], ref[i]);
	}
}
TEST(ms_blas_mkl, cATB_8)
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
	ms::math::blas::mkl::cATB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 18.470224531309999,  16.712738084116999,  10.418514118810000,  22.167911283120002,  21.030398669553001,  14.150019875622000,  27.174477241769999,  26.434492089979003,  18.454029295989997 };

	for (size_t i = 0; i < result.size(); ++i)
	{
		EXPECT_DOUBLE_EQ(result[i], ref[i]);
	}
}
TEST(ms_blas_mkl, cATB_9)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 2;
	constexpr auto A_leading_dimeension = 2;
	std::vector<double> A_vals = { 1,2,3,4 };

	//constexpr auto B_num_rows = 2;
	constexpr auto B_num_columns = 2;
	constexpr auto B_leading_dimeension = 2;
	std::vector<double> B_vals = { 1,3,4,5 };

	constexpr double c = -1.0;
	std::vector<double> result(A_num_columns * B_num_columns);
	ms::math::blas::mkl::cATB(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_columns, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { -13,-18,-18,-26 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cABT_1)
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
	ms::math::blas::mkl::cABT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 7,14,15,32 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cABT_2)
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
	ms::math::blas::mkl::cABT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 14,10,32,28 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cABT_3)
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
	ms::math::blas::mkl::cABT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 18.515714565372999,  17.342737125037932, 36.932522712599997 , 39.438937931479998 };
	blas_test_API::compare_considering_4ULP(result, ref);
}
TEST(ms_blas_mkl, cABT_4)
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
	ms::math::blas::mkl::cABT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 14,10,32,28 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cABT_5)
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
	ms::math::blas::mkl::cABT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 18.515714565372999,  17.342737125037932, 36.932522712599997 , 39.438937931479998 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cABT_6)
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
	ms::math::blas::mkl::cABT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 18.515714565372999,  17.342737125037932, 36.932522712599997 , 39.438937931479998 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cABT_7)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 2;
	constexpr auto A_leading_dimeension = 2;
	std::vector<double> A_vals = { 1,2,3,4 };

	constexpr auto B_num_rows = 2;
	//constexpr auto B_num_columns = 2;
	constexpr auto B_leading_dimeension = 2;
	std::vector<double> B_vals = { 1,3,4,5 };

	constexpr double c = -1.0;
	std::vector<double> result(A_num_rows * B_num_rows);
	ms::math::blas::mkl::cABT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { -7,-14,-15,-32 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cATBT_1)
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
	ms::math::blas::mkl::cATBT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 10,19,14,28 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cATBT_2)
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
	ms::math::blas::mkl::cATBT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 9,15,6,12,21,9,15,27,12 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cATBT_3)
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
	ms::math::blas::mkl::cATBT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { 9,15,6,12,21,9,15,27,12 };
	EXPECT_EQ(result, ref);
}
TEST(ms_blas_mkl, cATBT_4)
{
	constexpr auto A_num_rows = 2;
	constexpr auto A_num_columns = 2;
	constexpr auto A_leading_dimeension = 2;
	std::vector<double> A_vals = { 1,2,3,4 };

	constexpr auto B_num_rows = 2;
	//constexpr auto B_num_columns = 2;
	constexpr auto B_leading_dimeension = 2;
	std::vector<double> B_vals = { 1,3,4,5 };

	constexpr double c = -1.0;
	std::vector<double> result(A_num_columns * B_num_rows);
	ms::math::blas::mkl::cATBT(result.data(), c, A_vals.data(), B_vals.data(), A_num_rows, A_num_columns, B_num_rows, A_leading_dimeension, B_leading_dimeension);

	std::vector<double> ref = { -10,-19,-14,-28 };
	EXPECT_EQ(result, ref);
}