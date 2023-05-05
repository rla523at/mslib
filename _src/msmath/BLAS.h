#pragma once

namespace ms::math::blas::manual
{
void   abs_x(double* result_ptr, const double* x_ptr, const int n, const int incr = 1, const int incx = 1);
void   abs_x(double* x_ptr, const int n, const int incx = 1);
double abs_sum_x(const double* x_ptr, const int n, const int incx = 1);
void   copy(double* result_ptr, const double* x_ptr, const int n, const int incr = 1, const int incx = 1);
void   cx(double* result_ptr, const double c, const double* x_ptr, const int n, const int incr = 1, const int incx = 1);
void   cx(const double c, double* x_ptr, const int n, const int incx = 1);
double x_dot_y(const double* x_ptr, const double* y_ptr, const int n, const int incx = 1, const int incy = 1);
void   x_plus_y(double* result_ptr, const double* x_ptr, const double* y_ptr, const int n, const int incr = 1, const int incx = 1, const int incy = 1);
void   x_plus_assign_y(double* x_ptr, const double* y_ptr, const int n, const int incx = 1, const int incy = 1);
void   x_plus_assign_cy(double* x_ptr, const double c, const double* y_ptr, const int n, const int incx = 1, const int incy = 1);
void   x_minus_y(double* result_ptr, const double* x_ptr, const double* y_ptr, const int n, const int incr = 1, const int incx = 1, const int incy = 1);
void   x_minus_assign_y(double* x_ptr, const double* y_ptr, const int n, const int incx = 1, const int incy = 1);

// method for row major matrix
void cA(const double c, double* A_ptr, const int num_rows, const int num_columns, const int A_leading_dimension);
void cA(double* result_ptr, const double c, const double* A_ptr, const int num_rows, const int num_columns, const int A_leading_dimension);
void A_plus_cB(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int result_num_rows, const int result_num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension);
void AT_plus_cB(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int result_num_rows, const int result_num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension);
void A_plus_cBT(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int result_num_rows, const int result_num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension);
void AT_plus_cBT(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int result_num_rows, const int result_num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension);
void cAB(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_columns, const int A_leading_dimension, const int B_leading_dimension);
void cATB(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_columns, const int A_leading_dimension, const int B_leading_dimension);
void cABT(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_rows, const int A_leading_dimension, const int B_leading_dimension);
void cATBT(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_rows, const int A_leading_dimension, const int B_leading_dimension);

void Ax(double* result_ptr, const double* A_ptr, const double* x_ptr, const int num_rows, const int num_columns, const int leading_dimension, const int incr = 1, const int incv = 1);
void ATx(double* result_ptr, const double* A_ptr, const double* x_ptr, const int num_rows, const int num_columns, const int leading_dimension, const int incr = 1, const int incv = 1);
} // namespace ms::math::blas::manual

namespace ms::math::blas::mkl
{
double abs_sum_x(const double* x_ptr, const int n, const int incx = 1);
void   copy(double* result_ptr, const double* x_ptr, const int n, const int incr = 1, const int incx = 1);
void   cx(double* result_ptr, const double c, const double* x_ptr, const int n, const int incr = 1, const int incx = 1);
void   cx(const double c, double* x_ptr, const int n, const int incx = 1);
double x_dot_y(const double* x_ptr, const double* y_ptr, const int n, const int incx = 1, const int incy = 1);
void   x_plus_y(double* result_ptr, const double* x_ptr, const double* y_ptr, const int n, const int incr = 1, const int incx = 1, const int incy = 1);
void   x_plus_assign_y(double* x_ptr, const double* y_ptr, const int n, const int incx = 1, const int incy = 1);
void   x_plus_assign_cy(double* x_ptr, const double c, const double* y_ptr, const int n, const int incx = 1, const int incy = 1);
void   x_minus_y(double* result_ptr, const double* x_ptr, const double* y_ptr, const int n, const int incr = 1, const int incx = 1, const int incy = 1);
void   x_minus_assign_y(double* x_ptr, const double* y_ptr, const int n, const int incx = 1, const int incy = 1);

// method for row major matrix
void cA(const double c, double* A_ptr, const int num_rows, const int num_columns, const int A_leading_dimension);
void cA(double* result_ptr, const double c, const double* A_ptr, const int num_rows, const int num_columns, const int A_leading_dimension);
void A_plus_cB(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int result_num_rows, const int result_num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension);
void AT_plus_cB(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int result_num_rows, const int result_num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension);
void A_plus_cBT(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int result_num_rows, const int result_num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension);
void AT_plus_cBT(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int result_num_rows, const int result_num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension);
void cAB(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_columns, const int A_leading_dimension, const int B_leading_dimension);
void cATB(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_columns, const int A_leading_dimension, const int B_leading_dimension);
void cABT(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_rows, const int A_leading_dimension, const int B_leading_dimension);
void cATBT(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_rows, const int A_leading_dimension, const int B_leading_dimension);

void Ax(double* result_ptr, const double* A_ptr, const double* x_ptr, const int num_rows, const int num_columns, const int leading_dimension, const int incr = 1, const int incx = 1);
void ATx(double* result_ptr, const double* A_ptr, const double* x_ptr, const int num_rows, const int num_columns, const int leading_dimension, const int incr = 1, const int incx = 1);
void invA(double* A_ptr, const int num_rows, const int num_columns, const int leading_dimension);
} // namespace ms::math::blas::mkl

namespace ms::math::blas
{
inline constexpr int copy_criteria            = 50;
inline constexpr int scal_criteria            = 10;
inline constexpr int asum_criteria            = 10;
inline constexpr int axpy_criteria            = 20;
inline constexpr int dot_criteria             = 15;
inline constexpr int gemv_criteria            = 50;
inline constexpr int gemm_criteria            = 25;
inline constexpr int cA_criteria              = 25;
inline constexpr int A_plus_cB_criteria       = 25;
inline constexpr int AT_plus_cB_criteria      = 50;
inline constexpr int A_plus_cBT_criteria      = 50;
inline constexpr int AT_plus_cBT_criteria     = 65;
inline constexpr int cAB_A_rows_criteria      = 7;
inline constexpr int cAB_size_criteria        = 50;
inline constexpr int cABT_A_rows_criteria     = 7;
inline constexpr int cABT_size_criteria       = 50;
inline constexpr int cATB_A_columns_criteria  = 7;
inline constexpr int cATB_size_criteria       = 50;
inline constexpr int cATBT_A_columns_criteria = 7;
inline constexpr int cATBT_size_criteria      = 50;

void   abs_x(double* result_ptr, const double* x_ptr, const int n, const int incr = 1, const int incx = 1);
void   abs_x(double* x_ptr, const int n, const int incx = 1);
double abs_sum_x(const double* x_ptr, const int n, const int incx = 1);
void   copy(double* result_ptr, const double* x_ptr, const int n, const int incr = 1, const int incx = 1);
void   cx(double* result_ptr, const double c, const double* x_ptr, const int n, const int incr = 1, const int incx = 1);
void   cx(const double c, double* x_ptr, const int n, const int incx = 1);
size_t find_maximum_element_pos(const double* x_ptr, const int n, const int incx = 1);
double x_dot_y(const double* x_ptr, const double* y_ptr, const int n, const int incx = 1, const int incy = 1);
void   x_plus_y(double* result_ptr, const double* x_ptr, const double* y_ptr, const int n, const int incr = 1, const int incx = 1, const int incy = 1);
void   x_plus_assign_y(double* x_ptr, const double* y_ptr, const int n, const int incx = 1, const int incy = 1);
void   x_plus_assign_cy(double* x_ptr, const double c, const double* y_ptr, const int n, const int incx = 1, const int incy = 1);
void   x_minus_y(double* result_ptr, const double* x_ptr, const double* y_ptr, const int n, const int incr = 1, const int incx = 1, const int incy = 1);
void   x_minus_assign_y(double* x_ptr, const double* y_ptr, const int n, const int incx = 1, const int incy = 1);

// method for row major matrix
void cA(const double c, double* A_ptr, const int num_rows, const int num_columns, const int A_leading_dimension);
void cA(double* result_ptr, const double c, const double* A_ptr, const int num_rows, const int num_columns, const int A_leading_dimension);
void A_plus_cB(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int result_num_rows, const int result_num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension);
void AT_plus_cB(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int result_num_rows, const int result_num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension);
void A_plus_cBT(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int result_num_rows, const int result_num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension);
void AT_plus_cBT(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int result_num_rows, const int result_num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension);
void cAB(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_columns, const int A_leading_dimension, const int B_leading_dimension);
void cATB(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_columns, const int A_leading_dimension, const int B_leading_dimension);
void cABT(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_rows, const int A_leading_dimension, const int B_leading_dimension);
void cATBT(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_rows, const int A_leading_dimension, const int B_leading_dimension);

void Ax(double* result_ptr, const double* A_ptr, const double* x_ptr, const int num_rows, const int num_columns, const int leading_dimension, const int incr = 1, const int incv = 1);
void ATx(double* result_ptr, const double* A_ptr, const double* x_ptr, const int num_rows, const int num_columns, const int leading_dimension, const int incr = 1, const int incv = 1);
void invA(double* A_ptr, const int num_rows, const int num_columns, const int leading_dimension);
} // namespace ms::math::blas