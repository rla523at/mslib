#include "msmath/BLAS.h"

#include "mkl/mkl.h"
#include "msexception/Exception.h"
#include <cmath>
#include <vector>

namespace ms::math::blas::manual
{
void abs_x(double* result_ptr, const double* x_ptr, const int n, const int incr, const int incx)
{
  REQUIRE(0 < n, "n should be greater than 0");

  for (int i = 0; i < n; ++i)
  {
    result_ptr[i * incr] = std::abs(x_ptr[i * incx]);
  }
}

void abs_x(double* x_ptr, const int n, const int incx)
{
  REQUIRE(0 < n, "n should be greater than 0");

  for (int i = 0; i < n; ++i)
  {
    x_ptr[i * incx] = std::abs(x_ptr[i * incx]);
  }
}

double abs_sum_x(const double* x_ptr, const int n, const int incx)
{
  REQUIRE(0 < n, "n should be greater than 0");

  double result = 0.0;

  for (int i = 0; i < n; ++i)
  {
    result += std::abs(x_ptr[i * incx]);
  }

  return result;
}

void copy(double* result_ptr, const double* x_ptr, const int n, const int incr, const int incx)
{
  REQUIRE(0 < n, "n should be greater than 0");

  for (int i = 0; i < n; ++i)
  {
    result_ptr[i * incr] = x_ptr[i * incx];
  }
}

void cx(double* result_ptr, const double c, const double* x_ptr, const int n, const int incr, const int incx)
{
  REQUIRE(0 < n, "n should be greater than 0");

  for (int i = 0; i < n; ++i)
  {
    result_ptr[i * incr] = c * x_ptr[i * incx];
  }
}

void cx(const double c, double* x_ptr, const int n, const int incx)
{
  REQUIRE(0 < n, "n should be greater than 0");

  for (int i = 0; i < n; ++i)
  {
    x_ptr[i * incx] *= c;
  }
}

double x_dot_y(const double* x_ptr, const double* y_ptr, const int n, const int incx, const int incy)
{
  REQUIRE(0 < n, "n should be greater than 0");

  double result = 0.0;

  for (int i = 0; i < n; ++i)
  {
    result += x_ptr[i * incx] * y_ptr[i * incy];
  }

  return result;
}

void x_plus_y(double* result_ptr, const double* x_ptr, const double* y_ptr, const int n, const int incr, const int incx, const int incy)
{
  REQUIRE(0 < n, "n should be greater than 0");

  for (int i = 0; i < n; ++i)
  {
    result_ptr[i * incr] = x_ptr[i * incx] + y_ptr[i * incy];
  }
}

void x_plus_assign_y(double* x_ptr, const double* y_ptr, const int n, const int incx, const int incy)
{
  REQUIRE(0 < n, "n should be greater than 0");

  for (int i = 0; i < n; ++i)
  {
    x_ptr[i * incx] += y_ptr[i * incy];
  }
}

void x_plus_assign_cy(double* x_ptr, const double c, const double* y_ptr, const int n, const int incx, const int incy)
{
  REQUIRE(0 < n, "n should be greater than 0");

  for (int i = 0; i < n; ++i)
  {
    x_ptr[i * incx] += c * y_ptr[i * incy];
  }
}

void x_minus_y(double* result_ptr, const double* x_ptr, const double* y_ptr, const int n, const int incr, const int incx, const int incy)
{
  REQUIRE(0 < n, "n should be greater than 0");

  for (int i = 0; i < n; ++i)
  {
    result_ptr[i * incr] = x_ptr[i * incx] - y_ptr[i * incy];
  }
}

void x_minus_assign_y(double* x_ptr, const double* y_ptr, const int n, const int incx, const int incy)
{
  REQUIRE(0 < n, "n should be greater than 0");

  for (int i = 0; i < n; ++i)
  {
    x_ptr[i * incx] -= y_ptr[i * incy];
  }
}

void cA(const double c, double* A_ptr, const int num_rows, const int num_columns, const int leading_dimension)
{
  REQUIRE(0 < num_rows, "result_num_rows should be greater than 0");
  REQUIRE(0 < num_columns, "result_num_columns should be greater than 0");
  REQUIRE(num_columns <= leading_dimension, "leading_dimension error");

  for (int i = 0; i < num_rows; ++i)
  {
    const auto start_index = i * leading_dimension;
    for (int j = 0; j < num_columns; ++j)
    {
      A_ptr[start_index + j] *= c;
    }
  }
}

void cA(double* result_ptr, const double c, const double* A_ptr, const int num_rows, const int num_columns, const int A_leading_dimension)
{
  REQUIRE(0 < num_rows, "result_num_rows should be greater than 0");
  REQUIRE(0 < num_columns, "result_num_columns should be greater than 0");
  REQUIRE(num_columns <= A_leading_dimension, "leading_dimension error");

  for (int i = 0; i < num_rows; ++i)
  {
    const auto A_start_index      = i * A_leading_dimension;
    const auto result_start_index = i * num_columns;
    for (int j = 0; j < num_columns; ++j)
    {
      result_ptr[result_start_index + j] = c * A_ptr[A_start_index + j];
    }
  }
}

void A_plus_cB(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int result_num_rows, const int result_num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension)
{
  for (int i = 0; i < result_num_rows; ++i)
  {
    const auto A_start_index      = i * A_leading_dimension;
    const auto B_start_index      = i * B_leading_dimension;
    const auto result_start_index = i * result_leading_dimension;
    for (int j = 0; j < result_num_columns; ++j)
    {
      result_ptr[result_start_index + j] = A_ptr[A_start_index + j] + c * B_ptr[B_start_index + j];
    }
  }
}

void AT_plus_cB(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int result_num_rows, const int result_num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension)
{
  for (int i = 0; i < result_num_rows; ++i)
  {
    const auto result_start_index = i * result_leading_dimension;
    const auto B_start_index      = i * B_leading_dimension;
    for (int j = 0; j < result_num_columns; ++j)
    {
      result_ptr[result_start_index + j] = A_ptr[j * A_leading_dimension + i] + c * B_ptr[B_start_index + j];
    }
  }
}

void A_plus_cBT(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int result_num_rows, const int result_num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension)
{
  for (int i = 0; i < result_num_rows; ++i)
  {
    const auto A_start_index      = i * A_leading_dimension;
    const auto result_start_index = i * result_leading_dimension;
    for (int j = 0; j < result_num_columns; ++j)
    {
      result_ptr[result_start_index + j] = A_ptr[A_start_index + j] + c * B_ptr[j * B_leading_dimension + i];
    }
  }
}

void AT_plus_cBT(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int result_num_rows, const int result_num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension)
{
  for (int i = 0; i < result_num_rows; ++i)
  {
    const auto result_start_index = i * result_leading_dimension;
    for (int j = 0; j < result_num_columns; ++j)
    {
      result_ptr[result_start_index + j] = A_ptr[j * A_leading_dimension + i] + c * B_ptr[j * B_leading_dimension + i];
    }
  }
}

void cAB(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_columns, const int A_leading_dimension, const int B_leading_dimension)
{
  for (int i = 0; i < A_num_rows; ++i)
  {
    const auto A_start_index      = i * A_leading_dimension;
    const auto result_start_index = i * B_num_columns;
    for (int j = 0; j < B_num_columns; ++j)
    {
      for (int l = 0; l < A_num_columns; ++l)
      {
        result_ptr[result_start_index + j] += c * A_ptr[A_start_index + l] * B_ptr[l * B_leading_dimension + j];
      }
    }
  }
}

void cATB(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_columns, const int A_leading_dimension, const int B_leading_dimension)
{
  const auto AT_num_rows    = A_num_columns;
  const auto AT_num_columns = A_num_rows;

  for (int i = 0; i < AT_num_rows; ++i)
  {
    const auto result_start_index = i * B_num_columns;
    for (int j = 0; j < B_num_columns; ++j)
    {
      for (int l = 0; l < AT_num_columns; ++l)
      {
        result_ptr[result_start_index + j] += c * A_ptr[l * A_leading_dimension + i] * B_ptr[l * B_leading_dimension + j];
      }
    }
  }
}

void cABT(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_rows, const int A_leading_dimension, const int B_leading_dimension)
{
  const auto BT_num_columns = B_num_rows;

  for (int i = 0; i < A_num_rows; ++i)
  {
    const auto A_start_index      = i * A_leading_dimension;
    const auto result_start_index = i * BT_num_columns;
    for (int j = 0; j < BT_num_columns; ++j)
    {
      const auto B_start_index = j * B_leading_dimension;
      for (int l = 0; l < A_num_columns; ++l)
      {
        result_ptr[result_start_index + j] += c * A_ptr[A_start_index + l] * B_ptr[B_start_index + l];
      }
    }
  }
}

void cATBT(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_rows, const int A_leading_dimension, const int B_leading_dimension)
{
  const auto AT_num_rows    = A_num_columns;
  const auto AT_num_columns = A_num_rows;
  const auto BT_num_columns = B_num_rows;

  for (int i = 0; i < AT_num_rows; ++i)
  {
    const auto result_start_index = i * BT_num_columns;
    for (int j = 0; j < BT_num_columns; ++j)
    {
      const auto B_start_index = j * B_leading_dimension;
      for (int l = 0; l < AT_num_columns; ++l)
      {
        result_ptr[result_start_index + j] += c * A_ptr[l * A_leading_dimension + i] * B_ptr[B_start_index + l];
      }
    }
  }
}

void Ax(double* result_ptr, const double* A_ptr, const double* x_ptr, const int num_rows, const int num_columns, const int leading_dimension, const int incr, const int incx)
{
  REQUIRE(0 < num_rows, "result_num_rows should be greater than 0");
  REQUIRE(0 < num_columns, "result_num_columns should be greater than 0");
  REQUIRE(num_columns <= leading_dimension, "leading_dimension error");

  for (int i = 0; i < num_rows; ++i)
  {
    result_ptr[i * incr] = 0.0;
    for (int j = 0; j < num_columns; ++j)
    {
      result_ptr[i * incr] += A_ptr[i * leading_dimension + j] * x_ptr[j * incx];
    }
  }
}

void ATx(double* result_ptr, const double* A_ptr, const double* x_ptr, const int num_rows, const int num_columns, const int leading_dimension, const int incr, const int incx)
{
  REQUIRE(0 < num_rows, "result_num_rows should be greater than 0");
  REQUIRE(0 < num_columns, "result_num_columns should be greater than 0");
  REQUIRE(num_columns <= leading_dimension, "leading_dimension error");

  for (int i = 0; i < num_columns; ++i)
  {
    result_ptr[i * incr] = 0.0;
    for (int j = 0; j < num_rows; ++j)
    {
      result_ptr[i * incr] += A_ptr[j * leading_dimension + i] * x_ptr[j * incx];
    }
  }
}

void cAx(double* result_ptr, const double c, const double* A_ptr, const double* x_ptr, const int num_rows, const int num_columns, const int leading_dimension, const int incr, const int incx)
{
  REQUIRE(0 < num_rows, "result_num_rows should be greater than 0");
  REQUIRE(0 < num_columns, "result_num_columns should be greater than 0");
  REQUIRE(num_columns <= leading_dimension, "leading_dimension error");

  for (int i = 0; i < num_rows; ++i)
  {
    result_ptr[i * incr] = 0.0;
    for (int j = 0; j < num_columns; ++j)
    {
      result_ptr[i * incr] += A_ptr[i * leading_dimension + j] * x_ptr[j * incx];
    }

    result_ptr[i * incr] *= c;
  }
}

void cATx(double* result_ptr, const double c, const double* A_ptr, const double* x_ptr, const int num_rows, const int num_columns, const int leading_dimension, const int incr, const int incx)
{
  REQUIRE(0 < num_rows, "result_num_rows should be greater than 0");
  REQUIRE(0 < num_columns, "result_num_columns should be greater than 0");
  REQUIRE(num_columns <= leading_dimension, "leading_dimension error");

  for (int i = 0; i < num_columns; ++i)
  {
    result_ptr[i * incr] = 0.0;
    for (int j = 0; j < num_rows; ++j)
    {
      result_ptr[i * incr] += A_ptr[j * leading_dimension + i] * x_ptr[j * incx];
    }

    result_ptr[i * incr] *= c;
  }
}

} // namespace ms::math::blas::manual

namespace ms::math::blas::mkl
{
double abs_sum_x(const double* x_ptr, const int n, const int incx)
{
  REQUIRE(0 < n, "n should be greater than 0");

  return cblas_dasum(n, x_ptr, incx);
}

void copy(double* result_ptr, const double* x_ptr, const int n, const int incr, const int incx)
{
  REQUIRE(0 < n, "n should be greater than 0");

  cblas_dcopy(n, x_ptr, incx, result_ptr, incr);
}

void cx(double* result_ptr, const double c, const double* x_ptr, const int n, const int incr, const int incx)
{
  REQUIRE(0 < n, "n should be greater than 0");

  ms::math::blas::mkl::copy(result_ptr, x_ptr, n, incr, incx);
  ms::math::blas::mkl::cx(c, result_ptr, n, incr);
}

void cx(const double c, double* x_ptr, const int n, const int incx)
{
  REQUIRE(0 < n, "n should be greater than 0");

  cblas_dscal(n, c, x_ptr, incx);
}

double x_dot_y(const double* x_ptr, const double* y_ptr, const int n, const int incx, const int incy)
{
  REQUIRE(0 < n, "n should be greater than 0");

  return cblas_ddot(n, x_ptr, incx, y_ptr, incy);
}

void x_plus_y(double* result_ptr, const double* x_ptr, const double* y_ptr, const int n, const int incr, const int incx, const int incy)
{
  ms::math::blas::copy(result_ptr, x_ptr, n, incr, incx);
  ms::math::blas::x_plus_assign_y(result_ptr, y_ptr, n, incr, incy);
}

void x_plus_assign_y(double* x_ptr, const double* y_ptr, const int n, const int incx, const int incy)
{
  REQUIRE(0 < n, "n should be greater than 0");

  constexpr auto a = 1.0;
  cblas_daxpy(n, a, y_ptr, incx, x_ptr, incy);
}

void x_plus_assign_cy(double* x_ptr, const double c, const double* y_ptr, const int n, const int incx, const int incy)
{
  REQUIRE(0 < n, "n should be greater than 0");

  cblas_daxpy(n, c, y_ptr, incy, x_ptr, incx);
}

void x_minus_y(double* result_ptr, const double* x_ptr, const double* y_ptr, const int n, const int incr, const int incx, const int incy)
{
  ms::math::blas::copy(result_ptr, x_ptr, n, incr, incx);
  ms::math::blas::x_minus_assign_y(result_ptr, y_ptr, n, incr, incy);
}

void x_minus_assign_y(double* x_ptr, const double* y_ptr, const int n, const int incx, const int incy)
{
  REQUIRE(0 < n, "n should be greater than 0");

  constexpr auto a = -1.0;
  cblas_daxpy(n, a, y_ptr, incy, x_ptr, incx);
}

void cA(double* result_ptr, const double c, const double* A_ptr, const int num_rows, const int num_columns, const int A_leading_dimension)
{
  REQUIRE(0 < num_rows, "result_num_rows should be greater than 0");
  REQUIRE(0 < num_columns, "result_num_columns should be greater than 0");
  REQUIRE(num_columns <= A_leading_dimension, "leading_dimension error");

  if (num_columns == A_leading_dimension)
  {
    ms::math::blas::cx(result_ptr, c, A_ptr, num_rows * A_leading_dimension);
  }
  else
  {
    auto start_ptr = A_ptr;
    for (int i = 0; i < num_rows; i++)
    {
      ms::math::blas::cx(result_ptr, c, start_ptr, num_columns);
      result_ptr += num_columns;
      start_ptr += A_leading_dimension;
    }
  }
}

void cA(const double c, double* A_ptr, const int num_rows, const int num_columns, const int leading_dimension)
{
  REQUIRE(0 < num_rows, "result_num_rows should be greater than 0");
  REQUIRE(0 < num_columns, "result_num_columns should be greater than 0");
  REQUIRE(num_columns <= leading_dimension, "leading_dimension error");

  // MKL_Dimatcopy('R', 'N', num_rows, num_columns, c, A_ptr, leading_dimensions, leading_dimensions);

  if (num_columns == leading_dimension)
  {
    ms::math::blas::cx(c, A_ptr, num_rows * leading_dimension);
  }
  else
  {
    auto start_ptr = A_ptr;
    for (int i = 0; i < num_rows; i++)
    {
      ms::math::blas::cx(c, start_ptr, num_columns);
      start_ptr += leading_dimension;
    }
  }
}

void A_plus_cB(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int num_rows, const int num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension)
{
  REQUIRE(0 < num_rows, "result_num_rows should be greater than 0");
  REQUIRE(0 < num_columns, "result_num_columns should be greater than 0");
  REQUIRE(num_columns <= A_leading_dimension && num_columns <= B_leading_dimension, "leading_dimension error");

  MKL_Domatadd('R', 'N', 'N', num_rows, num_columns, 1.0, A_ptr, A_leading_dimension, c, B_ptr, B_leading_dimension, result_ptr, result_leading_dimension);
}

void AT_plus_cB(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int num_rows, const int num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension)
{
  REQUIRE(0 < num_rows, "result_num_rows should be greater than 0");
  REQUIRE(0 < num_columns, "result_num_columns should be greater than 0");
  REQUIRE(num_rows <= A_leading_dimension && num_columns <= B_leading_dimension, "leading_dimension error");

  MKL_Domatadd('R', 'T', 'N', num_rows, num_columns, 1.0, A_ptr, A_leading_dimension, c, B_ptr, B_leading_dimension, result_ptr, result_leading_dimension);
}

void A_plus_cBT(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int num_rows, const int num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension)
{
  REQUIRE(0 < num_rows, "result_num_rows should be greater than 0");
  REQUIRE(0 < num_columns, "result_num_columns should be greater than 0");
  REQUIRE(num_columns <= A_leading_dimension && num_rows <= B_leading_dimension, "leading_dimension error");

  MKL_Domatadd('R', 'N', 'T', num_rows, num_columns, 1.0, A_ptr, A_leading_dimension, c, B_ptr, B_leading_dimension, result_ptr, result_leading_dimension);
}

void AT_plus_cBT(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int num_rows, const int num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension)
{
  REQUIRE(0 < num_rows, "result_num_rows should be greater than 0");
  REQUIRE(0 < num_columns, "result_num_columns should be greater than 0");
  REQUIRE(num_rows <= A_leading_dimension && num_rows <= B_leading_dimension, "leading_dimension error");

  MKL_Domatadd('R', 'T', 'T', num_rows, num_columns, 1.0, A_ptr, A_leading_dimension, c, B_ptr, B_leading_dimension, result_ptr, result_leading_dimension);
}

void cAB(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_columns, const int A_leading_dimension, const int B_leading_dimension)
{
  const auto layout = CBLAS_LAYOUT::CblasRowMajor;
  const auto transA = CBLAS_TRANSPOSE::CblasNoTrans;
  const auto transB = CBLAS_TRANSPOSE::CblasNoTrans;

  cblas_dgemm(layout, transA, transB, A_num_rows, B_num_columns, A_num_columns, c, A_ptr, A_leading_dimension, B_ptr, B_leading_dimension, 0.0, result_ptr, B_num_columns);
}

void cATB(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_columns, const int A_leading_dimension, const int B_leading_dimension)
{
  const auto layout = CBLAS_LAYOUT::CblasRowMajor;
  const auto transA = CBLAS_TRANSPOSE::CblasTrans;
  const auto transB = CBLAS_TRANSPOSE::CblasNoTrans;

  const auto AT_num_rows    = A_num_columns;
  const auto AT_num_columns = A_num_rows;

  cblas_dgemm(layout, transA, transB, AT_num_rows, B_num_columns, AT_num_columns, c, A_ptr, A_leading_dimension, B_ptr, B_leading_dimension, 0.0, result_ptr, B_num_columns);
}
void cABT(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_rows, const int A_leading_dimension, const int B_leading_dimension)
{
  const auto layout = CBLAS_LAYOUT::CblasRowMajor;
  const auto transA = CBLAS_TRANSPOSE::CblasNoTrans;
  const auto transB = CBLAS_TRANSPOSE::CblasTrans;

  const auto BT_num_columns = B_num_rows;

  cblas_dgemm(layout, transA, transB, A_num_rows, BT_num_columns, A_num_columns, c, A_ptr, A_leading_dimension, B_ptr, B_leading_dimension, 0.0, result_ptr, BT_num_columns);
}
void cATBT(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_rows, const int A_leading_dimension, const int B_leading_dimension)
{
  const auto layout = CBLAS_LAYOUT::CblasRowMajor;
  const auto transA = CBLAS_TRANSPOSE::CblasTrans;
  const auto transB = CBLAS_TRANSPOSE::CblasTrans;

  const auto AT_num_rows    = A_num_columns;
  const auto AT_num_columns = A_num_rows;
  const auto BT_num_columns = B_num_rows;

  cblas_dgemm(layout, transA, transB, AT_num_rows, BT_num_columns, AT_num_columns, c, A_ptr, A_leading_dimension, B_ptr, B_leading_dimension, 0.0, result_ptr, BT_num_columns);
}

void Ax(double* result_ptr, const double* A_ptr, const double* x_ptr, const int num_rows, const int num_columns, const int leading_dimension, const int incr, const int incx)
{
  REQUIRE(0 < num_rows, "result_num_rows should be greater than 0");
  REQUIRE(0 < num_columns, "result_num_columns should be greater than 0");
  REQUIRE(num_columns <= leading_dimension, "leading_dimension error");

  const auto Layout = CBLAS_LAYOUT::CblasRowMajor;
  const auto trans  = CBLAS_TRANSPOSE::CblasNoTrans;

  cblas_dgemv(Layout, trans, num_rows, num_columns, 1.0, A_ptr, leading_dimension, x_ptr, incx, 0.0, result_ptr, incr);
}

void ATx(double* result_ptr, const double* A_ptr, const double* x_ptr, const int num_rows, const int num_columns, const int leading_dimension, const int incr, const int incx)
{
  REQUIRE(0 < num_rows, "result_num_rows should be greater than 0");
  REQUIRE(0 < num_columns, "result_num_columns should be greater than 0");
  REQUIRE(num_columns <= leading_dimension, "leading_dimension error");

  const auto Layout = CBLAS_LAYOUT::CblasRowMajor;
  const auto trans  = CBLAS_TRANSPOSE::CblasTrans;

  cblas_dgemv(Layout, trans, num_rows, num_columns, 1.0, A_ptr, leading_dimension, x_ptr, incx, 0.0, result_ptr, incr);
}

void cAx(double* result_ptr, const double c, const double* A_ptr, const double* x_ptr, const int num_rows, const int num_columns, const int leading_dimension, const int incr, const int incx)
{
  REQUIRE(0 < num_rows, "result_num_rows should be greater than 0");
  REQUIRE(0 < num_columns, "result_num_columns should be greater than 0");
  REQUIRE(num_columns <= leading_dimension, "leading_dimension error");

  const auto Layout = CBLAS_LAYOUT::CblasRowMajor;
  const auto trans  = CBLAS_TRANSPOSE::CblasNoTrans;

  cblas_dgemv(Layout, trans, num_rows, num_columns, c, A_ptr, leading_dimension, x_ptr, incx, 0.0, result_ptr, incr);
}

void cATx(double* result_ptr, const double c, const double* A_ptr, const double* x_ptr, const int num_rows, const int num_columns, const int leading_dimension, const int incr, const int incx)
{
  REQUIRE(0 < num_rows, "result_num_rows should be greater than 0");
  REQUIRE(0 < num_columns, "result_num_columns should be greater than 0");
  REQUIRE(num_columns <= leading_dimension, "leading_dimension error");

  const auto Layout = CBLAS_LAYOUT::CblasRowMajor;
  const auto trans  = CBLAS_TRANSPOSE::CblasTrans;

  cblas_dgemv(Layout, trans, num_rows, num_columns, c, A_ptr, leading_dimension, x_ptr, incx, 0.0, result_ptr, incr);
}

void invA(double* A_ptr, const int num_rows, const int num_columns, const int leading_dimension)
{
  const int        matrix_layout = LAPACK_ROW_MAJOR;
  const lapack_int m             = num_rows;
  const lapack_int n             = num_columns;
  const lapack_int lda           = leading_dimension;

  // PLU decomposition
  std::vector<int> ipiv(std::min(m, n));
  const lapack_int PLU_decomposition_routine_info = LAPACKE_dgetrf(matrix_layout, m, n, A_ptr, lda, ipiv.data());

  REQUIRE(0 <= PLU_decomposition_routine_info, "info should be greater than 0 when sucess PLU decomposition");

  // Inverse
  const lapack_int inverse_routine_info = LAPACKE_dgetri(matrix_layout, n, A_ptr, lda, ipiv.data());

  REQUIRE(inverse_routine_info == 0, "info should be 0 when success matrix get_inverse");
  // info > 0 "U is singular matrix in L-U decomposition"
  // info < 0 "fail to inverse the matrix"
}

} // namespace ms::math::blas::mkl

namespace ms::math::blas
{
void abs_x(double* result_ptr, const double* x_ptr, const int n, const int incr, const int incx)
{
  ms::math::blas::manual::abs_x(result_ptr, x_ptr, n, incr, incx);
}

void abs_x(double* x_ptr, const int n, const int incx)
{
  ms::math::blas::manual::abs_x(x_ptr, n, incx);
}

double abs_sum_x(const double* x_ptr, const int n, const int incx)
{
  if (n < ms::math::blas::asum_criteria)
  {
    return ms::math::blas::manual::abs_sum_x(x_ptr, n, incx);
  }
  else
  {
    return ms::math::blas::mkl::abs_sum_x(x_ptr, n, incx);
  }
}

void copy(double* result_ptr, const double* x_ptr, const int n, const int incr, const int incx)
{
  if (n < ms::math::blas::copy_criteria)
  {
    ms::math::blas::manual::copy(result_ptr, x_ptr, n, incr, incx);
  }
  else
  {
    ms::math::blas::mkl::copy(result_ptr, x_ptr, n, incr, incx);
  }
}

void cx(double* result_ptr, const double c, const double* x_ptr, const int n, const int incr, const int incx)
{
  if (n < ms::math::blas::scal_criteria)
  {
    ms::math::blas::manual::cx(result_ptr, c, x_ptr, n, incr, incx);
  }
  else
  {
    ms::math::blas::mkl::cx(result_ptr, c, x_ptr, n, incr, incx);
  }
}

void cx(const double c, double* x_ptr, const int n, const int incx)
{
  if (n < ms::math::blas::scal_criteria)
  {
    ms::math::blas::manual::cx(c, x_ptr, n, incx);
  }
  else
  {
    ms::math::blas::mkl::cx(c, x_ptr, n, incx);
  }
}

size_t find_maximum_element_pos(const double* x_ptr, const int n, const int incx)
{
  return cblas_idamax(n, x_ptr, incx);
}

double x_dot_y(const double* x_ptr, const double* y_ptr, const int n, const int incx, const int incy)
{
  if (n < ms::math::blas::dot_criteria)
  {
    return ms::math::blas::manual::x_dot_y(x_ptr, y_ptr, n, incx, incy);
  }
  else
  {
    return ms::math::blas::mkl::x_dot_y(x_ptr, y_ptr, n, incx, incy);
  }
}

void x_plus_y(double* result_ptr, const double* x_ptr, const double* y_ptr, const int n, const int incr, const int incx, const int incy)
{
  if (n < ms::math::blas::axpy_criteria)
  {
    ms::math::blas::manual::x_plus_y(result_ptr, x_ptr, y_ptr, n, incr, incx, incy);
  }
  else
  {
    ms::math::blas::mkl::x_plus_y(result_ptr, x_ptr, y_ptr, n, incr, incx, incy);
  }
}

void x_plus_assign_y(double* x_ptr, const double* y_ptr, const int n, const int incx, const int incy)
{
  if (n < ms::math::blas::axpy_criteria)
  {
    ms::math::blas::manual::x_plus_assign_y(x_ptr, y_ptr, n, incx, incy);
  }
  else
  {
    ms::math::blas::mkl::x_plus_assign_y(x_ptr, y_ptr, n, incx, incy);
  }
}

void x_plus_assign_cy(double* x_ptr, const double c, const double* y_ptr, const int n, const int incx, const int incy)
{
  if (n < ms::math::blas::axpy_criteria)
  {
    ms::math::blas::x_plus_assign_cy(x_ptr, c, y_ptr, n, incx, incy);
  }
  else
  {
    ms::math::blas::mkl::x_plus_assign_cy(x_ptr, c, y_ptr, n, incx, incy);
  }
}

void x_minus_y(double* result_ptr, const double* x_ptr, const double* y_ptr, const int n, const int incr, const int incx, const int incy)
{
  if (n < ms::math::blas::axpy_criteria)
  {
    ms::math::blas::manual::x_minus_y(result_ptr, x_ptr, y_ptr, n, incr, incx, incy);
  }
  else
  {
    ms::math::blas::mkl::x_minus_y(result_ptr, x_ptr, y_ptr, n, incr, incx, incy);
  }
}

void x_minus_assign_y(double* x_ptr, const double* y_ptr, const int n, const int incx, const int incy)
{
  if (n < ms::math::blas::axpy_criteria)
  {
    ms::math::blas::manual::x_minus_assign_y(x_ptr, y_ptr, n, incx, incy);
  }
  else
  {
    ms::math::blas::mkl::x_minus_assign_y(x_ptr, y_ptr, n, incx, incy);
  }
}

void cA(const double c, double* A_ptr, const int num_rows, const int num_columns, const int leading_dimension)
{
  if (num_rows * num_columns < ms::math::blas::cA_criteria)
  {
    ms::math::blas::manual::cA(c, A_ptr, num_rows, num_columns, leading_dimension);
  }
  else
  {
    ms::math::blas::mkl::cA(c, A_ptr, num_rows, num_columns, leading_dimension);
  }
}

void cA(double* result_ptr, const double c, const double* A_ptr, const int num_rows, const int num_columns, const int A_leading_dimension)
{
  if (num_rows * num_columns < ms::math::blas::cA_criteria)
  {
    ms::math::blas::manual::cA(result_ptr, c, A_ptr, num_rows, num_columns, A_leading_dimension);
  }
  else
  {
    ms::math::blas::mkl::cA(result_ptr, c, A_ptr, num_rows, num_columns, A_leading_dimension);
  }
}

void A_plus_cB(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int num_rows, const int num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension)
{
  if (num_columns < ms::math::blas::A_plus_cB_criteria)
  {
    ms::math::blas::manual::A_plus_cB(result_ptr, A_ptr, c, B_ptr, num_rows, num_columns, result_leading_dimension, A_leading_dimension, B_leading_dimension);
  }
  else
  {
    ms::math::blas::mkl::A_plus_cB(result_ptr, A_ptr, c, B_ptr, num_rows, num_columns, result_leading_dimension, A_leading_dimension, B_leading_dimension);
  }
}

void AT_plus_cB(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int num_rows, const int num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension)
{
  if (num_columns < ms::math::blas::AT_plus_cB_criteria)
  {
    ms::math::blas::manual::AT_plus_cB(result_ptr, A_ptr, c, B_ptr, num_rows, num_columns, result_leading_dimension, A_leading_dimension, B_leading_dimension);
  }
  else
  {
    ms::math::blas::mkl::AT_plus_cB(result_ptr, A_ptr, c, B_ptr, num_rows, num_columns, result_leading_dimension, A_leading_dimension, B_leading_dimension);
  }
}

void A_plus_cBT(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int num_rows, const int num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension)
{
  if (num_columns < ms::math::blas::A_plus_cBT_criteria)
  {
    ms::math::blas::manual::A_plus_cBT(result_ptr, A_ptr, c, B_ptr, num_rows, num_columns, result_leading_dimension, A_leading_dimension, B_leading_dimension);
  }
  else
  {
    ms::math::blas::mkl::A_plus_cBT(result_ptr, A_ptr, c, B_ptr, num_rows, num_columns, result_leading_dimension, A_leading_dimension, B_leading_dimension);
  }
}

void AT_plus_cBT(double* result_ptr, const double* A_ptr, const double c, const double* B_ptr, const int result_num_rows, const int result_num_columns, const int result_leading_dimension, const int A_leading_dimension, const int B_leading_dimension)
{
  if (result_num_columns < ms::math::blas::AT_plus_cBT_criteria)
  {
    ms::math::blas::manual::AT_plus_cBT(result_ptr, A_ptr, c, B_ptr, result_num_rows, result_num_columns, result_leading_dimension, A_leading_dimension, B_leading_dimension);
  }
  else
  {
    ms::math::blas::mkl::AT_plus_cBT(result_ptr, A_ptr, c, B_ptr, result_num_rows, result_num_columns, result_leading_dimension, A_leading_dimension, B_leading_dimension);
  }
}

void cAB(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_columns, const int A_leading_dimension, const int B_leading_dimension)
{
  if (A_num_rows < ms::math::blas::cAB_A_rows_criteria && A_num_rows * A_num_columns < ms::math::blas::cAB_size_criteria)
  {
    ms::math::blas::manual::cAB(result_ptr, c, A_ptr, B_ptr, A_num_rows, A_num_columns, B_num_columns, A_leading_dimension, B_leading_dimension);
  }
  else
  {
    ms::math::blas::mkl::cAB(result_ptr, c, A_ptr, B_ptr, A_num_rows, A_num_columns, B_num_columns, A_leading_dimension, B_leading_dimension);
  }
}

void cATB(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_columns, const int A_leading_dimension, const int B_leading_dimension)
{
  if (A_num_rows < ms::math::blas::cATB_A_columns_criteria && A_num_rows * A_num_columns < ms::math::blas::cATB_size_criteria)
  {
    ms::math::blas::manual::cATB(result_ptr, c, A_ptr, B_ptr, A_num_rows, A_num_columns, B_num_columns, A_leading_dimension, B_leading_dimension);
  }
  else
  {
    ms::math::blas::mkl::cATB(result_ptr, c, A_ptr, B_ptr, A_num_rows, A_num_columns, B_num_columns, A_leading_dimension, B_leading_dimension);
  }
}

void cABT(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_columns, const int A_leading_dimension, const int B_leading_dimension)
{
  if (A_num_rows < ms::math::blas::cABT_A_rows_criteria && A_num_rows * A_num_columns < ms::math::blas::cABT_size_criteria)
  {
    ms::math::blas::manual::cABT(result_ptr, c, A_ptr, B_ptr, A_num_rows, A_num_columns, B_num_columns, A_leading_dimension, B_leading_dimension);
  }
  else
  {
    ms::math::blas::mkl::cABT(result_ptr, c, A_ptr, B_ptr, A_num_rows, A_num_columns, B_num_columns, A_leading_dimension, B_leading_dimension);
  }
}

void cATBT(double* result_ptr, const double c, const double* A_ptr, const double* B_ptr, const int A_num_rows, const int A_num_columns, const int B_num_columns, const int A_leading_dimension, const int B_leading_dimension)
{
  if (A_num_rows < ms::math::blas::cATBT_A_columns_criteria && A_num_rows * A_num_columns < ms::math::blas::cATBT_size_criteria)
  {
    ms::math::blas::manual::cATBT(result_ptr, c, A_ptr, B_ptr, A_num_rows, A_num_columns, B_num_columns, A_leading_dimension, B_leading_dimension);
  }
  else
  {
    ms::math::blas::mkl::cATBT(result_ptr, c, A_ptr, B_ptr, A_num_rows, A_num_columns, B_num_columns, A_leading_dimension, B_leading_dimension);
  }
}

void Ax(double* result_ptr, const double* A_ptr, const double* x_ptr, const int num_rows, const int num_columns, const int leading_dimension, const int incr, const int incx)
{
  if (num_rows * num_columns < ms::math::blas::gemv_criteria)
  {
    ms::math::blas::manual::Ax(result_ptr, A_ptr, x_ptr, num_rows, num_columns, leading_dimension, incr, incx);
  }
  else
  {
    ms::math::blas::mkl::Ax(result_ptr, A_ptr, x_ptr, num_rows, num_columns, leading_dimension, incr, incx);
  }
}

void ATx(double* result_ptr, const double* A_ptr, const double* x_ptr, const int num_rows, const int num_columns, const int leading_dimension, const int incr, const int incx)
{
  if (num_rows * num_columns < ms::math::blas::gemv_criteria)
  {
    ms::math::blas::manual::ATx(result_ptr, A_ptr, x_ptr, num_rows, num_columns, leading_dimension, incr, incx);
  }
  else
  {
    ms::math::blas::mkl::ATx(result_ptr, A_ptr, x_ptr, num_rows, num_columns, leading_dimension, incr, incx);
  }
}

void cAx(double* result_ptr, const double c, const double* A_ptr, const double* x_ptr, const int num_rows, const int num_columns, const int leading_dimension, const int incr, const int incx)
{
  if (num_rows * num_columns < ms::math::blas::gemv_criteria)
  {
    ms::math::blas::manual::cAx(result_ptr, c, A_ptr, x_ptr, num_rows, num_columns, leading_dimension, incr, incx);
  }
  else
  {
    ms::math::blas::mkl::cAx(result_ptr, c, A_ptr, x_ptr, num_rows, num_columns, leading_dimension, incr, incx);
  }
}

void cATx(double* result_ptr, const double c, const double* A_ptr, const double* x_ptr, const int num_rows, const int num_columns, const int leading_dimension, const int incr, const int incx)
{
  if (num_rows * num_columns < ms::math::blas::gemv_criteria)
  {
    ms::math::blas::manual::cATx(result_ptr, c, A_ptr, x_ptr, num_rows, num_columns, leading_dimension, incr, incx);
  }
  else
  {
    ms::math::blas::mkl::cATx(result_ptr, c, A_ptr, x_ptr, num_rows, num_columns, leading_dimension, incr, incx);
  }
}

void invA(double* A_ptr, const int num_rows, const int num_columns, const int leading_dimension)
{
  ms::math::blas::mkl::invA(A_ptr, num_rows, num_columns, leading_dimension);
}

} // namespace ms::math::blas


