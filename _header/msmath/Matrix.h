#pragma once
#include <concepts>
#include <string>
#include <vector>

#include "Vector.h"

// forward declaration
namespace ms::math
{

class Matrix;

} // namespace ms::math

/*










*/

// Matrix_View class declaration
namespace ms::math
{
class Matrix_View
{
public:
  Matrix_View(void) = default;
  template <const_span T>
  Matrix_View(const int num_rows, const int num_columns, const T& values, const int leading_dimension = -1);

public:
  Matrix    operator+(const Matrix_View other) const;
  Matrix    operator-(const Matrix_View other) const;
  Matrix    operator*(const double constant) const;
  Vector<0> operator*(const Vector_View vec) const;
  Matrix    operator*(const Matrix_View other) const;
  bool      operator==(const Matrix_View other) const;

public:
  double              at(const int row_index, const int column_index) const;
  Vector_View         column_vector_view(const int index) const;
  const double*       data(void) const;
  bool                has_compact_data(void) const;
  bool                is_transposed(void) const;
  Matrix              inverse(void) const;
  int                 leading_dimension(void) const;
  int                 num_columns(void) const;
  int                 num_rows(void) const;
  int                 num_values(void) const;
  Matrix_View         sub_view(const int start_row_index, const int end_row_index, const int start_column_index, const int end_column_index) const;
  Vector_View         row_vector_view(const int index) const;
  std::pair<int, int> size(void) const;
  Matrix_View         transpose(void) const;
  std::string         to_string(void) const;
  Matrix_View         view(void) const;

protected:
  bool                    is_square_matrix(void) const;
  bool                    is_valid_column_index(const int index) const;
  bool                    is_valid_column_range(const int start_index, const int end_index) const;
  bool                    is_valid_row_index(const int index) const;
  bool                    is_valid_row_range(const int start_index, const int end_index) const;
  int                     index_to_position(const int row_index, const int column_index) const;
  std::span<const double> values_part_view(const int start_row_index, const int end_row_index, const int start_column_index, const int end_column_index) const;
  std::vector<double>     values_vector(void) const;

protected:
  int                     _num_rows          = 0;
  int                     _num_columns       = 0;
  int                     _leading_dimension = 0;
  bool                    _is_transposed     = false;
  std::span<const double> _values_view;
};

} // namespace ms::math

/*










*/

// Matrix_Wrap class declaration
namespace ms::math
{

class Matrix_Wrap : public Matrix_View
{
public:
  Matrix_Wrap(void) = default;

  template <const_span T>
  Matrix_Wrap(const int num_rows, const int num_columns, T& values, const int leading_dimension = -1)
      : Matrix_View(num_rows, num_columns, values, leading_dimension),
        _values_wrap(values){};

public:
  void operator*=(const double constant);
  void operator+=(const Matrix_View other);
  void operator-=(const Matrix_View other);

public:
  double&     at(const int row, const int column);
  void        be_inverse(void);
  void        be_transpose(void);
  void        change_column(const int index, const double* values);
  void        change_row(const int index, const double* values);
  void        change_value(const Matrix_View other);
  Vector_Wrap column_vector_wrap(const int index);
  double*     data(void);
  Matrix_Wrap sub_wrap(const int start_row_index, const int end_row_index, const int start_col_index, const int end_col_index);
  Vector_Wrap row_vector_wrap(const int index);
  Matrix_Wrap wrap(void);

  // resolve parent method hiding problem (overloading across scope problem)
public:
  using Matrix_View::at;
  using Matrix_View::data;

private:
  std::span<double> values_part_wrap(const int start_row_index, const int end_row_index, const int start_column_index, const int end_column_index);

protected:
  std::span<double> _values_wrap;
};

} // namespace ms::math

/*










*/

// Matrix class declaration
namespace ms::math
{
class Matrix : public Matrix_Wrap
{
public:
  Matrix(void) = default;
  Matrix(const int matrix_order);
  Matrix(const int num_row, const int num_column);
  Matrix(const int matrix_order, const std::vector<double>& value);
  Matrix(const int num_row, const int num_column, std::vector<double>&& values);
  Matrix(const Matrix& other);
  Matrix(Matrix&& other) noexcept;

public:
  void operator=(const Matrix& other);
  void operator=(Matrix&& other) noexcept;

public:
  void resize(const int row, const int column);

private:
  void reallocate_values(void);

private:
  std::vector<double> _values;
};

} // namespace ms::math

/*










*/

// free function declaration
namespace ms::math
{

Matrix        operator*(const double constant, const Matrix_View M);
std::ostream& operator<<(std::ostream& os, const Matrix_View m);

} // namespace ms::math

namespace ms::math::blas
{
void mpm(Matrix_Wrap result, const Matrix_View m1, const Matrix_View m2);
void cmv(Vector_Wrap result, const double c, const Matrix_View m, const Vector_View v);
} // namespace ms::math::blas

/*










*/

// template definition for class Matrix_Values_Const_Wrapper
namespace ms::math
{

template <const_span T>
Matrix_View::Matrix_View(const int num_rows, const int num_columns, const T& values, const int leading_dimension)
    : _num_rows(num_rows),
      _num_columns(num_columns),
      _leading_dimension(leading_dimension),
      _values_view(values)
{
  REQUIRE(num_rows > 0 && num_columns > 0, "number of rows or number of columns should be positive");
  REQUIRE(this->num_values() <= this->_values_view.size(), "given values can't be smaller then given size");

  if (this->_leading_dimension == -1)
  {
    this->_leading_dimension = this->_num_columns;
  }

  REQUIRE(this->_num_columns <= this->_leading_dimension, "leading dimension shold be larger than num_cols");
};

} // namespace ms::math