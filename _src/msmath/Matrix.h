#pragma once
#include <concepts>
#include <string>
#include <vector>

// forward declaration
namespace ms::math
{
class Vector_Const_Wrapper;
class Vector_Wrapper;

template <int Dim>
class Vector;

class Matrix;
} // namespace ms::math

// class declaration
namespace ms::math
{

enum class Transpose_Type
{
  not_transposed,
  transposed
};

class Matrix_Const_Wrapper
{
public:
  Matrix_Const_Wrapper(void) = default;
  Matrix_Const_Wrapper(const int num_rows, const int num_columns, const double* ptr, const Transpose_Type transpose_type);
  Matrix_Const_Wrapper(const int num_rows, const int num_columns, const double* ptr, const int leading_dimension = -1, const Transpose_Type transpose_type = Transpose_Type::not_transposed);

public:
  Matrix    operator+(const Matrix_Const_Wrapper& other) const;
  Matrix    operator-(const Matrix_Const_Wrapper& other) const;
  Matrix    operator*(const double constant) const;
  Vector<0> operator*(const Vector_Const_Wrapper& vec) const;
  Matrix    operator*(const Matrix_Const_Wrapper& other) const;
  bool      operator==(const Matrix_Const_Wrapper& other) const;

public:
  double                    at(const int row_index, const int column_index) const;
  Vector_Const_Wrapper      const_column_vector_wrapper(const int index) const;
  Vector_Const_Wrapper      const_row_vector_wrapper(const int index) const;
  Matrix_Const_Wrapper      const_part(const int start_row_index, const int end_row_index, const int start_column_index, const int end_column_index);
  const double*             data(void) const;
  bool                      has_compact_data(void) const;
  bool                      is_transposed(void) const;
  Matrix                    inverse_matrix(void) const;
  size_t                    num_cols(void) const;
  size_t                    num_rows(void) const;
  size_t                    num_values(void) const;
  int                       leading_dimension(void) const;
  std::pair<size_t, size_t> size(void) const;
  std::string               to_string(void) const;
  Matrix_Const_Wrapper      transposed_const_matrix_wrapper(void) const;
  Matrix                    transposed_matrix(void) const;

protected:
  bool                is_square_matrix(void) const;
  bool                is_valid_column_index(const int index) const;
  bool                is_valid_column_range(const int start_index, const int end_index) const;
  bool                is_valid_row_index(const int index) const;
  bool                is_valid_row_range(const int start_index, const int end_index) const;
  Transpose_Type      transpose_type(void) const;
  std::vector<double> values_vector(void) const;

protected:
  int           num_rows_          = 0;
  int           num_columns_       = 0;
  int           leading_dimension_ = 0;
  const double* const_data_ptr_    = nullptr;
  bool          is_transposed_     = false;
};

class Matrix_Wrapper : public Matrix_Const_Wrapper
{
public:
  Matrix_Wrapper(void) = default;
  Matrix_Wrapper(const int num_row, const int num_column, double* ptr, const Transpose_Type transpose_type)
      : Matrix_Const_Wrapper(num_row, num_column, ptr, transpose_type),
        data_ptr_(ptr){};
  Matrix_Wrapper(const int num_row, const int num_column, double* ptr, const int leading_dimension = -1, const Transpose_Type transpose_type = Transpose_Type::not_transposed)
      : Matrix_Const_Wrapper(num_row, num_column, ptr, leading_dimension, transpose_type),
        data_ptr_(ptr){};

public:
  void operator*=(const double constant);
  void operator+=(const Matrix_Const_Wrapper& other);
  void operator-=(const Matrix_Const_Wrapper& other);

public:
  double&        at(const int row, const int column);
  void           change_column(const int index, const double* values);
  void           change_row(const int index, const double* values);
  void           change_value(const Matrix_Const_Wrapper& cmw);
  Vector_Wrapper column_wrapper(const int index);
  double*        data(void);
  void           inverse(void);
  Matrix_Wrapper part(const int start_row_index, const int end_row_index, const int start_col_index, const int end_col_index);
  Vector_Wrapper row_wrapper(const int index);
  void           transpose(void);

  // resolve parent method hiding problem (overloading across scope problem)
public:
  using Matrix_Const_Wrapper::at;
  using Matrix_Const_Wrapper::data;

protected:
  double* data_ptr_ = nullptr;
};

class Matrix : public Matrix_Wrapper
{
public:
  Matrix(void) = default;
  Matrix(const int matrix_order);
  Matrix(const int matrix_order, const std::vector<double>& value);
  Matrix(const int num_row, const int num_column, const Transpose_Type transpose_type = Transpose_Type::not_transposed);
  Matrix(const int num_row, const int num_column, const double* value_ptr);
  Matrix(const int num_row, const int num_column, std::vector<double>&& values, const Transpose_Type transpose_type = Transpose_Type::not_transposed);
  Matrix(const Matrix& other);
  Matrix(Matrix&& other) noexcept;

public:
  void operator=(const Matrix& other);
  void operator=(Matrix&& other) noexcept;

private:
  void essential_code_for_constructor(void);

private:
  std::vector<double> values_;
};

} // namespace ms::math

// free function declaration
namespace ms::math
{
Matrix        operator*(const double constant, const Matrix_Const_Wrapper& M);
std::ostream& operator<<(std::ostream& os, const Matrix_Const_Wrapper& m);
std::ostream& operator<<(std::ostream& os, const Matrix& m);
} // namespace ms::math