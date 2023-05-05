#include "Matrix.h"

#include "BLAS.h"
#include "Vector.h"

#include "msexception/Exception.h"
#include <iomanip>

namespace ms::math
{

Matrix_Const_Wrapper::Matrix_Const_Wrapper(const int num_rows, const int num_columns, const double* ptr, const Transpose_Type transpose_type)
    : num_rows_(num_rows), num_columns_(num_columns), const_data_ptr_(ptr), leading_dimension_(num_columns)
{
  REQUIRE(num_rows > 0 && num_columns > 0, "number of rows or number of columns should be positive");

  if (transpose_type == Transpose_Type::transposed)
  {
    this->is_transposed_ = true;
  }
};

Matrix_Const_Wrapper::Matrix_Const_Wrapper(const int num_rows, const int num_columns, const double* ptr, const int leading_dimension, const Transpose_Type transpose_type)
    : num_rows_(num_rows), num_columns_(num_columns), const_data_ptr_(ptr), leading_dimension_(leading_dimension)
{
  REQUIRE(num_rows > 0 && num_columns > 0, "number of rows or number of columns should be positive");

  if (transpose_type == Transpose_Type::transposed)
  {
    this->is_transposed_ = true;
  }

  if (this->leading_dimension_ == -1)
  {
    this->leading_dimension_ = this->num_columns_;
  }

  REQUIRE(this->num_columns_ <= this->leading_dimension_, "leading dimension shold be larger than num_cols");
};

Matrix Matrix_Const_Wrapper::operator+(const Matrix_Const_Wrapper& other) const
{
  REQUIRE(this->size() == other.size(), "two matrix should be same size");

  const auto result_num_rows    = static_cast<int>(this->num_rows());
  const auto result_num_columns = static_cast<int>(this->num_cols());
  Matrix     result(result_num_rows, result_num_columns);

  constexpr double c = 1.0;

  if (this->is_transposed_)
  {
    if (other.is_transposed_)
    {
      ms::math::blas::AT_plus_cBT(result.data(), this->const_data_ptr_, c, other.const_data_ptr_, result_num_rows, result_num_columns, result_num_columns, this->leading_dimension_, other.leading_dimension_);
    }
    else
    {
      ms::math::blas::AT_plus_cB(result.data(), this->const_data_ptr_, c, other.const_data_ptr_, result_num_rows, result_num_columns, result_num_columns, this->leading_dimension_, other.leading_dimension_);
    }
  }
  else
  {
    if (other.is_transposed_)
    {
      ms::math::blas::A_plus_cBT(result.data(), this->const_data_ptr_, c, other.const_data_ptr_, result_num_rows, result_num_columns, result_num_columns, this->leading_dimension_, other.leading_dimension_);
    }
    else
    {
      ms::math::blas::A_plus_cB(result.data(), this->const_data_ptr_, c, other.const_data_ptr_, result_num_rows, result_num_columns, result_num_columns, this->leading_dimension_, other.leading_dimension_);
    }
  }

  return result;
}

Matrix Matrix_Const_Wrapper::operator-(const Matrix_Const_Wrapper& other) const
{
  REQUIRE(this->size() == other.size(), "two matrix should be same size");

  const auto result_num_rows    = static_cast<int>(this->num_rows());
  const auto result_num_columns = static_cast<int>(this->num_cols());
  Matrix     result(result_num_rows, result_num_columns);

  constexpr double c = -1.0;

  if (this->is_transposed_)
  {
    if (other.is_transposed_)
    {
      ms::math::blas::AT_plus_cBT(result.data(), this->const_data_ptr_, c, other.const_data_ptr_, result_num_rows, result_num_columns, result_num_columns, this->leading_dimension_, other.leading_dimension_);
    }
    else
    {
      ms::math::blas::AT_plus_cB(result.data(), this->const_data_ptr_, c, other.const_data_ptr_, result_num_rows, result_num_columns, result_num_columns, this->leading_dimension_, other.leading_dimension_);
    }
  }
  else
  {
    if (other.is_transposed_)
    {
      ms::math::blas::A_plus_cBT(result.data(), this->const_data_ptr_, c, other.const_data_ptr_, result_num_rows, result_num_columns, result_num_columns, this->leading_dimension_, other.leading_dimension_);
    }
    else
    {
      ms::math::blas::A_plus_cB(result.data(), this->const_data_ptr_, c, other.const_data_ptr_, result_num_rows, result_num_columns, result_num_columns, this->leading_dimension_, other.leading_dimension_);
    }
  }

  return result;
}

Matrix Matrix_Const_Wrapper::operator*(const double constant) const
{
  Matrix result(this->num_rows_, this->num_columns_, this->transpose_type());

  ms::math::blas::cA(result.data(), constant, this->const_data_ptr_, this->num_rows_, this->num_columns_, this->leading_dimension_);

  return result;
}

Vector<0> Matrix_Const_Wrapper::operator*(const Vector_Const_Wrapper& vec) const
{
  if (this->is_transposed_)
  {
    REQUIRE(this->num_rows_ == vec.dimension(), "size should be mathced");

    Vector result(this->num_columns_);
    ms::math::blas::ATx(result.data(), this->const_data_ptr_, vec.data(), this->num_rows_, this->num_columns_, this->leading_dimension_, result.inc(), vec.inc());

    return result;
  }
  else
  {
    REQUIRE(this->num_columns_ == vec.dimension(), "size should be mathced");

    Vector result(this->num_rows_);
    ms::math::blas::Ax(result.data(), this->const_data_ptr_, vec.data(), this->num_rows_, this->num_columns_, this->leading_dimension_, result.inc(), vec.inc());

    return result;
  }
}

Matrix Matrix_Const_Wrapper::operator*(const Matrix_Const_Wrapper& other) const
{
  REQUIRE(this->num_cols() == other.num_rows(), "two matrix size should be mathched");

  const auto result_num_rows    = static_cast<int>(this->num_rows());
  const auto result_num_columns = static_cast<int>(other.num_cols());
  Matrix     result(result_num_rows, result_num_columns);

  constexpr double c = 1.0;

  if (this->is_transposed_)
  {
    if (other.is_transposed_)
    {
      ms::math::blas::cATBT(result.data(), c, this->const_data_ptr_, other.const_data_ptr_, this->num_rows_, this->num_columns_, other.num_rows_, this->leading_dimension_, other.leading_dimension_);
    }
    else
    {
      ms::math::blas::cATB(result.data(), c, this->const_data_ptr_, other.const_data_ptr_, this->num_rows_, this->num_columns_, other.num_columns_, this->leading_dimension_, other.leading_dimension_);
    }
  }
  else
  {
    if (other.is_transposed_)
    {
      ms::math::blas::cABT(result.data(), c, this->const_data_ptr_, other.const_data_ptr_, this->num_rows_, this->num_columns_, other.num_rows_, this->leading_dimension_, other.leading_dimension_);
    }
    else
    {
      ms::math::blas::cAB(result.data(), c, this->const_data_ptr_, other.const_data_ptr_, this->num_rows_, this->num_columns_, other.num_columns_, this->leading_dimension_, other.leading_dimension_);
    }
  }

  return result;
}

bool Matrix_Const_Wrapper::operator==(const Matrix_Const_Wrapper& other) const
{
  if (this->size() != other.size())
  {
    return false;
  }

  const auto [num_rows, num_columns] = this->size();

  for (int i = 0; i < num_rows; i++)
  {
    for (int j = 0; j < num_columns; j++)
    {
      if (this->at(i, j) != other.at(i, j))
        return false;
    }
  }

  return true;
}

double Matrix_Const_Wrapper::at(const int row_index, const int column_index) const
{
  REQUIRE(this->is_valid_row_index(row_index) && this->is_valid_column_index(column_index), "Given indexes should be valid");

  if (this->is_transposed_)
  {
    return this->const_data_ptr_[column_index * this->leading_dimension_ + row_index];
  }
  else
  {
    return this->const_data_ptr_[row_index * this->leading_dimension_ + column_index];
  }
}

Vector_Const_Wrapper Matrix_Const_Wrapper::const_column_vector_wrapper(const int index) const
{
  REQUIRE(this->is_valid_column_index(index), "Given index should be valid");

  if (this->is_transposed_)
  {
    const auto start_data_ptr = this->const_data_ptr_ + index * this->leading_dimension_;

    Vector_Const_Wrapper column_vector(start_data_ptr, this->num_columns_);
    return column_vector;
  }
  else
  {
    const auto start_data_ptr = this->const_data_ptr_ + index;

    Vector_Const_Wrapper column_vector(start_data_ptr, this->num_rows_, this->leading_dimension_);
    return column_vector;
  }
}

Vector_Const_Wrapper Matrix_Const_Wrapper::const_row_vector_wrapper(const int index) const
{
  REQUIRE(this->is_valid_row_index(index), "Given index should be valid");

  if (this->is_transposed_)
  {
    const auto start_data_ptr = this->const_data_ptr_ + index;

    Vector_Const_Wrapper row_vector(start_data_ptr, this->num_rows_, this->leading_dimension_);
    return row_vector;
  }
  else
  {
    const auto start_data_ptr = this->const_data_ptr_ + index * this->leading_dimension_;

    Vector_Const_Wrapper row_vector(start_data_ptr, this->num_columns_);
    return row_vector;
  }
}

Matrix_Const_Wrapper Matrix_Const_Wrapper::const_part(const int start_row_index, const int end_row_index, const int start_column_index, const int end_column_index)
{
  REQUIRE(this->is_valid_row_range(start_row_index, end_row_index), "Given row range should be valid");
  REQUIRE(this->is_valid_column_range(start_column_index, end_column_index), "Given column range should be valid");

  const auto num_row    = end_row_index - start_row_index;
  const auto num_column = end_column_index - start_column_index;

  if (this->is_transposed_)
  {
    const auto start_data_ptr = this->const_data_ptr_ + start_column_index * this->leading_dimension_ + start_row_index;

    Matrix_Const_Wrapper part(num_column, num_row, start_data_ptr, this->leading_dimension_, Transpose_Type::transposed);
    return part;
  }
  else
  {
    const auto start_data_ptr = this->const_data_ptr_ + start_row_index * this->leading_dimension_ + start_column_index;

    Matrix_Const_Wrapper part(num_row, num_column, start_data_ptr, this->leading_dimension_, Transpose_Type::not_transposed);
    return part;
  }
}

const double* Matrix_Const_Wrapper::data(void) const
{
  return this->const_data_ptr_;
}

bool Matrix_Const_Wrapper::has_compact_data(void) const
{
  return this->leading_dimension_ == this->num_columns_;
}

bool Matrix_Const_Wrapper::is_transposed(void) const
{
  return this->is_transposed_;
}

Matrix Matrix_Const_Wrapper::inverse_matrix(void) const
{
  REQUIRE(this->is_square_matrix(), "matrix should be square matrix");
  REQUIRE(!this->is_transposed_, "matrix should not be transposed");

  Matrix result(this->num_rows_, this->num_columns_, this->values_vector());
  result.inverse();

  return result;
}

size_t Matrix_Const_Wrapper::num_cols(void) const
{
  if (this->is_transposed_)
  {
    return this->num_rows_;
  }
  else
  {
    return this->num_columns_;
  }
}

size_t Matrix_Const_Wrapper::num_rows(void) const
{
  if (this->is_transposed_)
  {
    return this->num_columns_;
  }
  else
  {
    return this->num_rows_;
  }
}

size_t Matrix_Const_Wrapper::num_values(void) const
{
  return this->num_rows_ * this->num_columns_;
}

int Matrix_Const_Wrapper::leading_dimension(void) const
{
  return this->leading_dimension_;
}

std::pair<size_t, size_t> Matrix_Const_Wrapper::size(void) const
{
  if (this->is_transposed_)
  {
    return {this->num_columns_, this->num_rows_};
  }
  else
  {
    return {this->num_rows_, this->num_columns_};
  }
}

std::string Matrix_Const_Wrapper::to_string(void) const
{
  std::ostringstream oss;
  oss << std::setprecision(16) << std::showpoint << std::left;
  for (int i = 0; i < this->num_rows(); ++i)
  {
    for (int j = 0; j < this->num_cols(); ++j)
      oss << std::setw(25) << this->at(i, j);
    oss << "\n";
  }
  return oss.str();
}

Matrix_Const_Wrapper Matrix_Const_Wrapper::transposed_const_matrix_wrapper(void) const
{
  auto transpose_type = Transpose_Type::transposed;

  if (this->is_transposed_)
  {
    transpose_type = Transpose_Type::not_transposed;
  }

  Matrix_Const_Wrapper result(this->num_rows_, this->num_columns_, this->const_data_ptr_, this->leading_dimension_, transpose_type);

  return result;
}

Matrix Matrix_Const_Wrapper::transposed_matrix(void) const
{
  auto transpose_type = Transpose_Type::transposed;

  if (this->is_transposed_)
  {
    transpose_type = Transpose_Type::not_transposed;
  }

  Matrix result(this->num_rows_, this->num_columns_, this->values_vector(), transpose_type);
  return result;
}

bool Matrix_Const_Wrapper::is_square_matrix(void) const
{
  return this->num_rows_ == this->num_columns_;
}

bool Matrix_Const_Wrapper::is_valid_column_index(const int index) const
{
  if (index < 0)
  {
    return false;
  }

  if (this->is_transposed_)
  {
    return index < this->num_rows_;
  }
  else
  {
    return index < this->num_columns_;
  }
}

bool Matrix_Const_Wrapper::is_valid_column_range(const int start_index, const int end_index) const
{
  if (end_index <= start_index)
  {
    return false;
  }

  if (start_index < 0)
  {
    return false;
  }

  if (this->is_transposed_)
  {
    return end_index <= this->num_rows_;
  }
  else
  {
    return end_index <= this->num_columns_;
  }
}

bool Matrix_Const_Wrapper::is_valid_row_index(const int index) const
{
  if (index < 0)
  {
    return false;
  }

  if (this->is_transposed_)
  {
    return index < this->num_columns_;
  }
  else
  {
    return index < this->num_rows_;
  }
}

bool Matrix_Const_Wrapper::is_valid_row_range(const int start_index, const int end_index) const
{
  if (end_index <= start_index)
  {
    return false;
  }

  if (start_index < 0)
  {
    return false;
  }

  if (this->is_transposed_)
  {
    return end_index <= this->num_columns_;
  }
  else
  {
    return end_index <= this->num_rows_;
  }
}

Transpose_Type Matrix_Const_Wrapper::transpose_type(void) const
{
  if (this->is_transposed_)
  {
    return Transpose_Type::transposed;
  }
  else
  {
    return Transpose_Type::not_transposed;
  }
}

std::vector<double> Matrix_Const_Wrapper::values_vector(void) const
{
  const auto          n = static_cast<int>(this->num_values());
  std::vector<double> result(n);

  if (this->has_compact_data())
  {
    ms::math::blas::copy(result.data(), this->const_data_ptr_, n);
  }
  else
  {
    auto result_ptr = result.data();
    auto value_ptr  = this->const_data_ptr_;
    for (int i = 0; i < this->num_rows_; i++)
    {
      ms::math::blas::copy(result_ptr, value_ptr, this->num_columns_);
      result_ptr += this->num_columns_;
      value_ptr += this->leading_dimension_;
    }
  }

  return result;
}

void Matrix_Wrapper::operator*=(const double constant)
{
  ms::math::blas::cA(constant, this->data_ptr_, this->num_rows_, this->num_columns_, this->leading_dimension_);
}

void Matrix_Wrapper::operator+=(const Matrix_Const_Wrapper& other)
{
  REQUIRE(this->size() == other.size(), "two matrix should be same size");
  REQUIRE(!this->is_transposed_, "matrix should not be transposed");

  constexpr double c = 1.0;

  if (other.is_transposed())
  {
    ms::math::blas::A_plus_cBT(this->data_ptr_, this->const_data_ptr_, c, other.data(), this->num_rows_, this->num_columns_, this->leading_dimension_, this->leading_dimension_, other.leading_dimension());
  }
  else
  {
    ms::math::blas::A_plus_cB(this->data_ptr_, this->const_data_ptr_, c, other.data(), this->num_rows_, this->num_columns_, this->leading_dimension_, this->leading_dimension_, other.leading_dimension());
  }
}

void Matrix_Wrapper::operator-=(const Matrix_Const_Wrapper& other)
{
  REQUIRE(this->size() == other.size(), "two matrix should be same size");
  REQUIRE(!this->is_transposed_, "matrix should not be transposed");

  constexpr double c = -1.0;

  if (other.is_transposed())
  {
    ms::math::blas::A_plus_cBT(this->data_ptr_, this->const_data_ptr_, c, other.data(), this->num_rows_, this->num_columns_, this->leading_dimension_, this->leading_dimension_, other.leading_dimension());
  }
  else
  {
    ms::math::blas::A_plus_cB(this->data_ptr_, this->const_data_ptr_, c, other.data(), this->num_rows_, this->num_columns_, this->leading_dimension_, this->leading_dimension_, other.leading_dimension());
  }
}

double& Matrix_Wrapper::at(const int row_index, const int column_index)
{
  REQUIRE(this->is_valid_row_index(row_index) && this->is_valid_column_index(column_index), "Given indexes should be valid");

  if (this->is_transposed_)
  {
    return this->data_ptr_[column_index * this->leading_dimension_ + row_index];
  }
  else
  {
    return this->data_ptr_[row_index * this->leading_dimension_ + column_index];
  }
}

void Matrix_Wrapper::change_column(const int index, const double* values)
{
  const auto num_cols = this->num_cols();
  auto       column   = this->column_wrapper(index);

  for (int i = 0; i < num_cols; ++i)
  {
    column[i] = values[i];
  }
}

void Matrix_Wrapper::change_row(const int index, const double* values)
{
  const auto num_rows = this->num_rows();
  auto       row      = this->row_wrapper(index);

  for (int i = 0; i < num_rows; ++i)
  {
    row[i] = values[i];
  }
}

void Matrix_Wrapper::change_value(const Matrix_Const_Wrapper& cmw)
{
  REQUIRE(this->size() == cmw.size(), "size should be matched");

  for (int i = 0; i < this->num_rows(); i++)
  {
    for (int j = 0; j < this->num_cols(); j++)
    {
      this->at(i, j) = cmw.at(i, j);
    }
  }
}

Vector_Wrapper Matrix_Wrapper::column_wrapper(const int index)
{
  REQUIRE(this->is_valid_column_index(index), "Given index should be valid");

  if (this->is_transposed_)
  {
    const auto     start_data_ptr = this->data_ptr_ + index * this->leading_dimension_;
    Vector_Wrapper column_vector(start_data_ptr, this->num_columns_);

    return column_vector;
  }
  else
  {
    const auto     start_data_ptr = this->data_ptr_ + index;
    Vector_Wrapper column_vector(start_data_ptr, this->num_rows_, this->leading_dimension_);

    return column_vector;
  }
}

double* Matrix_Wrapper::data(void)
{
  return this->data_ptr_;
}

Vector_Wrapper Matrix_Wrapper::row_wrapper(const int index)
{
  REQUIRE(this->is_valid_row_index(index), "Given index should be valid");

  if (this->is_transposed_)
  {
    const auto     start_data_ptr = this->data_ptr_ + index;
    Vector_Wrapper row_vector(start_data_ptr, this->num_rows_, this->leading_dimension_);

    return row_vector;
  }
  else
  {
    const auto     start_data_ptr = this->data_ptr_ + index * this->leading_dimension_;
    Vector_Wrapper row_vector(start_data_ptr, this->num_columns_);

    return row_vector;
  }
}

Matrix_Wrapper Matrix_Wrapper::part(const int start_row_index, const int end_row_index, const int start_column_index, const int end_column_index)
{
  REQUIRE(this->is_valid_row_range(start_row_index, end_row_index), "Given row range should be valid");
  REQUIRE(this->is_valid_column_range(start_column_index, end_column_index), "Given column range should be valid");

  const auto num_row    = end_row_index - start_row_index;
  const auto num_column = end_column_index - start_column_index;

  if (this->is_transposed_)
  {
    const auto start_data_ptr = this->data_ptr_ + start_column_index * this->leading_dimension_ + start_row_index;

    Matrix_Wrapper part(num_column, num_row, start_data_ptr, this->leading_dimension_, Transpose_Type::transposed);
    return part;
  }
  else
  {
    const auto start_data_ptr = this->data_ptr_ + start_row_index * this->leading_dimension_ + start_column_index;

    Matrix_Wrapper part(num_row, num_column, start_data_ptr, this->leading_dimension_, Transpose_Type::not_transposed);
    return part;
  }
}

void Matrix_Wrapper::inverse(void)
{
  REQUIRE(this->is_square_matrix(), "invertable matrix should be square matrix");
  REQUIRE(!this->is_transposed_, "matrix should not be transposed");

  ms::math::blas::invA(this->data_ptr_, this->num_rows_, this->num_columns_, this->leading_dimension_);
}

void Matrix_Wrapper::transpose(void)
{
  this->is_transposed_ = !this->is_transposed_;
}

Matrix::Matrix(const int matrix_order)
    : Matrix_Wrapper(matrix_order, matrix_order, this->values_.data()), values_(matrix_order * matrix_order)
{
  this->essential_code_for_constructor();

  for (int i = 0; i < matrix_order; ++i)
  {
    this->at(i, i) = 1.0;
  }
}

Matrix::Matrix(const int matrix_order, const std::vector<double>& diagonal_values)
    : Matrix_Wrapper(matrix_order, matrix_order, this->values_.data()), values_(matrix_order * matrix_order)
{
  REQUIRE(matrix_order == diagonal_values.size(), "number of diagonal value of square matrix should be same with matrix order");

  this->essential_code_for_constructor();

  for (int i = 0; i < matrix_order; ++i)
  {
    this->at(i, i) = diagonal_values[i];
  }
}

Matrix::Matrix(const int num_row, const int num_column, const Transpose_Type transpose_type)
    : Matrix_Wrapper(num_row, num_column, this->values_.data(), transpose_type), values_(num_row * num_column)
{
  this->essential_code_for_constructor();
}

Matrix::Matrix(const int num_row, const int num_column, const double* value_ptr)
    : Matrix_Wrapper(num_row, num_column, this->values_.data()), values_(value_ptr, value_ptr + num_row * num_column)
{
  this->essential_code_for_constructor();
}

Matrix::Matrix(const int num_row, const int num_column, std::vector<double>&& values, const Transpose_Type transpose_type)
    : Matrix_Wrapper(num_row, num_column, this->values_.data(), transpose_type), values_(std::move(values))
{
  REQUIRE(num_row * num_column == this->values_.size(), "num value should be same with matrix size");

  this->essential_code_for_constructor();
}

Matrix::Matrix(const Matrix& other)
    : Matrix_Wrapper(other.num_rows_, other.num_columns_, this->values_.data(), other.leading_dimension_, other.transpose_type()), values_(other.values_)
{
  this->essential_code_for_constructor();
}

Matrix::Matrix(Matrix&& other) noexcept
    : Matrix_Wrapper(other.num_rows_, other.num_columns_, this->values_.data(), other.leading_dimension_, other.transpose_type()), values_(std::move(other.values_))
{
  this->essential_code_for_constructor();

  other.const_data_ptr_ = nullptr;
  other.data_ptr_       = nullptr;
}

void Matrix::operator=(const Matrix& other)
{
  this->is_transposed_  = other.is_transposed_;
  this->num_rows_       = other.num_rows_;
  this->num_columns_    = other.num_columns_;
  this->values_         = other.values_;
  this->const_data_ptr_ = this->values_.data();
  this->data_ptr_       = this->values_.data();
}

void Matrix::operator=(Matrix&& other) noexcept
{
  this->is_transposed_  = other.is_transposed_;
  this->num_rows_       = other.num_rows_;
  this->num_columns_    = other.num_columns_;
  this->values_         = std::move(other.values_);
  this->const_data_ptr_ = this->values_.data();
  this->data_ptr_       = this->values_.data();

  other.const_data_ptr_ = nullptr;
  other.data_ptr_       = nullptr;
}

void Matrix::essential_code_for_constructor(void)
{
  // This is code to prevent "data_ptr_" related member variables from becoming a dangling pointer.
  // In the constructor function, "values_" data() function call occurs before "values_" initialization.
  // "values_" reallocation may occur during initialization
  // Thus, "data_ptr_" related member variables reassigned here

  this->const_data_ptr_ = this->values_.data();
  this->data_ptr_       = this->values_.data();
}

Matrix operator*(const double constant, const Matrix_Const_Wrapper& M)
{
  return M * constant;
}

std::ostream& operator<<(std::ostream& os, const Matrix_Const_Wrapper& m)
{
  return os << m.to_string();
}

std::ostream& operator<<(std::ostream& os, const Matrix& m)
{
  return os << m.to_string();
}

} // namespace ms::math
