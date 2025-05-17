#include "msmath/Matrix.h"

#include "msexception/Exception.h"
#include "msmath/BLAS.h"
#include <iomanip>
#include <sstream>

namespace ms::math
{

  Matrix Matrix_View::operator+( const Matrix_View other ) const
  {
    REQUIRE( this->size() == other.size(), "two matrix should be same size" );

    const auto result_num_rows    = static_cast<int>( this->num_rows() );
    const auto result_num_columns = static_cast<int>( this->num_columns() );
    Matrix     result( result_num_rows, result_num_columns );

    constexpr double c = 1.0;

    const auto this_data  = this->data();
    const auto other_data = other.data();

    if ( this->_is_transposed )
    {
      if ( other._is_transposed )
      {
        ms::math::blas::AT_plus_cBT( result.data(), this_data, c, other_data, result_num_rows, result_num_columns, result_num_columns, this->_leading_dimension, other._leading_dimension );
      }
      else
      {
        ms::math::blas::AT_plus_cB( result.data(), this_data, c, other_data, result_num_rows, result_num_columns, result_num_columns, this->_leading_dimension, other._leading_dimension );
      }
    }
    else
    {
      if ( other._is_transposed )
      {
        ms::math::blas::A_plus_cBT( result.data(), this_data, c, other_data, result_num_rows, result_num_columns, result_num_columns, this->_leading_dimension, other._leading_dimension );
      }
      else
      {
        ms::math::blas::A_plus_cB( result.data(), this_data, c, other_data, result_num_rows, result_num_columns, result_num_columns, this->_leading_dimension, other._leading_dimension );
      }
    }

    return result;
  }

  Matrix Matrix_View::operator-( const Matrix_View other ) const
  {
    REQUIRE( this->size() == other.size(), "two matrix should be same size" );

    const auto result_num_rows    = static_cast<int>( this->num_rows() );
    const auto result_num_columns = static_cast<int>( this->num_columns() );
    Matrix     result( result_num_rows, result_num_columns );

    constexpr double c = -1.0;

    const auto this_data  = this->data();
    const auto other_data = other.data();

    if ( this->_is_transposed )
    {
      if ( other._is_transposed )
      {
        ms::math::blas::AT_plus_cBT( result.data(), this_data, c, other_data, result_num_rows, result_num_columns, result_num_columns, this->_leading_dimension, other._leading_dimension );
      }
      else
      {
        ms::math::blas::AT_plus_cB( result.data(), this_data, c, other_data, result_num_rows, result_num_columns, result_num_columns, this->_leading_dimension, other._leading_dimension );
      }
    }
    else
    {
      if ( other._is_transposed )
      {
        ms::math::blas::A_plus_cBT( result.data(), this_data, c, other_data, result_num_rows, result_num_columns, result_num_columns, this->_leading_dimension, other._leading_dimension );
      }
      else
      {
        ms::math::blas::A_plus_cB( result.data(), this_data, c, other_data, result_num_rows, result_num_columns, result_num_columns, this->_leading_dimension, other._leading_dimension );
      }
    }

    return result;
  }

  Matrix Matrix_View::operator*( const double constant ) const
  {
    Matrix result( this->_num_rows, this->_num_columns );

    ms::math::blas::cA( result.data(), constant, this->_values_view.data(), this->_num_rows, this->_num_columns, this->_leading_dimension );

    if ( this->_is_transposed )
    {
      result.be_transpose();
    }

    return result;
  }

  Matrix Matrix_View::operator*( const Matrix_View other ) const
  {
    REQUIRE( this->num_columns() == other.num_rows(), "two matrix size should be mathched" );

    const auto result_num_rows    = static_cast<int>( this->num_rows() );
    const auto result_num_columns = static_cast<int>( other.num_columns() );
    Matrix     result( result_num_rows, result_num_columns );

    constexpr double c = 1.0;

    const auto this_data  = this->data();
    const auto other_data = other.data();

    if ( this->_is_transposed )
    {
      if ( other._is_transposed )
      {
        ms::math::blas::cATBT( result.data(), c, this_data, other_data, this->_num_rows, this->_num_columns, other._num_rows, this->_leading_dimension, other._leading_dimension );
      }
      else
      {
        ms::math::blas::cATB( result.data(), c, this_data, other_data, this->_num_rows, this->_num_columns, other._num_columns, this->_leading_dimension, other._leading_dimension );
      }
    }
    else
    {
      if ( other._is_transposed )
      {
        ms::math::blas::cABT( result.data(), c, this_data, other_data, this->_num_rows, this->_num_columns, other._num_rows, this->_leading_dimension, other._leading_dimension );
      }
      else
      {
        ms::math::blas::cAB( result.data(), c, this_data, other_data, this->_num_rows, this->_num_columns, other._num_columns, this->_leading_dimension, other._leading_dimension );
      }
    }

    return result;
  }

  Vector<0> Matrix_View::operator*( const Vector_View vec ) const
  {
    const auto leading_dimension = this->leading_dimension();
    const auto data_ptr          = this->data();

    if ( this->_is_transposed )
    {
      REQUIRE( this->_num_rows == vec.dimension(), "size should be mathced" );

      Vector result( this->_num_columns );
      ms::math::blas::ATx( result.data(), data_ptr, vec.data(), this->_num_rows, this->_num_columns, leading_dimension, result.inc(), vec.inc() );

      return result;
    }
    else
    {
      REQUIRE( this->_num_columns == vec.dimension(), "size should be mathced" );

      Vector result( this->_num_rows );
      ms::math::blas::Ax( result.data(), data_ptr, vec.data(), this->_num_rows, this->_num_columns, leading_dimension, result.inc(), vec.inc() );

      return result;
    }
  }

  bool Matrix_View::operator==( const Matrix_View other ) const
  {
    if ( this->size() != other.size() )
    {
      return false;
    }

    const auto [num_rows, num_columns] = this->size();

    for ( int i = 0; i < num_rows; i++ )
    {
      for ( int j = 0; j < num_columns; j++ )
      {
        if ( this->at( i, j ) != other.at( i, j ) )
          return false;
      }
    }

    return true;
  }

  double Matrix_View::at( const int row_index, const int column_index ) const
  {
    const auto pos = this->index_to_position( row_index, column_index );
    return this->_values_view[pos];
  }

  Vector_View Matrix_View::column_vector_view( const int index ) const
  {
    const auto values_part_view = this->values_part_view( 0, this->num_rows(), index, index + 1 );

    if ( this->_is_transposed )
    {
      Vector_View column_vector( values_part_view );
      return column_vector;
    }
    else
    {
      Vector_View column_vector( values_part_view, this->_leading_dimension );
      return column_vector;
    }
  }

  inline const double* Matrix_View::data( void ) const
  {
    return this->_values_view.data();
  }

  inline bool Matrix_View::has_compact_data( void ) const
  {
    return this->_leading_dimension == this->_num_columns;
  }

  inline bool Matrix_View::is_transposed( void ) const
  {
    return this->_is_transposed;
  }

  Matrix Matrix_View::inverse( void ) const
  {
    REQUIRE( this->is_square_matrix(), "matrix should be square matrix" );
    REQUIRE( !this->_is_transposed, "matrix should not be transposed" );

    Matrix result( this->_num_rows, this->_num_columns, this->values_vector() );
    result.be_inverse();
    return result;
  }

  inline int Matrix_View::leading_dimension( void ) const
  {
    return this->_leading_dimension;
  }

  inline int Matrix_View::num_columns( void ) const
  {
    if ( this->_is_transposed )
    {
      return this->_num_rows;
    }
    else
    {
      return this->_num_columns;
    }
  }

  inline int Matrix_View::num_rows( void ) const
  {
    if ( this->_is_transposed )
    {
      return this->_num_columns;
    }
    else
    {
      return this->_num_rows;
    }
  }

  inline int Matrix_View::num_values( void ) const
  {
    return this->_num_rows * this->_num_columns;
  }

  Matrix_View Matrix_View::sub_view( const int start_row_index, const int end_row_index, const int start_column_index, const int end_column_index ) const
  {
    REQUIRE( this->is_valid_row_range( start_row_index, end_row_index ), "Given row range should be valid" );
    REQUIRE( this->is_valid_column_range( start_column_index, end_column_index ), "Given column range should be valid" );

    const auto values_part_view = this->values_part_view( start_row_index, end_row_index, start_column_index, end_column_index );

    const auto num_row    = end_row_index - start_row_index;
    const auto num_column = end_column_index - start_column_index;

    if ( this->_is_transposed )
    {
      Matrix_View part( num_column, num_row, values_part_view, this->_leading_dimension );
      part._is_transposed = true;
      return part;
    }
    else
    {
      Matrix_View part( num_row, num_column, values_part_view, this->_leading_dimension );
      return part;
    }
  }

  Vector_View Matrix_View::row_vector_view( const int index ) const
  {
    const auto values_part_view = this->values_part_view( index, index + 1, 0, this->num_columns() );

    if ( this->_is_transposed )
    {
      Vector_View row_vector( values_part_view, this->_leading_dimension );
      return row_vector;
    }
    else
    {
      Vector_View row_vector( values_part_view );
      return row_vector;
    }
  }

  std::pair<int, int> Matrix_View::size( void ) const
  {
    if ( this->_is_transposed )
    {
      return { this->_num_columns, this->_num_rows };
    }
    else
    {
      return { this->_num_rows, this->_num_columns };
    }
  }

  Matrix_View Matrix_View::transpose( void ) const
  {
    auto result           = *this;
    result._is_transposed = !this->_is_transposed;

    return result;
  }

  std::string Matrix_View::to_string( void ) const
  {
    std::ostringstream oss;
    oss << std::setprecision( 16 ) << std::showpoint << std::left;
    for ( int i = 0; i < this->num_rows(); ++i )
    {
      for ( int j = 0; j < this->num_columns(); ++j )
        oss << std::setw( 25 ) << this->at( i, j );
      oss << "\n";
    }
    return oss.str();
  }

  Matrix_View Matrix_View::view( void ) const
  {
    return *this;
  }

  bool Matrix_View::is_square_matrix( void ) const
  {
    return this->_num_rows == this->_num_columns;
  }

  bool Matrix_View::is_valid_column_index( const int index ) const
  {
    if ( index < 0 )
    {
      return false;
    }

    if ( this->_is_transposed )
    {
      return index < this->_num_rows;
    }
    else
    {
      return index < this->_num_columns;
    }
  }

  bool Matrix_View::is_valid_column_range( const int start_index, const int end_index ) const
  {
    if ( end_index <= start_index )
    {
      return false;
    }

    if ( start_index < 0 )
    {
      return false;
    }

    if ( this->_is_transposed )
    {
      return end_index <= this->_num_rows;
    }
    else
    {
      return end_index <= this->_num_columns;
    }
  }

  bool Matrix_View::is_valid_row_index( const int index ) const
  {
    if ( index < 0 )
    {
      return false;
    }

    if ( this->_is_transposed )
    {
      return index < this->_num_columns;
    }
    else
    {
      return index < this->_num_rows;
    }
  }

  bool Matrix_View::is_valid_row_range( const int start_index, const int end_index ) const
  {
    if ( end_index <= start_index )
    {
      return false;
    }

    if ( start_index < 0 )
    {
      return false;
    }

    if ( this->_is_transposed )
    {
      return end_index <= this->_num_columns;
    }
    else
    {
      return end_index <= this->_num_rows;
    }
  }

  int Matrix_View::index_to_position( const int row_index, const int column_index ) const
  {
    REQUIRE( this->is_valid_row_index( row_index ), "Given row index should be valid" );
    REQUIRE( this->is_valid_column_index( column_index ), "Given column index should be valid" );

    if ( this->_is_transposed )
    {
      return column_index * this->_leading_dimension + row_index;
    }
    else
    {
      return row_index * this->_leading_dimension + column_index;
    }
  }

  std::span<const double> Matrix_View::values_part_view( const int start_row_index, const int end_row_index, const int start_column_index, const int end_column_index ) const
  {
    const auto start_pos  = this->index_to_position( start_row_index, start_column_index );
    const auto end_pos    = this->index_to_position( end_row_index - 1, end_column_index - 1 );
    const auto num_values = end_pos - start_pos + 1;
    return this->_values_view.subspan( start_pos, num_values );
  }

  std::vector<double> Matrix_View::values_vector( void ) const
  {
    const auto          n = static_cast<int>( this->num_values() );
    std::vector<double> result( n );

    if ( this->has_compact_data() )
    {
      ms::math::blas::copy( result.data(), this->_values_view.data(), n );
    }
    else
    {
      auto result_ptr = result.data();
      auto value_ptr  = this->_values_view.data();
      for ( int i = 0; i < this->_num_rows; i++ )
      {
        ms::math::blas::copy( result_ptr, value_ptr, this->_num_columns );
        result_ptr += this->_num_columns;
        value_ptr  += this->_leading_dimension;
      }
    }

    return result;
  }

} // namespace ms::math

/*










*/

namespace ms::math
{

  void Matrix_Wrap::operator*=( const double constant )
  {
    ms::math::blas::cA( constant, this->data(), this->_num_rows, this->_num_columns, this->_leading_dimension );
  }

  void Matrix_Wrap::operator+=( const Matrix_View other )
  {
    REQUIRE( this->size() == other.size(), "two matrix should be same size" );
    REQUIRE( !this->_is_transposed, "matrix should not be transposed" );

    constexpr double c = 1.0;

    const auto this_data = this->data();

    if ( other.is_transposed() )
    {
      ms::math::blas::A_plus_cBT( this_data, this_data, c, other.data(), this->_num_rows, this->_num_columns, this->_leading_dimension, this->_leading_dimension, other.leading_dimension() );
    }
    else
    {
      ms::math::blas::A_plus_cB( this_data, this_data, c, other.data(), this->_num_rows, this->_num_columns, this->_leading_dimension, this->_leading_dimension, other.leading_dimension() );
    }
  }

  void Matrix_Wrap::operator-=( const Matrix_View other )
  {
    REQUIRE( this->size() == other.size(), "two matrix should be same size" );
    REQUIRE( !this->_is_transposed, "matrix should not be transposed" );

    constexpr double c = -1.0;

    const auto this_data = this->data();

    if ( other.is_transposed() )
    {
      ms::math::blas::A_plus_cBT( this_data, this_data, c, other.data(), this->_num_rows, this->_num_columns, this->_leading_dimension, this->_leading_dimension, other.leading_dimension() );
    }
    else
    {
      ms::math::blas::A_plus_cB( this_data, this_data, c, other.data(), this->_num_rows, this->_num_columns, this->_leading_dimension, this->_leading_dimension, other.leading_dimension() );
    }
  }

  double& Matrix_Wrap::at( const int row_index, const int column_index )
  {
    REQUIRE( this->is_valid_row_index( row_index ) && this->is_valid_column_index( column_index ), "Given numbers should be valid" );

    if ( this->_is_transposed )
    {
      return this->_values_wrap[column_index * this->_leading_dimension + row_index];
    }
    else
    {
      return this->_values_wrap[row_index * this->_leading_dimension + column_index];
    }
  }

  void Matrix_Wrap::be_inverse( void )
  {
    REQUIRE( this->is_square_matrix(), "invertable matrix should be square matrix" );
    REQUIRE( !this->_is_transposed, "matrix should not be transposed" );

    ms::math::blas::invA( this->data(), this->_num_rows, this->_num_columns, this->_leading_dimension );
  }

  void Matrix_Wrap::be_transpose( void )
  {
    this->_is_transposed = !this->_is_transposed;
  }

  void Matrix_Wrap::change_column( const int index, const double* values )
  {
    const auto num_cols = this->num_columns();
    auto       column   = this->column_vector_wrap( index );

    for ( int i = 0; i < num_cols; ++i )
    {
      column[i] = values[i];
    }
  }

  void Matrix_Wrap::change_row( const int index, const double* values )
  {
    const auto num_rows = this->num_rows();
    auto       row      = this->row_vector_wrap( index );

    for ( int i = 0; i < num_rows; ++i )
    {
      row[i] = values[i];
    }
  }

  void Matrix_Wrap::change_value( const Matrix_View other )
  {
    REQUIRE( this->size() == other.size(), "size should be matched" );

    for ( int i = 0; i < this->num_rows(); i++ )
    {
      for ( int j = 0; j < this->num_columns(); j++ )
      {
        this->at( i, j ) = other.at( i, j );
      }
    }
  }

  Vector_Wrap Matrix_Wrap::column_vector_wrap( const int index )
  {
    const auto values_part_wrap = this->values_part_wrap( 0, this->num_rows(), index, index + 1 );

    if ( this->_is_transposed )
    {
      Vector_Wrap column_vector( values_part_wrap );
      return column_vector;
    }
    else
    {
      Vector_Wrap column_vector( values_part_wrap, this->_leading_dimension );
      return column_vector;
    }

    REQUIRE( this->is_valid_column_index( index ), "Given index should be valid" );

    if ( this->_is_transposed )
    {
      const auto start_index      = index * this->_leading_dimension;
      const auto num_values       = this->_num_columns;
      const auto values_wrap_part = this->_values_wrap.subspan( start_index, num_values );

      Vector_Wrap column_vector( values_wrap_part );
      return column_vector;
    }
    else
    {
      const auto start_index      = index;
      const auto num_values       = this->_num_rows * this->_leading_dimension;
      const auto values_wrap_part = this->_values_wrap.subspan( start_index, num_values );

      Vector_Wrap column_vector( values_wrap_part, this->_leading_dimension );
      return column_vector;
    }
  }

  double* Matrix_Wrap::data( void )
  {
    return this->_values_wrap.data();
  }

  Vector_Wrap Matrix_Wrap::row_vector_wrap( const int index )
  {
    const auto values_wrap_part = this->values_part_wrap( index, index + 1, 0, this->num_columns() );

    REQUIRE( this->is_valid_row_index( index ), "Given index should be valid" );

    if ( this->_is_transposed )
    {
      Vector_Wrap row_vector( values_wrap_part, this->_leading_dimension );
      return row_vector;
    }
    else
    {
      Vector_Wrap row_vector( values_wrap_part );
      return row_vector;
    }
  }

  Matrix_Wrap Matrix_Wrap::sub_wrap( const int start_row_index, const int end_row_index, const int start_column_index, const int end_column_index )
  {
    REQUIRE( this->is_valid_row_range( start_row_index, end_row_index ), "Given row range should be valid" );
    REQUIRE( this->is_valid_column_range( start_column_index, end_column_index ), "Given column range should be valid" );

    const auto values_part_wrap = this->values_part_wrap( start_row_index, end_row_index, start_column_index, end_column_index );

    const auto num_row    = end_row_index - start_row_index;
    const auto num_column = end_column_index - start_column_index;

    if ( this->_is_transposed )
    {
      Matrix_Wrap result( num_column, num_row, values_part_wrap, this->_leading_dimension );
      result._is_transposed = true;
      return result;
    }
    else
    {
      Matrix_Wrap result( num_row, num_column, values_part_wrap, this->_leading_dimension );
      return result;
    }
  }

  Matrix_Wrap Matrix_Wrap::wrap( void )
  {
    return *this;
  }

  std::span<double> Matrix_Wrap::values_part_wrap( const int start_row_index, const int end_row_index, const int start_column_index, const int end_column_index )
  {
    const auto start_pos  = this->index_to_position( start_row_index, start_column_index );
    const auto end_pos    = this->index_to_position( end_row_index - 1, end_column_index - 1 );
    const auto num_values = end_pos - start_pos + 1;
    return this->_values_wrap.subspan( start_pos, num_values );
  }

} // namespace ms::math

/*










*/

namespace ms::math
{

  Matrix::Matrix( const int matrix_order )
    : _values( matrix_order * matrix_order )
  {
    this->_num_rows          = matrix_order;
    this->_num_columns       = matrix_order;
    this->_leading_dimension = this->_num_columns;
    this->reallocate_values();

    for ( int i = 0; i < matrix_order; ++i )
    {
      this->at( i, i ) = 1.0;
    }
  }

  Matrix::Matrix( const int matrix_order, const std::vector<double>& diagonal_values )
    : _values( matrix_order * matrix_order )
  {
    REQUIRE( matrix_order == diagonal_values.size(), "number of diagonal value of square matrix should be same with matrix order" );

    this->_num_rows          = matrix_order;
    this->_num_columns       = matrix_order;
    this->_leading_dimension = this->_num_columns;
    this->reallocate_values();

    for ( int i = 0; i < matrix_order; ++i )
    {
      this->at( i, i ) = diagonal_values[i];
    }
  }

  Matrix::Matrix( const int num_rows, const int num_columns )
    : _values( num_rows * num_columns )
  {
    this->_num_rows          = num_rows;
    this->_num_columns       = num_columns;
    this->_leading_dimension = this->_num_columns;
    this->reallocate_values();
  }

  Matrix::Matrix( const int num_rows, const int num_columns, std::vector<double>&& values )
    : _values( std::move( values ) )
  {
    REQUIRE( num_rows * num_columns == this->_values.size(), "num value should be same with matrix size" );

    this->_num_rows          = num_rows;
    this->_num_columns       = num_columns;
    this->_leading_dimension = this->_num_columns;
    this->reallocate_values();
  }

  Matrix::Matrix( const Matrix& other )
    : _values( other._values )
  {
    this->_is_transposed     = other._is_transposed;
    this->_num_rows          = other._num_rows;
    this->_num_columns       = other._num_columns;
    this->_leading_dimension = this->_num_columns;
    this->reallocate_values();
  }

  Matrix::Matrix( Matrix&& other ) noexcept
    : _values( std::move( other._values ) )
  {
    this->_is_transposed     = other._is_transposed;
    this->_num_rows          = other._num_rows;
    this->_num_columns       = other._num_columns;
    this->_leading_dimension = this->_num_columns;
    this->reallocate_values();

    other = Matrix();
  }

  void Matrix::operator=( const Matrix& other )
  {
    this->_is_transposed     = other._is_transposed;
    this->_num_rows          = other._num_rows;
    this->_num_columns       = other._num_columns;
    this->_leading_dimension = this->_num_columns;
    this->_values            = other._values;
    this->reallocate_values();
  }

  void Matrix::operator=( Matrix&& other ) noexcept
  {
    this->_is_transposed     = other._is_transposed;
    this->_num_rows          = other._num_rows;
    this->_num_columns       = other._num_columns;
    this->_leading_dimension = this->_num_columns;
    this->_values            = std::move( other._values );
    this->reallocate_values();

    other._leading_dimension = 0;
    other._num_columns       = 0;
    other._num_rows          = 0;
    other._values.clear();
    other._values_view = std::span<const double>();
    other._values_wrap = std::span<double>();
  }

  void Matrix::resize( const int row, const int column )
  {
    this->_num_rows    = row;
    this->_num_columns = column;
    this->_values.resize( row, column );
    this->reallocate_values();
  }

  void Matrix::reallocate_values( void )
  {
    this->_values_view = this->_values;
    this->_values_wrap = this->_values;
  }

} // namespace ms::math

/*










*/

namespace ms::math
{

  Matrix operator*( const double constant, const Matrix_View M )
  {
    return M * constant;
  }

  std::ostream& operator<<( std::ostream& os, const Matrix_View m )
  {
    return os << m.to_string();
  }

} // namespace ms::math

namespace ms::math::blas
{
  void mpm( Matrix_Wrap result, const Matrix_View m1, const Matrix_View m2 )
  {
    REQUIRE( result.size() == m1.size(), "two matrix should be same size" );
    REQUIRE( m1.size() == m2.size(), "two matrix should be same size" );

    const auto num_rows   = static_cast<int>( m1.num_rows() );
    const auto num_cols   = static_cast<int>( m1.num_columns() );
    const auto result_ptr = result.data();
    const auto m1_ptr     = m1.data();
    const auto m2_ptr     = m2.data();
    const auto result_ld  = result.leading_dimension();
    const auto m1_ld      = m1.leading_dimension();
    const auto m2_ld      = m2.leading_dimension();

    constexpr double c = 1.0;

    if ( m1.is_transposed() )
    {
      if ( m2.is_transposed() )
      {
        ms::math::blas::AT_plus_cBT( result_ptr, m1_ptr, c, m2_ptr, num_rows, num_cols, result_ld, m1_ld, m2_ld );
      }
      else
      {
        ms::math::blas::AT_plus_cB( result_ptr, m1_ptr, c, m2_ptr, num_rows, num_cols, result_ld, m1_ld, m2_ld );
      }
    }
    else
    {
      if ( m2.is_transposed() )
      {
        ms::math::blas::A_plus_cBT( result_ptr, m1_ptr, c, m2_ptr, num_rows, num_cols, result_ld, m1_ld, m2_ld );
      }
      else
      {
        ms::math::blas::A_plus_cB( result_ptr, m1_ptr, c, m2_ptr, num_rows, num_cols, result_ld, m1_ld, m2_ld );
      }
    }
  }

  void cmv( Vector_Wrap result, const double c, const Matrix_View m, const Vector_View v )
  {
    const auto num_rows   = static_cast<int>( m.num_rows() );
    const auto num_cols   = static_cast<int>( m.num_columns() );
    const auto result_dim = result.dimension();
    const auto result_ptr = result.data();
    const auto result_inc = result.inc();
    const auto m_ptr      = m.data();
    const auto m_ld       = m.leading_dimension();
    const auto v_dim      = v.dimension();
    const auto v_ptr      = v.data();
    const auto v_inc      = v.inc();

    if ( m.is_transposed() )
    {
      REQUIRE( num_rows == v_dim, "size should be mathced" );
      REQUIRE( num_cols == result_dim, "size should be matched" );

      ms::math::blas::cATx( result_ptr, c, m_ptr, v_ptr, num_rows, num_cols, m_ld, result_inc, v_inc );
    }
    else
    {
      REQUIRE( num_cols == v_dim, "size should be mathced" );
      REQUIRE( num_rows == result_dim, "size should be matched" );

      ms::math::blas::cAx( result_ptr, c, m_ptr, v_ptr, num_rows, num_cols, m_ld, result_inc, v_inc );
    }
  }

} // namespace ms::math::blas
