#include "msmath/Vector.h"

#include "msexception/Exception.h"
#include "msmath/BLAS.h"
#include <cmath>
#include <iomanip>
#include <sstream>

namespace
{
bool compare_double(const double d1, const double d2, const int ULP_precision = 4)
{
  const auto lower_ULP = d1 - std::nextafter(d1, std::numeric_limits<double>::lowest());
  const auto upper_ULP = std::nextafter(d1, std::numeric_limits<double>::max()) - d1;

  return d1 - ULP_precision * lower_ULP <= d2 && d2 <= d1 + ULP_precision * upper_ULP;
}
} // namespace

/*










*/

namespace ms::math
{

Vector_Const_Iterator& Vector_Const_Iterator::operator++(void)
{
  REQUIRE(this->_const_data_ptr != this->_last_data_ptr, "can't increment iterator past end");
  this->_const_data_ptr += _inc;
  return *this;
}

Vector_Const_Iterator Vector_Const_Iterator::operator++(int)
{
  Vector_Const_Iterator tmp = *this;
  ++(*this);
  return tmp;
}

Vector_Const_Iterator::reference Vector_Const_Iterator::operator*(void) const
{
  REQUIRE(this->_const_data_ptr != this->_last_data_ptr, "can't dereference end iterator");
  return *this->_const_data_ptr;
}

Vector_Const_Iterator::pointer Vector_Const_Iterator::operator->(void) const
{
  REQUIRE(this->_const_data_ptr != this->_last_data_ptr, "can't reference end iterator");
  return this->_const_data_ptr;
}

bool Vector_Const_Iterator::operator==(const Vector_Const_Iterator other) const
{
  return this->_const_data_ptr == other._const_data_ptr && this->_inc == other._inc && this->_last_data_ptr == other._last_data_ptr;
}

bool Vector_Const_Iterator::operator!=(const Vector_Const_Iterator other) const
{
  return this->_const_data_ptr != other._const_data_ptr;
}

/*










*/

Vector<0> Vector_View::operator*(const double scalar) const
{
  constexpr auto incr = 1;

  const auto dim  = this->dimension();
  const auto data = this->data();
  const auto inc  = this->inc();

  Vector result(dim);

  ms::math::blas::cx(result.data(), scalar, data, dim, incr, inc);

  return result;
}

Vector<0> Vector_View::operator+(const Vector_View other) const
{
  constexpr auto incr = 1;

  const auto dim        = this->dimension();
  const auto data       = this->data();
  const auto inc        = this->inc();
  const auto other_dim  = other.dimension();
  const auto other_data = other._values_view.data();
  const auto other_inc  = other.inc();
  REQUIRE(dim == other_dim, "other vector should be same dimension");

  Vector result(dim);
  ms::math::blas::x_plus_y(result.data(), data, other_data, dim, incr, inc, other_inc);

  return result;
}

Vector<0> Vector_View::operator-(const Vector_View other) const
{
  constexpr auto incr = 1;

  const auto dim        = this->dimension();
  const auto data       = this->data();
  const auto inc        = this->inc();
  const auto other_dim  = other.dimension();
  const auto other_data = other._values_view.data();
  const auto other_inc  = other.inc();
  REQUIRE(dim == other_dim, "other vector should be same dimension");

  Vector result(dim);
  ms::math::blas::x_minus_y(result.data(), data, other_data, dim, incr, inc, other_inc);

  return result;
}

double Vector_View::operator[](const int position) const
{
  REQUIRE(0 <= position, "position should be positive number");
  REQUIRE(position < this->_dimension, "position should be less then size");
  return this->_values_view[position * this->_inc];
}

bool Vector_View::operator==(const Vector_View other) const
{
  if (this->_dimension != other._dimension)
    return false;

  for (int i = 0; i < this->_dimension; ++i)
  {
    if (this->_values_view[i * this->_inc] != other._values_view[i * other._inc])
      return false;
  }

  return true;
}

bool Vector_View::operator!=(const Vector_View other) const
{
  return !((*this) == other);
}

double Vector_View::at(const int position) const
{
  REQUIRE(position < this->_dimension, "position can't exceed given range");
  return this->_values_view[position * this->_inc];
}

Vector_Const_Iterator Vector_View::begin(void) const
{
  Vector_Const_Iterator begin(this->data(), this->_inc, this->end_ptr());
  return begin;
}

Vector_Const_Iterator Vector_View::end(void) const
{
  const auto end_ptr = this->end_ptr();

  Vector_Const_Iterator end(end_ptr, this->_inc, end_ptr);
  return end;
}

double Vector_View::cosine(const Vector_View other) const
{
  const auto dim       = this->dimension();
  const auto other_dim = other.dimension();
  REQUIRE(dim == other_dim, "other vector should be same dimension");

  return this->inner_product(other) / (this->L2_norm() * other.L2_norm());
}

inline int Vector_View::dimension(void) const
{
  return this->_dimension;
}
const double* Vector_View::data(void) const
{
  return this->_values_view.data();
}

bool Vector_View::empty(void) const
{
  return this->_values_view.empty();
}

Vector_View Vector_View::get_view(void) const
{
  return *this;
}

inline int Vector_View::inc(void) const
{
  return this->_inc;
}

double Vector_View::inner_product(const Vector_View other) const
{
  const auto dim        = this->dimension();
  const auto data       = this->data();
  const auto inc        = this->inc();
  const auto other_dim  = other.dimension();
  const auto other_data = other._values_view.data();
  const auto other_inc  = other._inc;
  REQUIRE(dim == other_dim, "other vector should be same dimension");

  return ms::math::blas::x_dot_y(data, other_data, dim, inc, other_inc);
}

bool Vector_View::is_parallel(const Vector_View other) const
{
  const auto abs_cosine = std::abs(this->cosine(other));
  return compare_double(abs_cosine, 1.0);
}

double Vector_View::L1_norm(void) const
{
  const auto dim  = this->dimension();
  const auto data = this->data();
  const auto inc  = this->inc();
  return ms::math::blas::abs_sum_x(data, dim, inc);
}

inline double Vector_View::L2_norm(void) const
{
  return std::sqrt(this->inner_product(*this));
}

double Vector_View::Linf_norm(void) const
{
  const auto dim  = this->dimension();
  const auto data = this->data();
  const auto inc  = this->inc();

  const auto pos = static_cast<int>(ms::math::blas::find_maximum_element_pos(data, dim, inc));
  return this->_values_view[pos];
}

Vector_View Vector_View::sub_view(const int start_pos, const int end_pos) const
{
  REQUIRE(start_pos < end_pos, "start position should be smaller than end position");
  REQUIRE(end_pos <= this->_dimension, "end position should be smaller or equal to dimension");

  const auto num_values = end_pos - start_pos - 1; // [start_pos, end_pos)

  const auto start_index     = start_pos * this->_inc;
  const auto value_count     = 1 + this->_inc * num_values;
  const auto sub_values_view = this->_values_view.subspan(start_index, value_count);

  const auto result = Vector_View(sub_values_view, this->_inc);
  return result;
}

Vector_View Vector_View::sub_view(const int start_position, const int inc, const int num_values) const
{
  REQUIRE(0 <= start_position, "start position should not be negative number");
  REQUIRE(0 < inc, "increment should be positive number");
  REQUIRE(0 < num_values, "number of values should be positive number");

  const auto new_inc = this->_inc * inc;

  const auto start_index     = start_position * this->_inc;
  const auto value_count     = 1 + new_inc * (num_values - 1);
  const auto sub_values_view = this->_values_view.subspan(start_index, value_count);

  const auto result = Vector_View(sub_values_view, new_inc);
  return result;
}

std::string Vector_View::to_string(void) const
{
  std::ostringstream oss;
  oss << std::setprecision(16) << std::showpoint << std::left;
  for (int i = 0; i < this->_dimension; ++i)
  {
    oss << std::setw(25) << this->at(i);
  }
  return oss.str();
}

std::vector<double> Vector_View::to_vector(void) const
{
  std::vector<double> result(this->_dimension);
  ms::math::blas::copy(result.data(), this->data(), this->_dimension, 1, this->_inc);
  return result;
}

template <>
Vector<3> Vector_View::cross_product<2>(const Vector_View other) const
{
  Vector<3> result;
  result[2] = (*this)[0] * other[1] - (*this)[1] * other[0];

  return result;
}

template <>
Vector<3> Vector_View::cross_product<3>(const Vector_View other) const
{
  Vector<3> result;
  result[0] = (*this)[1] * other[2] - (*this)[2] * other[1];
  result[1] = (*this)[2] * other[0] - (*this)[0] * other[2];
  result[2] = (*this)[0] * other[1] - (*this)[1] * other[0];

  return result;
}

void Vector_View::initialize(void)
{
  REQUIRE(0 < this->_inc, "inc should greater than 0");

  this->_dimension = static_cast<int>(std::ceil(this->_values_view.size() / static_cast<double>(this->_inc)));
  REQUIRE(0 < this->_dimension, "number of values should greater than 0");
}

inline const double* Vector_View::end_ptr(void) const
{
  const auto last_index = static_cast<int>(std::ceil(this->_values_view.size() / static_cast<double>(this->_inc)) * this->_inc);
  return this->data() + last_index;
}

/*










*/

void Vector_Wrap::operator*=(const double constant)
{
  const auto dim  = this->dimension();
  auto       data = this->data();
  const auto inc  = this->inc();

  ms::math::blas::cx(constant, data, dim, inc);
}

void Vector_Wrap::operator+=(const Vector_View other)
{
  const auto dim        = this->dimension();
  auto       data       = this->data();
  const auto inc        = this->inc();
  const auto other_dim  = other.dimension();
  const auto other_data = other.data();
  const auto other_inc  = other.inc();

  REQUIRE(dim == other.dimension(), "other vector should be same size");
  ms::math::blas::x_plus_assign_y(data, other.data(), dim, inc, other.inc());
}

void Vector_Wrap::operator-=(const Vector_View other)
{
  const auto dim        = this->dimension();
  auto       data       = this->data();
  const auto inc        = this->inc();
  const auto other_dim  = other.dimension();
  const auto other_data = other.data();
  const auto other_inc  = other.inc();

  REQUIRE(dim == other.dimension(), "other vector should be same size");
  ms::math::blas::x_minus_assign_y(data, other.data(), dim, inc, other.inc());
}

double& Vector_Wrap::operator[](const int index)
{
  REQUIRE(0 <= index, "index should be positive");
  REQUIRE(index < this->_dimension, "index should be less then size");
  return this->_values_wrap[index * this->_inc];
}

double& Vector_Wrap::at(const int index)
{
  REQUIRE(0 <= index, "index should be positive");
  REQUIRE(index < this->_dimension, "index should be less then size");
  return this->_values_wrap[index * this->_inc];
}

void Vector_Wrap::change_value(const Vector_View other)
{
  REQUIRE(this->_dimension == other.dimension(), "should be same size");
  ms::math::blas::copy(this->data(), other.data(), this->_dimension, this->_inc, other.inc());
}

double* Vector_Wrap::data(void)
{
  return this->_values_wrap.data();
}

void Vector_Wrap::normalize(void)
{
  *this *= 1.0 / this->L2_norm();
}

Vector_Wrap Vector_Wrap::sub_wrap(const int start_position, const int end_position) const
{
  REQUIRE(start_position < end_position, "start position should be smaller than end position");
  REQUIRE(end_position <= this->_dimension, "end position should be smaller or equal to dimension");

  const auto num_values = end_position - start_position - 1; // [start_position, end_pos)

  const auto start_index     = start_position * this->_inc;
  const auto value_count     = 1 + this->_inc * num_values;
  const auto sub_values_wrap = this->_values_wrap.subspan(start_index, value_count);

  auto result = Vector_Wrap(sub_values_wrap, this->_inc);
  return result;
}

Vector_Wrap Vector_Wrap::sub_wrap(const int start_position, const int inc, const int num_values) const
{
  REQUIRE(0 <= start_position, "start position should not be negative number");
  REQUIRE(0 < inc, "increment should be positive number");
  REQUIRE(0 < num_values, "number of values should be positive number");

  const auto new_inc = this->_inc * inc;

  const auto start_index     = start_position * this->_inc;
  const auto value_count     = 1 + new_inc * (num_values - 1);
  const auto sub_values_wrap = this->_values_wrap.subspan(start_index, value_count);

  auto result = Vector_Wrap(sub_values_wrap, new_inc);
  return result;
}

/*










*/

Vector<0>::Vector(const int dimension)
    : _coordinates(dimension)
{
  this->reallocate_values();
}

Vector<0>::Vector(const std::initializer_list<double> list)
    : _coordinates(list)
{
  this->reallocate_values();
}

Vector<0>::Vector(const Vector_View vector_view)
    : _coordinates(vector_view.to_vector())
{
  this->reallocate_values();
}

Vector<0>::Vector(const Vector& other)
    : _coordinates(other._coordinates)
{
  this->reallocate_values();
}

Vector<0>::Vector(Vector&& other) noexcept
    : _coordinates(std::move(other._coordinates))
{
  this->reallocate_values();
}

void Vector<0>::operator=(const Vector& other)
{
  this->_coordinates = other._coordinates;
  this->reallocate_values();
}

void Vector<0>::operator=(Vector&& other) noexcept
{
  this->_coordinates = std::move(other._coordinates);
  this->reallocate_values();
}

double* Vector<0>::data(void)
{
  return this->_coordinates.data();
}

void Vector<0>::resize(const int size)
{
  REQUIRE(0 < size, "size should be positive");

  if (size <= this->dimension()) return;

  this->_coordinates.resize(size);
  this->reallocate_values();
}

void Vector<0>::reallocate_values(void)
{
  this->_values_view = this->_coordinates;
  this->_dimension   = static_cast<int>(this->_coordinates.size());
  this->_inc         = 1;
  this->_values_wrap = this->_coordinates;
}

/*










*/

Vector<0> operator*(const double constant, const Vector_View x)
{
  return x * constant;
}

} // namespace ms::math



