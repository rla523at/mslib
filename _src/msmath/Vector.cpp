#include "Vector.h"

#include "BLAS.h"

#include "msexception/Exception.h"
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

Vector_Values_Const_Iterator& Vector_Values_Const_Iterator::operator++(void)
{
  REQUIRE(this->_const_data_ptr != this->_last_data_ptr, "can't increment iterator past end");
  this->_const_data_ptr += _inc;
  return *this;
}

Vector_Values_Const_Iterator Vector_Values_Const_Iterator::operator++(int)
{
  Vector_Values_Const_Iterator tmp = *this;
  ++(*this);
  return tmp;
}

Vector_Values_Const_Iterator::reference Vector_Values_Const_Iterator::operator*(void) const
{
  REQUIRE(this->_const_data_ptr != this->_last_data_ptr, "can't dereference end iterator");
  return *this->_const_data_ptr;
}

Vector_Values_Const_Iterator::pointer Vector_Values_Const_Iterator::operator->(void) const
{
  REQUIRE(this->_const_data_ptr != this->_last_data_ptr, "can't reference end iterator");
  return this->_const_data_ptr;
}

bool Vector_Values_Const_Iterator::operator==(const Vector_Values_Const_Iterator other) const
{
  return this->_const_data_ptr == other._const_data_ptr && this->_inc == other._inc && this->_last_data_ptr == other._last_data_ptr;
}

bool Vector_Values_Const_Iterator::operator!=(const Vector_Values_Const_Iterator other) const
{
  return this->_const_data_ptr != other._const_data_ptr;
}

/*










*/

double Vector_Values_Const_Wrapper::operator[](const int position) const
{
  REQUIRE(0 <= position, "position should be positive number");
  REQUIRE(position < this->_num_values, "position should be less then size");
  return this->_const_values[position * this->_inc];
}

bool Vector_Values_Const_Wrapper::operator==(const Vector_Values_Const_Wrapper& other) const
{
  if (this->_num_values != other._num_values)
    return false;

  for (int i = 0; i < this->_num_values; ++i)
  {
    if (this->_const_values[i * this->_inc] != other._const_values[i * other._inc])
      return false;
  }

  return true;
}

bool Vector_Values_Const_Wrapper::operator!=(const Vector_Values_Const_Wrapper& other) const
{
  return !((*this) == other);
}

double Vector_Values_Const_Wrapper::at(const int position) const
{
  REQUIRE(position < this->_num_values, "position can't exceed given range");
  return this->_const_values[position * this->_inc];
}

Vector_Values_Const_Iterator Vector_Values_Const_Wrapper::begin(void) const
{
  Vector_Values_Const_Iterator begin(this->data(), this->_inc, this->end_ptr());
  return begin;
}

Vector_Values_Const_Iterator Vector_Values_Const_Wrapper::end(void) const
{
  const auto end_ptr = this->end_ptr();

  Vector_Values_Const_Iterator end(end_ptr, this->_inc, end_ptr);
  return end;
}

const double* Vector_Values_Const_Wrapper::data(void) const
{
  return this->_const_values.data();
}

bool Vector_Values_Const_Wrapper::empty(void) const
{
  return this->_const_values.empty();
}

int Vector_Values_Const_Wrapper::num_values(void) const
{
  return this->_num_values;
}

Vector_Values_Const_Wrapper Vector_Values_Const_Wrapper::part(const int start_pos, const int end_pos) const
{
  REQUIRE(start_pos < end_pos, "start position should be smaller than end position");
  REQUIRE(end_pos <= this->_num_values, "end position should be smaller or equal to dimension");

  const auto start_index = start_pos * this->_inc;
  const auto count       = (end_pos - start_pos - 1) * this->_inc + 1;

  const auto sub_const_values = this->_const_values.subspan(start_index, count);
  return {sub_const_values, this->_inc};
}

int Vector_Values_Const_Wrapper::inc(void) const
{
  return this->_inc;
}

std::string Vector_Values_Const_Wrapper::to_string(void) const
{
  std::ostringstream oss;
  oss << std::setprecision(16) << std::showpoint << std::left;
  for (int i = 0; i < this->_num_values; ++i)
  {
    oss << std::setw(25) << this->at(i);
  }
  return oss.str();
}

std::vector<double> Vector_Values_Const_Wrapper::to_vector(void) const
{
  std::vector<double> result(this->_num_values);
  ms::math::blas::copy(result.data(), this->data(), this->_num_values, 1, this->_inc);
  return result;
}

void Vector_Values_Const_Wrapper::initialize(void)
{
  REQUIRE(0 < this->_inc, "inc should greater than 0");

  this->_num_values = static_cast<int>(std::ceil(this->_const_values.size() / static_cast<double>(this->_inc)));
  REQUIRE(0 < this->_num_values, "number of values should greater than 0");
}

inline const double* Vector_Values_Const_Wrapper::end_ptr(void) const
{
  const auto last_index = static_cast<int>(std::ceil(this->_const_values.size() / static_cast<double>(this->_inc)) * this->_inc);
  return this->data() + last_index;
}

/*










*/

double& Vector_Values_Wrapper::operator[](const int index)
{
  REQUIRE(0 <= index, "index should be positive");
  REQUIRE(index < this->_num_values, "index should be less then size");
  return this->_values[index * this->_inc];
}

double& Vector_Values_Wrapper::at(const int index)
{
  REQUIRE(0 <= index, "index should be positive");
  REQUIRE(index < this->_num_values, "index should be less then size");
  return this->_values[index * this->_inc];
}

void Vector_Values_Wrapper::change_value(const Vector_Values_Const_Wrapper other)
{
  REQUIRE(this->_num_values == other.num_values(), "should be same size");
  ms::math::blas::copy(this->data(), other.data(), this->_num_values, this->_inc, other.inc());
}

double* Vector_Values_Wrapper::data(void)
{
  return this->_values.data();
}

} // namespace ms::math

namespace ms::math
{

Vector<0> Vector_Const_Wrapper::operator*(const double scalar) const
{
  constexpr auto incr = 1;

  const auto dim  = this->dimension();
  const auto data = this->data();
  const auto inc  = this->inc();

  Vector result(dim);

  ms::math::blas::cx(result.data(), scalar, data, dim, incr, inc);

  return result;
}

Vector<0> Vector_Const_Wrapper::operator+(const Vector_Const_Wrapper& other) const
{
  constexpr auto incr = 1;

  const auto dim        = this->dimension();
  const auto data       = this->data();
  const auto inc        = this->inc();
  const auto other_dim  = other.dimension();
  const auto other_data = other._values_cwrap.data();
  const auto other_inc  = other._values_cwrap.inc();
  REQUIRE(dim == other_dim, "other vector should be same dimension");

  Vector result(dim);
  ms::math::blas::x_plus_y(result.data(), data, other_data, dim, incr, inc, other_inc);

  return result;
}

Vector<0> Vector_Const_Wrapper::operator-(const Vector_Const_Wrapper& other) const
{
  constexpr auto incr = 1;

  const auto dim        = this->dimension();
  const auto data       = this->data();
  const auto inc        = this->inc();
  const auto other_dim  = other.dimension();
  const auto other_data = other._values_cwrap.data();
  const auto other_inc  = other._values_cwrap.inc();
  REQUIRE(dim == other_dim, "other vector should be same dimension");

  Vector result(dim);
  ms::math::blas::x_minus_y(result.data(), data, other_data, dim, incr, inc, other_inc);

  return result;
}

double Vector_Const_Wrapper::operator[](const int position) const
{
  return this->_values_cwrap[position];
}

bool Vector_Const_Wrapper::operator==(const Vector_Const_Wrapper& other) const
{
  return this->_values_cwrap == other._values_cwrap;
}

bool Vector_Const_Wrapper::operator!=(const Vector_Const_Wrapper& other) const
{
  return !((*this) == other);
}

// Vector_Const_Wrapper::operator std::pair<const double*, int>(void) const
//{
//   return {this->_coordinate_const_ptr, this->_inc};
// }

// double Vector_Const_Wrapper::at(const int position) const
//{
//   REQUIRE(position < this->_dimension, "position can't exceed given range");
//   return this->_coordinate_const_ptr[position * this->_inc];
// }

// const double* Vector_Const_Wrapper::begin(void) const
//{
//   return this->_coordinate_const_ptr;
// }

double Vector_Const_Wrapper::cosine(const Vector_Const_Wrapper& other) const
{
  const auto dim       = this->dimension();
  const auto other_dim = other.dimension();
  REQUIRE(dim == other_dim, "other vector should be same dimension");

  return this->inner_product(other) / (this->L2_norm() * other.L2_norm());
}

Vector_Const_Wrapper Vector_Const_Wrapper::const_wrapper(void) const
{
  return *this;
}

// const double* Vector_Const_Wrapper::data(void) const
//{
//   return this->_coordinate_const_ptr;
// }

// const double* Vector_Const_Wrapper::end(void) const
//{
//   const auto end_index = this->_inc * (this->_dimension - 1) + 1;
//   return this->_coordinate_const_ptr + end_index;
// }

bool Vector_Const_Wrapper::is_parallel(const Vector_Const_Wrapper& other) const
{
  const auto abs_cosine = std::abs(this->cosine(other));
  return compare_double(abs_cosine, 1.0);
}

double Vector_Const_Wrapper::L1_norm(void) const
{
  const auto dim  = this->dimension();
  const auto data = this->data();
  const auto inc  = this->inc();
  return ms::math::blas::abs_sum_x(data, dim, inc);
}

inline double Vector_Const_Wrapper::L2_norm(void) const
{
  return std::sqrt(this->inner_product(*this));
}

double Vector_Const_Wrapper::Linf_norm(void) const
{
  const auto dim  = this->dimension();
  const auto data = this->data();
  const auto inc  = this->inc();

  const auto pos = static_cast<int>(ms::math::blas::find_maximum_element_pos(data, dim, inc));
  return this->_values_cwrap[pos];
}

double Vector_Const_Wrapper::inner_product(const Vector_Const_Wrapper& other) const
{
  const auto dim        = this->dimension();
  const auto data       = this->data();
  const auto inc        = this->inc();
  const auto other_dim  = other.dimension();
  const auto other_data = other._values_cwrap.data();
  const auto other_inc  = other._values_cwrap.inc();
  REQUIRE(dim == other_dim, "other vector should be same dimension");

  return ms::math::blas::x_dot_y(data, other_data, dim, inc, other_inc);
}

inline int Vector_Const_Wrapper::dimension(void) const
{
  return this->_values_cwrap.num_values();
}

Vector_Const_Wrapper Vector_Const_Wrapper::part(const int start_pos, const int end_pos) const
{
  REQUIRE(start_pos < end_pos, "start position should be smaller than end position");
  REQUIRE(end_pos <= this->dimension(), "end position should be smaller or equal to dimension");

  return this->_values_cwrap.part(start_pos, end_pos);
}

std::string Vector_Const_Wrapper::to_string(void) const
{
  return this->_values_cwrap.to_string();
}

template <>
Vector<3> Vector_Const_Wrapper::cross_product<2>(const Vector_Const_Wrapper& other) const
{
  Vector<3> result;
  result[2] = (*this)[0] * other[1] - (*this)[1] * other[0];

  return result;
}

template <>
Vector<3> Vector_Const_Wrapper::cross_product<3>(const Vector_Const_Wrapper& other) const
{
  Vector<3> result;
  result[0] = (*this)[1] * other[2] - (*this)[2] * other[1];
  result[1] = (*this)[2] * other[0] - (*this)[0] * other[2];
  result[2] = (*this)[0] * other[1] - (*this)[1] * other[0];

  return result;
}

inline const double* Vector_Const_Wrapper::data(void) const
{
  return this->_values_cwrap.data();
}

inline int Vector_Const_Wrapper::inc(void) const
{
  return this->_values_cwrap.inc();
}

/*










*/

// Vector_Wrapper::Vector_Wrapper(Values_Wrapper values)
//{
//   this->_values_cwrap = values;
//   this->_values_wrap  = values;
// }

void Vector_Wrapper::operator*=(const double constant)
{
  const auto dim  = this->dimension();
  auto       data = this->data();
  const auto inc  = this->inc();

  ms::math::blas::cx(constant, data, dim, inc);
}

void Vector_Wrapper::operator+=(const Vector_Const_Wrapper& other)
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

void Vector_Wrapper::operator-=(const Vector_Const_Wrapper& other)
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

double& Vector_Wrapper::operator[](const int position)
{
  return this->_values_wrap[position];
}

// Vector_Wrapper::operator double*(void)
//{
//   return data;
// }

// double& Vector_Wrapper::at(const int position)
//{
//   REQUIRE(position < dim, "position should be less then size");
//   return data[position];
// }

void Vector_Wrapper::change_value(const Vector_Const_Wrapper other)
{
  const auto dim        = this->dimension();
  auto       data       = this->data();
  const auto inc        = this->inc();
  const auto other_dim  = other.dimension();
  const auto other_data = other.data();
  const auto other_inc  = other.inc();

  REQUIRE(dim == other.dimension(), "other vector should be same size");
  ms::math::blas::copy(data, other.data(), dim, inc, other.inc());
}

// void Vector_Wrapper::initalize(void)
//{
//   for (int i = 0; i < dim; ++i)
//   {
//     data[i] = 0.0;
//   }
// }

void Vector_Wrapper::normalize(void)
{
  *this *= 1.0 / this->L2_norm();
}

Vector_Wrapper Vector_Wrapper::wrapper(void)
{
  return *this;
}

double* Vector_Wrapper::data(void)
{
  return this->_values_wrap.data();
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

Vector<0>::Vector(Vector_Values_Const_Wrapper values_cwrap)
    : _coordinates(values_cwrap.to_vector())
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
  other._values_wrap  = Vector_Values_Wrapper();
  other._values_cwrap = Vector_Values_Const_Wrapper();
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
  this->_values_cwrap = Vector_Values_Const_Wrapper(this->_coordinates);
  this->_values_wrap  = Vector_Values_Wrapper(this->_coordinates);
}

/*










*/

Vector<0> operator*(const double constant, const Vector_Const_Wrapper& x)
{
  return x * constant;
}

} // namespace ms::math
