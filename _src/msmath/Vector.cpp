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

namespace ms::math
{

Vector_Const_Wrapper::Vector_Const_Wrapper(const double* coordinate_const_ptr, const int dimension, const int incx)
{
  REQUIRE(0 < dimension, "dimension should greater than 0");
   REQUIRE(coordinate_const_ptr != nullptr, "data ptr should not be nullptr");
  this->_dimension            = dimension;
  this->_coordinate_const_ptr = coordinate_const_ptr;
  this->_inc                  = incx;
}

Vector_Const_Wrapper::Vector_Const_Wrapper(const std::vector<double>& coordinates, const int incx)
{
  const auto dimension = static_cast<int>(std::ceil(coordinates.size() / static_cast<double>(incx)));

  REQUIRE(0 < dimension, "dimension should greater than 0");
  this->_dimension            = dimension;
  this->_coordinate_const_ptr = coordinates.data();
  this->_inc                  = incx;
}

Vector<0> Vector_Const_Wrapper::operator*(const double scalar) const
{
  Vector         result(this->_dimension);
  constexpr auto incr = 1;

  ms::math::blas::cx(result.data(), scalar, this->_coordinate_const_ptr, this->_dimension, incr, this->_inc);

  return result;
}

Vector<0> Vector_Const_Wrapper::operator+(const Vector_Const_Wrapper& other) const
{
  REQUIRE(this->_dimension == other._dimension, "other vector should be same size");

  Vector         result(this->_dimension);
  constexpr auto incr = 1;

  ms::math::blas::x_plus_y(result.data(), this->_coordinate_const_ptr, other._coordinate_const_ptr, this->_dimension, incr, this->_inc, other._inc);

  return result;
}

Vector<0> Vector_Const_Wrapper::operator-(const Vector_Const_Wrapper& other) const
{
  REQUIRE(this->_dimension == other._dimension, "other vector should be same size");

  Vector         result(this->_dimension);
  constexpr auto incr = 1;

  ms::math::blas::x_minus_y(result.data(), this->_coordinate_const_ptr, other._coordinate_const_ptr, this->_dimension, incr, this->_inc, other._inc);

  return result;
}

double Vector_Const_Wrapper::operator[](const size_t position) const
{
  REQUIRE(position < this->_dimension, "position should be less then size");
  return this->_coordinate_const_ptr[position * this->_inc];
}

bool Vector_Const_Wrapper::operator==(const Vector_Const_Wrapper& other) const
{
  if (this->_dimension != other._dimension)
    return false;

  for (int i = 0; i < this->_dimension; ++i)
  {
    if (this->_coordinate_const_ptr[i * this->_inc] != other._coordinate_const_ptr[i * other._inc])
      return false;
  }

  return true;
}

bool Vector_Const_Wrapper::operator!=(const Vector_Const_Wrapper& other) const
{
  return !((*this) == other);
}

Vector_Const_Wrapper::operator const double*(void) const
{
  return this->_coordinate_const_ptr;
}

double Vector_Const_Wrapper::at(const size_t position) const
{
  REQUIRE(position < this->_dimension, "position can't exceed given range");
  return this->_coordinate_const_ptr[position * this->_inc];
}

const double* Vector_Const_Wrapper::begin(void) const
{
  return this->_coordinate_const_ptr;
}

double Vector_Const_Wrapper::cosine(const Vector_Const_Wrapper& other) const
{
  REQUIRE(this->_dimension == other._dimension, "two vector should be same size");

  return this->inner_product(other) / (this->L2_norm() * other.L2_norm());
}

const double* Vector_Const_Wrapper::data(void) const
{
  return this->_coordinate_const_ptr;
}

const double* Vector_Const_Wrapper::end(void) const
{
  const auto end_index = this->_inc * (this->_dimension - 1) + 1;
  return this->_coordinate_const_ptr + end_index;
}

bool Vector_Const_Wrapper::is_parallel(const Vector_Const_Wrapper& other) const
{
  const auto abs_cosine = std::abs(this->cosine(other));
  return compare_double(abs_cosine, 1.0);
}

double Vector_Const_Wrapper::L1_norm(void) const
{
  return ms::math::blas::abs_sum_x(this->_coordinate_const_ptr, this->_dimension, this->_inc);
}

double Vector_Const_Wrapper::L2_norm(void) const
{
  return std::sqrt(this->inner_product(*this));
}

double Vector_Const_Wrapper::Linf_norm(void) const
{
  const auto pos = ms::math::blas::find_maximum_element_pos(this->_coordinate_const_ptr, this->_dimension, this->_inc);
  return this->at(pos);
}

double Vector_Const_Wrapper::inner_product(const Vector_Const_Wrapper& other) const
{
  REQUIRE(this->_dimension == other._dimension, "two vector should be same size");

  return ms::math::blas::x_dot_y(this->_coordinate_const_ptr, other._coordinate_const_ptr, this->_dimension);
}

size_t Vector_Const_Wrapper::dimension(void) const
{
  return this->_dimension;
}

Vector_Const_Wrapper Vector_Const_Wrapper::part(const int start_pos, const int end_pos) const
{
  REQUIRE(start_pos < end_pos, "start position should be smaller than end position");
  REQUIRE(end_pos <= this->_dimension, "end position should be smaller or equal to dimension");

  const auto start_index       = start_pos * this->_inc;
  const auto start_ptr         = this->_coordinate_const_ptr + start_index;
  const auto dimension_of_part = end_pos - start_pos;

  return {start_ptr, dimension_of_part, this->_inc};
}

std::string Vector_Const_Wrapper::to_string(void) const
{
  std::ostringstream oss;
  oss << std::setprecision(16) << std::showpoint << std::left;
  for (int i = 0; i < this->_dimension; ++i)
  {
    oss << std::setw(25) << this->_coordinate_const_ptr[i];
  }
  return oss.str();
}

template <>
Vector<3> Vector_Const_Wrapper::cross_product<2>(const Vector_Const_Wrapper& other) const
{
  Vector<3> result;
  result[2] = this->at(0) * other.at(1) - this->at(1) * other.at(0);

  return result;
}

template <>
Vector<3> Vector_Const_Wrapper::cross_product<3>(const Vector_Const_Wrapper& other) const
{
  Vector<3> result;
  result[0] = this->at(1) * other.at(2) - this->at(2) * other.at(1);
  result[1] = this->at(2) * other.at(0) - this->at(0) * other.at(2);
  result[2] = this->at(0) * other.at(1) - this->at(1) * other.at(0);

  return result;
}

void Vector_Wrapper::operator*=(const double constant)
{
  ms::math::blas::cx(constant, this->_coordinate_ptr, this->_dimension, this->_inc);
}

void Vector_Wrapper::operator+=(const Vector_Const_Wrapper& other)
{
  REQUIRE(this->_dimension == other.dimension(), "other vector should be same size");
  ms::math::blas::x_plus_assign_y(this->_coordinate_ptr, other.data(), this->_dimension, this->_inc, other.inc());
}

void Vector_Wrapper::operator-=(const Vector_Const_Wrapper& other)
{
  REQUIRE(this->_dimension == other.dimension(), "other vector should be same size");
  ms::math::blas::x_minus_assign_y(this->_coordinate_ptr, other.data(), this->_dimension, this->_inc, other.inc());
}

double& Vector_Wrapper::operator[](const size_t position)
{
  REQUIRE(position < this->_dimension, "position should be less then size");
  return this->_coordinate_ptr[position];
}

Vector_Wrapper::operator double*(void)
{
  return this->_coordinate_ptr;
}

double& Vector_Wrapper::at(const size_t position)
{
  REQUIRE(position < this->_dimension, "position should be less then size");
  return this->_coordinate_ptr[position];
}

void Vector_Wrapper::change_value(const Vector_Const_Wrapper& other)
{
  REQUIRE(this->_dimension == other.dimension(), "other vector should be same size");
  ms::math::blas::copy(this->_coordinate_ptr, other.data(), this->_dimension, this->_inc, other.inc());
}

double* Vector_Wrapper::data(void)
{
  return this->_coordinate_ptr;
}

void Vector_Wrapper::initalize(void)
{
  for (int i = 0; i < this->_dimension; ++i)
  {
    this->_coordinate_ptr[i] = 0.0;
  }
}

void Vector_Wrapper::normalize(void)
{
  *this *= 1.0 / this->L2_norm();
}

Vector<0>::Vector(const int dimension)
    : _coordinates(dimension)
{
  this->_dimension            = static_cast<int>(this->_coordinates.size());
  this->_coordinate_const_ptr = this->_coordinates.data();
  this->_coordinate_ptr       = this->_coordinates.data();
  this->_inc                  = 1;
}

Vector<0>::Vector(const std::initializer_list<double> list)
    : _coordinates(list)
{
  this->_dimension            = static_cast<int>(this->_coordinates.size());
  this->_coordinate_const_ptr = this->_coordinates.data();
  this->_coordinate_ptr       = this->_coordinates.data();
  this->_inc                  = 1;
}

Vector<0>::Vector(const std::vector<double>& values)
    : _coordinates(values)
{
  this->_dimension            = static_cast<int>(this->_coordinates.size());
  this->_coordinate_const_ptr = this->_coordinates.data();
  this->_coordinate_ptr       = this->_coordinates.data();
  this->_inc                  = 1;
}

Vector<0>::Vector(std::vector<double>&& values)
    : _coordinates(std::move(values))
{
  this->_dimension            = static_cast<int>(this->_coordinates.size());
  this->_coordinate_const_ptr = this->_coordinates.data();
  this->_coordinate_ptr       = this->_coordinates.data();
  this->_inc                  = 1;
};

Vector<0>::Vector(const Vector& other)
    : _coordinates(other._coordinates)
{
  this->_dimension            = other._dimension;
  this->_coordinate_const_ptr = this->_coordinates.data();
  this->_coordinate_ptr       = this->_coordinates.data();
  this->_inc                  = other._inc;
}

Vector<0>::Vector(Vector&& other) noexcept
    : _coordinates(std::move(other._coordinates))
{
  this->_dimension            = other._dimension;
  this->_coordinate_const_ptr = this->_coordinates.data();
  this->_coordinate_ptr       = this->_coordinates.data();
  this->_inc                  = other._inc;

  other._dimension            = 0;
  other._coordinate_const_ptr = nullptr;
  other._coordinate_ptr       = nullptr;
}

void Vector<0>::operator=(const Vector& other)
{
  this->_coordinates          = other._coordinates;
  this->_dimension            = other._dimension;
  this->_coordinate_const_ptr = this->_coordinates.data();
  this->_coordinate_ptr       = this->_coordinates.data();
  this->_inc                  = other._inc;
}

void Vector<0>::operator=(Vector&& other) noexcept
{
  this->_coordinates          = std::move(other._coordinates);
  this->_dimension            = other._dimension;
  this->_coordinate_const_ptr = this->_coordinates.data();
  this->_coordinate_ptr       = this->_coordinates.data();
  this->_inc                  = other._inc;
}

Vector<0> operator*(const double constant, const Vector_Const_Wrapper& x)
{
  return x * constant;
}

} // namespace ms::math
