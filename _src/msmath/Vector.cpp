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
  this->dimension_            = dimension;
  this->coordinate_const_ptr_ = coordinate_const_ptr;
  this->inc_                  = incx;
}

Vector_Const_Wrapper::Vector_Const_Wrapper(const std::vector<double>& coordinates, const int incx)
{
  const auto dimension = static_cast<int>(std::ceil(coordinates.size() / static_cast<double>(incx)));

  REQUIRE(0 < dimension, "dimension should greater than 0");
  this->dimension_            = dimension;
  this->coordinate_const_ptr_ = coordinates.data();
  this->inc_                  = incx;
}

Vector<0> Vector_Const_Wrapper::operator*(const double scalar) const
{
  Vector         result(this->dimension_);
  constexpr auto incr = 1;

  ms::math::blas::cx(result.data(), scalar, this->coordinate_const_ptr_, this->dimension_, incr, this->inc_);

  return result;
}

Vector<0> Vector_Const_Wrapper::operator+(const Vector_Const_Wrapper& other) const
{
  REQUIRE(this->dimension_ == other.dimension_, "other vector should be same size");

  Vector         result(this->dimension_);
  constexpr auto incr = 1;

  ms::math::blas::x_plus_y(result.data(), this->coordinate_const_ptr_, other.coordinate_const_ptr_, this->dimension_, incr, this->inc_, other.inc_);

  return result;
}

Vector<0> Vector_Const_Wrapper::operator-(const Vector_Const_Wrapper& other) const
{
  REQUIRE(this->dimension_ == other.dimension_, "other vector should be same size");

  Vector         result(this->dimension_);
  constexpr auto incr = 1;

  ms::math::blas::x_minus_y(result.data(), this->coordinate_const_ptr_, other.coordinate_const_ptr_, this->dimension_, incr, this->inc_, other.inc_);

  return result;
}

double Vector_Const_Wrapper::operator[](const size_t position) const
{
  REQUIRE(position < this->dimension_, "position should be less then size");
  return this->coordinate_const_ptr_[position * this->inc_];
}

bool Vector_Const_Wrapper::operator==(const Vector_Const_Wrapper& other) const
{
  if (this->dimension_ != other.dimension_)
    return false;

  for (int i = 0; i < this->dimension_; ++i)
  {
    if (this->coordinate_const_ptr_[i * this->inc_] != other.coordinate_const_ptr_[i * other.inc_])
      return false;
  }

  return true;
}

bool Vector_Const_Wrapper::operator!=(const Vector_Const_Wrapper& other) const
{
  return !((*this) == other);
}

double Vector_Const_Wrapper::at(const size_t position) const
{
  REQUIRE(position < this->dimension_, "position can't exceed given range");
  return this->coordinate_const_ptr_[position * this->inc_];
}

const double* Vector_Const_Wrapper::begin(void) const
{
  return this->coordinate_const_ptr_;
}

double Vector_Const_Wrapper::cosine(const Vector_Const_Wrapper& other) const
{
  REQUIRE(this->dimension_ == other.dimension_, "two vector should be same size");

  return this->inner_product(other) / (this->L2_norm() * other.L2_norm());
}

const double* Vector_Const_Wrapper::data(void) const
{
  return this->coordinate_const_ptr_;
}

const double* Vector_Const_Wrapper::end(void) const
{
  const auto end_index = this->inc_ * (this->dimension_ - 1) + 1;
  return this->coordinate_const_ptr_ + end_index;
}

bool Vector_Const_Wrapper::is_parallel(const Vector_Const_Wrapper& other) const
{
  const auto abs_cosine = std::abs(this->cosine(other));
  return compare_double(abs_cosine, 1.0);
}

double Vector_Const_Wrapper::L1_norm(void) const
{
  return ms::math::blas::abs_sum_x(this->coordinate_const_ptr_, this->dimension_, this->inc_);
}

double Vector_Const_Wrapper::L2_norm(void) const
{
  return std::sqrt(this->inner_product(*this));
}

double Vector_Const_Wrapper::Linf_norm(void) const
{
  const auto pos = ms::math::blas::find_maximum_element_pos(this->coordinate_const_ptr_, this->dimension_, this->inc_);
  return this->at(pos);
}

double Vector_Const_Wrapper::inner_product(const Vector_Const_Wrapper& other) const
{
  REQUIRE(this->dimension_ == other.dimension_, "two vector should be same size");

  return ms::math::blas::x_dot_y(this->coordinate_const_ptr_, other.coordinate_const_ptr_, this->dimension_);
}

size_t Vector_Const_Wrapper::dimension(void) const
{
  return this->dimension_;
}

Vector_Const_Wrapper Vector_Const_Wrapper::part(const int start_pos, const int end_pos) const
{
  REQUIRE(start_pos < end_pos, "start position should be smaller than end position");
  REQUIRE(end_pos <= this->dimension_, "end position should be smaller or equal to dimension");

  const auto start_index       = start_pos * this->inc_;
  const auto start_ptr         = this->coordinate_const_ptr_ + start_index;
  const auto dimension_of_part = end_pos - start_pos;

  return {start_ptr, dimension_of_part, this->inc_};
}

std::string Vector_Const_Wrapper::to_string(void) const
{
  std::ostringstream oss;
  oss << std::setprecision(16) << std::showpoint << std::left;
  for (int i = 0; i < this->dimension_; ++i)
  {
    oss << std::setw(25) << this->coordinate_const_ptr_[i];
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
  ms::math::blas::cx(constant, this->coordinate_ptr_, this->dimension_, this->inc_);
}

void Vector_Wrapper::operator+=(const Vector_Const_Wrapper& other)
{
  REQUIRE(this->dimension_ == other.dimension(), "other vector should be same size");
  ms::math::blas::x_plus_assign_y(this->coordinate_ptr_, other.data(), this->dimension_, this->inc_, other.inc());
}

void Vector_Wrapper::operator-=(const Vector_Const_Wrapper& other)
{
  REQUIRE(this->dimension_ == other.dimension(), "other vector should be same size");
  ms::math::blas::x_minus_assign_y(this->coordinate_ptr_, other.data(), this->dimension_, this->inc_, other.inc());
}

double& Vector_Wrapper::operator[](const size_t position)
{
  REQUIRE(position < this->dimension_, "position should be less then size");
  return this->coordinate_ptr_[position];
}

double& Vector_Wrapper::at(const size_t position)
{
  REQUIRE(position < this->dimension_, "position should be less then size");
  return this->coordinate_ptr_[position];
}

void Vector_Wrapper::change_value(const Vector_Const_Wrapper& other)
{
  REQUIRE(this->dimension_ == other.dimension(), "other vector should be same size");
  ms::math::blas::copy(this->coordinate_ptr_, other.data(), this->dimension_, this->inc_, other.inc());
}

double* Vector_Wrapper::data(void)
{
  return this->coordinate_ptr_;
}

void Vector_Wrapper::initalize(void)
{
  for (int i = 0; i < this->dimension_; ++i)
  {
    this->coordinate_ptr_[i] = 0.0;
  }
}

void Vector_Wrapper::normalize(void)
{
  *this *= 1.0 / this->L2_norm();
}

Vector<0>::Vector(const int dimension)
    : coordinates_(dimension)
{
  this->dimension_            = static_cast<int>(this->coordinates_.size());
  this->coordinate_const_ptr_ = this->coordinates_.data();
  this->coordinate_ptr_       = this->coordinates_.data();
  this->inc_                  = 1;
}

Vector<0>::Vector(const std::initializer_list<double> list)
    : coordinates_(list)
{
  this->dimension_            = static_cast<int>(this->coordinates_.size());
  this->coordinate_const_ptr_ = this->coordinates_.data();
  this->coordinate_ptr_       = this->coordinates_.data();
  this->inc_                  = 1;
}

Vector<0>::Vector(const std::vector<double>& values)
    : coordinates_(values)
{
  this->dimension_            = static_cast<int>(this->coordinates_.size());
  this->coordinate_const_ptr_ = this->coordinates_.data();
  this->coordinate_ptr_       = this->coordinates_.data();
  this->inc_                  = 1;
}

Vector<0>::Vector(std::vector<double>&& values)
    : coordinates_(std::move(values))
{
  this->dimension_            = static_cast<int>(this->coordinates_.size());
  this->coordinate_const_ptr_ = this->coordinates_.data();
  this->coordinate_ptr_       = this->coordinates_.data();
  this->inc_                  = 1;
};

Vector<0>::Vector(const Vector& other)
    : coordinates_(other.coordinates_)
{
  this->dimension_            = other.dimension_;
  this->coordinate_const_ptr_ = this->coordinates_.data();
  this->coordinate_ptr_       = this->coordinates_.data();
  this->inc_                  = other.inc_;
}

Vector<0>::Vector(Vector&& other) noexcept
    : coordinates_(std::move(other.coordinates_))
{
  this->dimension_            = other.dimension_;
  this->coordinate_const_ptr_ = this->coordinates_.data();
  this->coordinate_ptr_       = this->coordinates_.data();
  this->inc_                  = other.inc_;

  other.dimension_            = 0;
  other.coordinate_const_ptr_ = nullptr;
  other.coordinate_ptr_       = nullptr;
}

void Vector<0>::operator=(const Vector& other)
{
  this->coordinates_          = other.coordinates_;
  this->dimension_            = other.dimension_;
  this->coordinate_const_ptr_ = this->coordinates_.data();
  this->coordinate_ptr_       = this->coordinates_.data();
  this->inc_                  = other.inc_;
}

void Vector<0>::operator=(Vector&& other) noexcept
{
  this->coordinates_          = std::move(other.coordinates_);
  this->dimension_            = other.dimension_;
  this->coordinate_const_ptr_ = this->coordinates_.data();
  this->coordinate_ptr_       = this->coordinates_.data();
  this->inc_                  = other.inc_;
}

Vector<0> operator*(const double constant, const Vector_Const_Wrapper& x)
{
  return x * constant;
}

} // namespace ms::math
