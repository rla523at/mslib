#pragma once
#include <array>
#include <string>
#include <vector>

#define COMPILE_TIME_REQUIREMENT static_assert

namespace ms::math
{
template <typename... Args>
inline constexpr bool are_arithmetics = (... && std::is_arithmetic_v<Args>);

template <typename T>
using if_is_iterator = std::enable_if_t<std::_Is_iterator_v<T>, bool>;

template <typename... Args>
using if_is_more_than_two = std::enable_if_t<2 <= sizeof...(Args), bool>;
} // namespace ms::math

/*





*/

// forward declaration
namespace ms::math
{
template <int dim>
class Vector;
}

/*





*/

// class declaration
namespace ms::math
{

class Vector_Const_Wrapper
{
public:
  Vector_Const_Wrapper(void) = default;
  Vector_Const_Wrapper(const double* coordinate_ptr, const int dimension, const int incx = 1);
  Vector_Const_Wrapper(const std::vector<double>& coordinates, const int incx = 1);

public:
  Vector<0> operator*(const double scalar) const;
  Vector<0> operator-(const Vector_Const_Wrapper& other) const;
  Vector<0> operator+(const Vector_Const_Wrapper& other) const;
  double    operator[](const size_t position) const;
  bool      operator==(const Vector_Const_Wrapper& other) const;
  bool      operator!=(const Vector_Const_Wrapper& other) const;

public:
  operator const double*(void) const;

public:
  double               at(const size_t position) const;
  const double*        begin(void) const;
  double               cosine(const Vector_Const_Wrapper& other) const;
  const double*        data(void) const;
  size_t               dimension(void) const;
  const double*        end(void) const;
  double               inner_product(const Vector_Const_Wrapper& other) const;
  int                  inc(void) const { return this->_inc; };
  bool                 is_parallel(const Vector_Const_Wrapper& other) const;
  double               L1_norm(void) const;
  double               L2_norm(void) const;
  double               Linf_norm(void) const;
  Vector_Const_Wrapper part(const int start_pos, const int end_pos) const; // [start_pos, end_pos)
  std::string          to_string(void) const;

public:
  template <int dim>
  Vector<3> cross_product(const Vector_Const_Wrapper& other) const;

protected:
  int           _dimension            = 0;
  const double* _coordinate_const_ptr = nullptr;
  int           _inc                  = 0;
};

class Vector_Wrapper : public Vector_Const_Wrapper
{
public:
  Vector_Wrapper(void) = default;
  Vector_Wrapper(double* coordinate_ptr, const int dimension, const int incx = 1)
      : Vector_Const_Wrapper(coordinate_ptr, dimension, incx), _coordinate_ptr(coordinate_ptr){};
  Vector_Wrapper(std::vector<double>& coordinates, const int incx = 1)
      : Vector_Const_Wrapper(coordinates, incx), _coordinate_ptr(coordinates.data()){};

public:
  void    operator*=(const double constant);
  void    operator+=(const Vector_Const_Wrapper& other);
  void    operator-=(const Vector_Const_Wrapper& other);
  double& operator[](const size_t position);

public:
  operator double*(void);

public:
  double& at(const size_t position);
  void    change_value(const Vector_Const_Wrapper& other);
  double* data(void);
  void    initalize(void);
  void    normalize(void);

  // resolve base class method hiding problem (overloading across scope problem)
public:
  using Vector_Const_Wrapper::operator[];
  using Vector_Const_Wrapper::at;
  using Vector_Const_Wrapper::data;

protected:
  double* _coordinate_ptr = nullptr;
};

template <int Dim = 0>
class Vector : public Vector_Wrapper
{
  COMPILE_TIME_REQUIREMENT(0 <= Dim, "dimension shold be positive");

public:
  Vector(void)
      : Vector_Wrapper(this->coordinates_.data(), Dim){};
  Vector(const Vector& other)
      : Vector_Wrapper(this->coordinates_.data(), Dim), coordinates_{other.coordinates_} {};

  template <typename... Args>
  Vector(const Args&... args)
      : Vector_Wrapper(this->coordinates_.data(), Dim), coordinates_{static_cast<double>(args)...}
  {
    COMPILE_TIME_REQUIREMENT(sizeof...(Args) <= Dim, "Number of arguments can not exceed dimension");
    COMPILE_TIME_REQUIREMENT(ms::math::are_arithmetics<Args...>, "Every arguments should be arithmetics");
  }

public:
  void operator=(const Vector& other)
  {
    this->coordinates_ = other.coordinates_;
  }

private:
  std::array<double, Dim> coordinates_ = {0};
};

// user-defined deduction guides
template <typename Iter, ms::math::if_is_iterator<Iter> = true>
Vector(Iter iter1, Iter iter2) -> Vector<0>;

template <typename... Args, ms::math::if_is_more_than_two<Args...> = true>
Vector(Args... args) -> Vector<sizeof...(Args)>;

template <>
class Vector<0> : public Vector_Wrapper
{
public:
  Vector(void) = default;
  explicit Vector(const int dimension);
  Vector(const std::initializer_list<double> list);
  Vector(const std::vector<double>& values);
  Vector(std::vector<double>&& values);
  Vector(const Vector& other);
  Vector(Vector&& other) noexcept;

  template <typename Iter, ms::math::if_is_iterator<Iter> = true>
  Vector(Iter first, Iter last)
      : coordinates_(first, last)
  {
    this->_dimension            = static_cast<int>(this->coordinates_.size());
    this->_coordinate_const_ptr = this->coordinates_.data();
    this->_coordinate_ptr       = this->coordinates_.data();
    this->_inc                  = 1;
  };

public:
  void operator=(const Vector& other);
  void operator=(Vector&& other) noexcept;

private:
  std::vector<double> coordinates_;
};

Vector<0> operator*(const double constant, const Vector_Const_Wrapper& x);
} // namespace ms::math
