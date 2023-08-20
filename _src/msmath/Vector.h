#pragma once
#include "msexception/Exception.h"
#include <array>
#include <span>
#include <string>
#include <vector>

// forward declaration
namespace ms::math
{
class Vector_Wrapper;

template <int dim>
class Vector;
} // namespace ms::math

/*










*/

// miscellaneous definitions
namespace ms::math
{
#define COMPILE_TIME_REQUIREMENT static_assert

template <typename T>
concept const_span = requires(const T& t) { std::span<const double>(t); };

template <typename T>
concept span = requires(T& t) { std::span<double>(t); };

template <typename T>
concept contiguous_iterator_or_int = std::contiguous_iterator<T> || std::same_as<T, int>;

template <typename T>
concept output_iterator = std::output_iterator<T, typename T::value_type>;

template <typename T>
concept output_iterator_or_int = output_iterator<T> || std::same_as<T, int>;

template <typename... Args>
inline constexpr bool are_arithmetics = (... && std::is_arithmetic_v<Args>);

template <typename... Args>
concept more_than_two = (2 <= sizeof...(Args));

} // namespace ms::math

/*










*/

namespace ms::math
{
class Vector_Values_Const_Iterator
{
public:
  using iterator_category = std::input_iterator_tag;
  using difference_type   = std::ptrdiff_t;
  using value_type        = double;
  using pointer           = const double*;
  using reference         = const double&;

public:
  Vector_Values_Const_Iterator(const double* data_ptr, const int inc, const double* last_data_ptr)
      : _const_data_ptr(data_ptr),
        _inc(inc),
        _last_data_ptr(last_data_ptr){};

public:
  Vector_Values_Const_Iterator& operator++();
  Vector_Values_Const_Iterator  operator++(int); // Postfix increment

public:
  reference operator*(void) const;
  pointer   operator->(void) const;
  bool      operator==(const Vector_Values_Const_Iterator other) const;
  bool      operator!=(const Vector_Values_Const_Iterator other) const;

private:
  pointer _const_data_ptr = nullptr;
  int     _inc            = 1;
  pointer _last_data_ptr  = nullptr;
};
} // namespace ms::math

/*










*/

// Vector_Values_Const_Wrapper class declaration
namespace ms::math
{
class Vector_Values_Const_Wrapper
{
public:
  using const_iterator = Vector_Values_Const_Iterator;

public:
  Vector_Values_Const_Wrapper(void) = default;

  template <const_span T>
  Vector_Values_Const_Wrapper(const T& values, const int inc = 1);

  template <std::contiguous_iterator T, contiguous_iterator_or_int U>
  Vector_Values_Const_Wrapper(const T& t, const U& u, const int inc = 1);

public:
  double operator[](const int position) const;
  bool   operator==(const Vector_Values_Const_Wrapper& other) const;
  bool   operator!=(const Vector_Values_Const_Wrapper& other) const;

public:
  double                       at(const int position) const;
  Vector_Values_Const_Iterator begin(void) const;
  Vector_Values_Const_Iterator end(void) const;
  const double*                data(void) const;
  bool                         empty(void) const;
  int                          num_values(void) const;
  int                          inc(void) const;
  Vector_Values_Const_Wrapper  part(const int start_position, const int end_position) const;
  std::string                  to_string(void) const;
  std::vector<double>          to_vector(void) const;

private:
  void initialize(void);

private:
  const double* end_ptr(void) const;

protected:
  std::span<const double> _const_values;
  int                     _num_values = 0;
  int                     _inc        = 1;
};

} // namespace ms::math

/*










*/

// Vector_Values_Wrapper class declaration
namespace ms::math
{

class Vector_Values_Wrapper : public Vector_Values_Const_Wrapper
{
public:
  Vector_Values_Wrapper(void) = default;

  template <span T>
  Vector_Values_Wrapper(T& values, const int inc = 1)
      : Vector_Values_Const_Wrapper(values, inc),
        _values(values){};

  template <output_iterator T, output_iterator_or_int U>
  Vector_Values_Wrapper(T& t, const U& u, const int inc = 1)
      : Vector_Values_Const_Wrapper(t, u, inc),
        _values(t, u){};

public:
  double& operator[](const int position);

public:
  double& at(const int position);
  double* data(void);
  void    change_value(const Vector_Values_Const_Wrapper other);

public:
  using Vector_Values_Const_Wrapper::operator[];
  using Vector_Values_Const_Wrapper::at;
  using Vector_Values_Const_Wrapper::data;

private:
  std::span<double> _values;
};

} // namespace ms::math

/*










*/

// miscellaneous definitions
namespace ms::math
{

template <typename T>
concept vector_const_values = requires(const T& t) { ms::math::Vector_Values_Const_Wrapper(t); };

template <typename T>
concept values_wrapper_constructable = requires(T& t) { ms::math::Vector_Values_Wrapper(t); };

} // namespace ms::math

/*










*/

// Vector_Const_Wrapper class declaration
namespace ms::math
{
class Vector_Const_Wrapper
{
public:
  Vector_Const_Wrapper(void) = default;

  Vector_Const_Wrapper(const Vector_Values_Const_Wrapper values)
      : _values_cwrap(values){};

  template <const_span T>
  Vector_Const_Wrapper(const T& values, const int inc = 1)
      : _values_cwrap(values, inc){};

  template <std::input_iterator T, contiguous_iterator_or_int U>
  Vector_Const_Wrapper(const T& t, const U& u, const int inc = 1)
      : _values_cwrap(t, u, inc){};

public:
  Vector<0> operator*(const double scalar) const;
  Vector<0> operator-(const Vector_Const_Wrapper& other) const;
  Vector<0> operator+(const Vector_Const_Wrapper& other) const;
  double    operator[](const int position) const;
  bool      operator==(const Vector_Const_Wrapper& other) const;
  bool      operator!=(const Vector_Const_Wrapper& other) const;

public:
  double               cosine(const Vector_Const_Wrapper& other) const;
  Vector_Const_Wrapper const_wrapper(void) const;
  const double*        data(void) const;
  int                  dimension(void) const;
  double               inner_product(const Vector_Const_Wrapper& other) const;
  bool                 is_parallel(const Vector_Const_Wrapper& other) const;
  int                  inc(void) const;
  double               L1_norm(void) const;
  double               L2_norm(void) const;
  double               Linf_norm(void) const;
  Vector_Const_Wrapper part(const int start_pos, const int end_pos) const;
  std::string          to_string(void) const;

public:
  template <int dim>
  Vector<3> cross_product(const Vector_Const_Wrapper& other) const;

protected:
  Vector_Values_Const_Wrapper _values_cwrap;
};

} // namespace ms::math

/*










*/

// Vector_Wrapper class declaration
namespace ms::math
{

class Vector_Wrapper : public Vector_Const_Wrapper
{
public:
  Vector_Wrapper(void) = default;

  Vector_Wrapper(Vector_Values_Wrapper values)
      : Vector_Const_Wrapper(values),
        _values_wrap(values){};

  template <span T>
  Vector_Wrapper(T& values, const int inc = 1)
      : Vector_Const_Wrapper(values, inc),
        _values_wrap(values, inc){};

  template <output_iterator T, output_iterator_or_int U>
  Vector_Wrapper(const T& t, const U& u, const int inc = 1)
      : Vector_Const_Wrapper(t, u, inc),
        _values_wrap(t, u, inc){};

public:
  void    operator*=(const double constant);
  void    operator+=(const Vector_Const_Wrapper& other);
  void    operator-=(const Vector_Const_Wrapper& other);
  double& operator[](const int position);

public:
  void           change_value(const Vector_Const_Wrapper other);
  void           normalize(void);
  Vector_Wrapper wrapper(void);

  // resolve base class method hiding problem (overloading across scope problem)
public:
  using Vector_Const_Wrapper::operator[];

protected:
  double* data(void);

protected:
  Vector_Values_Wrapper _values_wrap;
};

} // namespace ms::math

/*










*/

// vector<dim> class declaration
namespace ms::math
{

template <int Dim = 0>
class Vector : public Vector_Wrapper
{
  COMPILE_TIME_REQUIREMENT(0 <= Dim, "dimension shold be positive");

public:
  Vector(void);
  Vector(const Vector& other);

  template <typename... Args>
  Vector(const Args&... args);

public:
  void operator=(const Vector& other);

private:
  void reallocate_values(void);

private:
  std::array<double, Dim> _coordinates = {0};
};

} // namespace ms::math

/*










*/

// Vector<0> class declaration
namespace ms::math
{

template <>
class Vector<0> : public Vector_Wrapper
{
public:
  Vector(void) = default;
  explicit Vector(const int dimension); // prevent implicit convsersion to int
  Vector(const std::initializer_list<double> list);
  Vector(Vector_Values_Const_Wrapper values_cwrap);
  Vector(const Vector& other);
  Vector(Vector&& other) noexcept;

public:
  template <std::input_iterator Iter>
  Vector(Iter first, Iter last);

public:
  void operator=(const Vector& other);
  void operator=(Vector&& other) noexcept;

public:
  double* data(void);
  void    resize(const int size);

private:
  void reallocate_values(void);

private:
  std::vector<double> _coordinates;
};

} // namespace ms::math

/*










*/

// user-defined deduction guides
namespace ms::math
{
template <std::input_iterator Iter>
Vector(Iter iter1, Iter iter2) -> Vector<0>;

template <typename... Args>
requires more_than_two<Args...>
Vector(Args... args) -> Vector<sizeof...(Args)>;

} // namespace ms::math

/*










*/

// template definitions
namespace ms::math
{
template <ms::math::const_span T>
Vector_Values_Const_Wrapper::Vector_Values_Const_Wrapper(const T& values, const int inc)
    : _const_values(values),
      _inc(inc) { this->initialize(); };

template <std::contiguous_iterator T, ms::math::contiguous_iterator_or_int U>
Vector_Values_Const_Wrapper::Vector_Values_Const_Wrapper(const T& t, const U& u, const int inc)
    : _const_values(t, u),
      _inc(inc)
{
  this->initialize();
};

template <int Dim>
Vector<Dim>::Vector(void) { this->reallocate_values(); };

template <int Dim>
Vector<Dim>::Vector(const Vector& other)
    : _coordinates{other._coordinates}
{
  this->reallocate_values();
};

template <int Dim>
template <typename... Args>
Vector<Dim>::Vector(const Args&... args)
    : _coordinates{static_cast<double>(args)...}
{
  COMPILE_TIME_REQUIREMENT(sizeof...(Args) <= Dim, "Number of arguments can not exceed dimension");
  COMPILE_TIME_REQUIREMENT(ms::math::are_arithmetics<Args...>, "Every arguments should be arithmetics");
  this->reallocate_values();
}

template <int Dim>
void Vector<Dim>::operator=(const Vector& other)
{
  this->_coordinates = other._coordinates;
}

template <int Dim>
void Vector<Dim>::reallocate_values(void)
{
  this->_values_cwrap = Vector_Values_Const_Wrapper(this->_coordinates);
  this->_values_wrap  = Vector_Values_Wrapper(this->_coordinates);
}

template <std::input_iterator Iter>
Vector<0>::Vector(Iter first, Iter last)
    : _coordinates(first, last)
{
  this->reallocate_values();
};

} // namespace ms::math

/*










*/

// free function decalarations
namespace ms::math
{
Vector<0> operator*(const double constant, const Vector_Const_Wrapper& x);
}
