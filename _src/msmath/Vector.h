#pragma once
#include "msexception/Exception.h"
#include <array>
#include <span>
#include <string>
#include <vector>

// forward declaration
namespace ms::math
{
class Vector_Wrap;

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

// Vector_Const_Iterator Class Declaration
namespace ms::math
{
class Vector_Const_Iterator
{
public:
  using iterator_category = std::input_iterator_tag;
  using difference_type   = std::ptrdiff_t;
  using value_type        = double;
  using pointer           = const double*;
  using reference         = const double&;

public:
  Vector_Const_Iterator(const double* data_ptr, const int inc, const double* last_data_ptr)
      : _const_data_ptr(data_ptr),
        _inc(inc),
        _last_data_ptr(last_data_ptr){};

public:
  Vector_Const_Iterator& operator++();
  Vector_Const_Iterator  operator++(int); // Postfix increment

public:
  reference operator*(void) const;
  pointer   operator->(void) const;
  bool      operator==(const Vector_Const_Iterator other) const;
  bool      operator!=(const Vector_Const_Iterator other) const;

private:
  pointer _const_data_ptr = nullptr;
  int     _inc            = 1;
  pointer _last_data_ptr  = nullptr;
};
} // namespace ms::math

/*










*/

// Vector_View Class Declaration
namespace ms::math
{

class Vector_View
{
public:
  Vector_View(void) = default;

  template <const_span T>
  Vector_View(const T& values, const int inc = 1);

  template <std::contiguous_iterator T, contiguous_iterator_or_int U>
  Vector_View(const T& t, const U& u, const int inc = 1);

public:
  Vector<0> operator*(const double scalar) const;
  Vector<0> operator-(const Vector_View other) const;
  Vector<0> operator+(const Vector_View other) const;
  double    operator[](const int position) const;
  bool      operator==(const Vector_View other) const;
  bool      operator!=(const Vector_View other) const;

public:
  double                at(const int position) const;
  Vector_Const_Iterator begin(void) const;
  Vector_Const_Iterator end(void) const;
  double                cosine(const Vector_View other) const;
  int                   dimension(void) const;
  const double*         data(void) const;
  bool                  empty(void) const;
  Vector_View           get_view(void) const;
  int                   inc(void) const;
  double                inner_product(const Vector_View other) const;
  bool                  is_parallel(const Vector_View other) const;
  double                L1_norm(void) const;
  double                L2_norm(void) const;
  double                Linf_norm(void) const;
  Vector_View           sub_view(const int start_position, const int end_position) const;
  Vector_View           sub_view(const int start_position, const int inc, const int num_values) const;
  std::string           to_string(void) const;
  std::vector<double>   to_vector(void) const;

public:
  template <int dim>
  Vector<3> cross_product(const Vector_View other) const;

private:
  void initialize(void);

private:
  const double* end_ptr(void) const;

protected:
  std::span<const double> _values_view;
  int                     _dimension = 0;
  int                     _inc       = 1;
};

} // namespace ms::math

/*










*/

// Vector_Wrap Class Declaration
namespace ms::math
{

class Vector_Wrap : public Vector_View
{
public:
  Vector_Wrap(void) = default;

  template <const_span T>
  Vector_Wrap(T& values, const int inc = 1)
      : Vector_View(values, inc),
        _values_wrap(values){};

  template <std::contiguous_iterator T, contiguous_iterator_or_int U>
  Vector_Wrap(const T& t, const U& u, const int inc = 1)
      : Vector_View(t, u, inc),
        _values_wrap(t, u){};

public:
  void    operator*=(const double constant);
  void    operator+=(const Vector_View other);
  void    operator-=(const Vector_View other);
  double& operator[](const int position);

public:
  double&     at(const int position);
  void        change_value(const Vector_View other);
  double*     data(void);
  void        normalize(void);
  Vector_Wrap sub_wrap(const int start_position, const int end_position) const;
  Vector_Wrap sub_wrap(const int start_position, const int inc, const int num_values) const;

  // resolve base class method hiding problem (overloading across scope problem)
public:
  using Vector_View::at;
  using Vector_View::data;
  using Vector_View::operator[];

protected:
  std::span<double> _values_wrap;
};

} // namespace ms::math

/*










*/

// vector<dim> class declaration
namespace ms::math
{

template <int Dim = 0>
class Vector : public Vector_Wrap
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
class Vector<0> : public Vector_Wrap
{
public:
  Vector(void) = default;
  explicit Vector(const int dimension); // prevent implicit convsersion to int
  Vector(const std::initializer_list<double> list);
  Vector(const Vector_View vector_view);
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

// free function decalarations
namespace ms::math
{
Vector<0> operator*(const double constant, const Vector_View x);
}

/*










*/

// template definitions
namespace ms::math
{

template <ms::math::const_span T>
Vector_View::Vector_View(const T& values, const int inc)
    : _values_view(values),
      _inc(inc) { this->initialize(); };

template <std::contiguous_iterator T, ms::math::contiguous_iterator_or_int U>
Vector_View::Vector_View(const T& t, const U& u, const int inc)
    : _values_view(t, u),
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
  this->_values_view = this->_coordinates;
  this->_dimension   = static_cast<int>(this->_coordinates.size());
  this->_inc         = 1;
  this->_values_wrap = this->_coordinates;
}

template <std::input_iterator Iter>
Vector<0>::Vector(Iter first, Iter last)
    : _coordinates(first, last)
{
  this->reallocate_values();
};

} // namespace ms::math