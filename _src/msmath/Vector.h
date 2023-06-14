#pragma once
#include <array>
#include <string>
#include <vector>

#define COMPILE_TIME_REQUIREMENT static_assert

// miscellaneous definitions
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
  operator std::pair<const double*, int>(void) const;

public:
  double               at(const size_t position) const;
  const double*        begin(void) const;
  double               cosine(const Vector_Const_Wrapper& other) const;
  Vector_Const_Wrapper const_wrapper(void) const;
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
  // Acess to the default constructor has been blocked, because it creates objects that do not function properly.
  // The reason why it was not deleted is that the constructor of the Vector<0> class uses the default constructor.
  Vector_Const_Wrapper(void) = default;

protected:
  int           _dimension            = 0;
  const double* _coordinate_const_ptr = nullptr;
  int           _inc                  = 0;
};

/*










*/

class Vector_Wrapper : public Vector_Const_Wrapper
{
public:
  Vector_Wrapper(double* coordinate_ptr, const int dimension, const int incx = 1)
      : Vector_Const_Wrapper(coordinate_ptr, dimension, incx),
        _coordinate_ptr(coordinate_ptr){};
  Vector_Wrapper(std::vector<double>& coordinates, const int incx = 1)
      : Vector_Const_Wrapper(coordinates, incx),
        _coordinate_ptr(coordinates.data()){};
  template <int dim>
  Vector_Wrapper(const Vector<dim>& vec) = delete;

public:
  void    operator*=(const double constant);
  void    operator+=(const Vector_Const_Wrapper& other);
  void    operator-=(const Vector_Const_Wrapper& other);
  double& operator[](const size_t position);

public:
  operator double*(void);

public:
  double&        at(const size_t position);
  void           change_value(const Vector_Const_Wrapper other);
  double*        data(void);
  void           initalize(void);
  void           normalize(void);
  Vector_Wrapper wrapper(void);

  // resolve base class method hiding problem (overloading across scope problem)
public:
  using Vector_Const_Wrapper::operator[];
  using Vector_Const_Wrapper::at;
  using Vector_Const_Wrapper::data;

protected:
  // The default constructor, which creates objects that do not function properly, has been blocked from access.
  // The reason why it was not deleted is that the constructor of the Vector<0> class uses the default constructor.
  Vector_Wrapper(void) = default;

protected:
  double* _coordinate_ptr = nullptr;
};

/*










*/

template <int Dim = 0>
class Vector : public Vector_Wrapper
{
  COMPILE_TIME_REQUIREMENT(0 <= Dim, "dimension shold be positive");

public:
  Vector(void)
      : Vector_Wrapper(this->coordinates_.data(), Dim){};
  Vector(const Vector& other)
      : Vector_Wrapper(this->coordinates_.data(), Dim),
        coordinates_{other.coordinates_} {};

  template <typename... Args>
  Vector(const Args&... args)
      : Vector_Wrapper(this->coordinates_.data(), Dim),
        coordinates_{static_cast<double>(args)...}
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

/*










*/

template <>
class Vector<0> : public Vector_Wrapper
{
  // The reason for using the default constructor of the Wrapper class is as follows:
  // The parent class, Wrapper class's constructor must be called before the coordinate is initialized.
  // If the coordinate is used as is, it will trigger the nullptr checking code in the Wrapper class constructor.
  // Therefore, instead of removing the nullptr checking code in the Wrapper class constructor,
  // we chose to use the inaccessible default constructor to solve the problem."

public:
  static Vector<0> null_vector(void);
  explicit Vector(const int dimension);
  Vector(const std::initializer_list<double> list);
  Vector(const std::vector<double>& values);
  Vector(std::vector<double>&& values);
  Vector(const Vector& other);
  Vector(Vector&& other) noexcept;

public:
  template <typename Iter, ms::math::if_is_iterator<Iter> = true>
  Vector(Iter first, Iter last)
      : _coordinates(first, last)
  {
    this->reallocate_ptr();
    this->_dimension = static_cast<int>(this->_coordinates.size());
    this->_inc                  = 1;
  };

public:
  void operator=(const Vector& other);
  void operator=(Vector&& other) noexcept;

public:
  void resize(const int size);

  private:
  void reallocate_ptr(void);

private:
  std::vector<double> _coordinates;

  Vector(void) = default;
};

} // namespace ms::math

/*





*/

// user-defined deduction guides
namespace ms::math
{
template <typename Iter, ms::math::if_is_iterator<Iter> = true>
Vector(Iter iter1, Iter iter2) -> Vector<0>;

template <typename... Args, ms::math::if_is_more_than_two<Args...> = true>
Vector(Args... args) -> Vector<sizeof...(Args)>;
} // namespace ms::math

/*





*/

// free function decalarations
namespace ms::math
{
Vector<0> operator*(const double constant, const Vector_Const_Wrapper& x);
}
