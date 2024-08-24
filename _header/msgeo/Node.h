#pragma once
#include <initializer_list>
#include <span>
#include <string>
#include <vector>

#include "msmath/BLAS.h"
#include "msmath/Vector.h"

// NODAL coords format : (X1,Y1,Z1) ... (Xn,Yn, Zn)
// BLOCK coords format : (X1,...,Xn) ... (Z1,...,Zn)

// miscellaneous definitions
namespace ms::geo
{

template <typename... T>
concept can_form_const_span = requires(T&&... args) {
  std::span<const double>(std::forward<T>(args)...);
};

template <typename T>
concept span = requires(T& t) { std::span<double>(t); };

} // namespace ms::geo

/*










*/

// class declaration
namespace ms::geo
{

class Node_View
{
public:
  Node_View(void) = default;

  template <typename... T>
    requires can_form_const_span<T...>
  Node_View(T&&... coordinates);

  template <typename... T>
    requires can_form_const_span<T...>
  Node_View(const int inc = 1, T&&... coordinates);

public:
  ms::math::Vector<0> operator-(const Node_View& other) const;
  double              operator[](const int index) const;
  bool                operator==(const Node_View& other) const;

public:
  double                at(const int index) const;
  int                   dimension(void) const;
  double                distance(const Node_View other) const;
  void                  other_to_this_vector(const Node_View& other, ms::math::Vector_Wrap vector_wrap) const;
  int                   size(void) const; // When a point is viewed as tuple, return a size of tuple
  std::string           to_string(void) const;
  ms::math::Vector_View to_vector_view(void) const;

protected:
  const double* data(void) const;

protected:
  std::span<const double> _coordinates_view;
  int                     _dimension = 0;
  int                     _inc       = 0;
};

} // namespace ms::geo

/*










*/

namespace ms::geo
{

class Node_Wrap : public Node_View
{
public:
  Node_Wrap(void) = default;

  template <span T>
  Node_Wrap(T& coordinates)
    : Node_View(coordinates),
      _coordinates_wrap(coordinates){};

  template <span T>
  Node_Wrap(T& coordsinates, const int dimension, const int inc = 1)
    : Node_View(coordsinates, dimension, inc),
      _coordinates_wrap(coordsinates){};

public:
  double& operator[](const int index);

public:
  double&               at(const int index);
  ms::math::Vector_Wrap to_vector_wrap(void);

  // resolve base class method hiding problem (overloading across scope problem)
public:
  using Node_View::at;
  using Node_View::operator[];

protected:
  std::span<double> _coordinates_wrap;
};

} // namespace ms::geo

/*










*/

namespace ms::geo
{

class Node : public Node_Wrap
{
public:
  explicit Node(const int dim);
  Node(std::initializer_list<double> coordinates);
  Node(const Node& other);
  Node(Node&& other) noexcept;

  template <typename... T>
    requires can_form_const_span<T...>
  Node(T&&... coordinates);

private:
  void reallocate_coordinates(void);

private:
  std::array<double, 3> _coordinates = {0.0};
};

} // namespace ms::geo

/*










*/

namespace ms::geo
{

struct Numbered_Node_View
{
  Node_View node_view;
  int       number;
};

struct Numbered_Node
{
  Node node;
  int  number;
};

// return ture if lhs < rhs
struct Node_Compare
{
  bool operator()(Node_View lhs, Node_View rhs) const;
};

} // namespace ms::geo

/*










*/

namespace ms::geo
{

class Nodes_View
{
public:
  Nodes_View(void) = default;

  template <typename... T>
    requires can_form_const_span<T...>
  Nodes_View(const int num_nodes, const int dimension, T&&... coordinates);

public:
  Node_View operator[](const int index) const;
  bool      operator==(const Nodes_View other) const;

public:
  Node_View              at(const int index) const;
  ms::math::Vector_View  axis_coordinates_vector_view(const int axis_index) const;
  void                   copy_coordinates(double* dest) const;
  int                    dimension(void) const;
  int                    num_coordinates(void) const;
  int                    num_nodes(void) const;
  std::vector<Node_View> nodes_at_indexes(const std::vector<int>& numbers) const;
  void                   nodes_at_indexes(std::vector<Node_View>& result, const std::vector<int>& numbers) const;

protected:
  std::span<const double> _coordinates_view;
  int                     _num_nodes = 0;
  int                     _dimension = 0;
};

} // namespace ms::geo

/*










*/

namespace ms::geo
{

class Nodes_Wrapper : public Nodes_View
{
public:
  Nodes_Wrapper(void) = default;

  template <span T>
  Nodes_Wrapper(T& coordinates, const int num_nodes, const int dimension)
    : Nodes_View(coordinates, num_nodes, dimension),
      _coordinates_wrap(coordinates){};

public:
  Node_Wrap operator[](const int index);

public:
  Node_Wrap at(const int index);

  // resolve base class method hiding problem (overloading across scope problem)
public:
  using Nodes_View::operator[];
  using Nodes_View::at;

protected:
  std::span<double> _coordinates_wrap;
};

} // namespace ms::geo

/*










*/

namespace ms::geo
{

class Nodes : public Nodes_Wrapper
{
public:
  Nodes(void) = default;
  Nodes(const int num_nodes, const int dimension);
  Nodes(const int num_nodes, const int dimension, std::vector<double>&& coordinates);
  Nodes(const Nodes& other);
  Nodes(Nodes&& other) noexcept;

public:
  void operator=(const Nodes& other);
  void operator=(Nodes&& other);

public:
  void add_node(const Node_View node);
  void add_nodes(Nodes&& other);

private:
  bool is_null_node(void) const;

private:
  std::vector<double> _coordinates; // it always be nodal coordinate type
};

} // namespace ms::geo

/*










*/

// template definitions
namespace ms::geo
{

template <typename... T>
  requires can_form_const_span<T...>
inline Node_View::Node_View(T&&... args)
  : _coordinates_view(std::forward<T>(args)...)
{
  this->_dimension = static_cast<int>(this->_coordinates_view.size());
  this->_inc       = 1;
}

template <typename... T>
  requires can_form_const_span<T...>
inline Node_View::Node_View(const int inc, T&&... coordinates)
  : _coordinates_view(std::forward<T>(coordinates)...)
{
  REQUIRE(0 < inc, "inc should be natural number");

  this->_dimension = static_cast<int>(this->_coordinates_view.size());
  this->_inc       = inc;
}

template <typename... T>
  requires can_form_const_span<T...>
inline Nodes_View::Nodes_View(const int num_nodes, const int dimension, T&&... coordinates)
  : _coordinates_view(std::forward<T>(coordinates)...),
    _num_nodes(num_nodes),
    _dimension(dimension)
{
  REQUIRE(0 < this->_num_nodes, "number of nodes should be natural number");
  REQUIRE(0 < this->_dimension, "dimension should be natural number");
  REQUIRE(this->_coordinates_view.size() == this->_num_nodes * this->_dimension, "Given coordinates size should be matched with given size");
};

template <typename... T>
  requires can_form_const_span<T...>
Node::Node(T&&... coordinates)
{
  std::span<const double> cooridnates_view(std::forward<T>(coordinates)...);

  this->_dimension = static_cast<int>(cooridnates_view.size());
  this->_inc       = 1;
  REQUIRE(0 < this->_dimension, "dimension should be natural number");
  REQUIRE(this->_dimension <= 3, "dimension can't not excced 3");

  ms::math::blas::copy(this->_coordinates.data(), cooridnates_view.data(), this->_dimension);
  this->reallocate_coordinates();
}

} // namespace ms::geo

/*










*/

namespace ms::geo
{
void A_to_B_vector(ms::math::Vector_Wrap vector_wrap, const Node_View A, const Node_View B);
}