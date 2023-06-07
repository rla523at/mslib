#pragma once
#include <vector>
#include <string>

// data definition
namespace ms::geo
{
// NODAL coords format : (X1,Y1,Z1) ... (Xn,Yn, Zn)
// BLOCK coords format : (X1,...,Xn) ... (Z1,...,Zn)
enum class Coordinates_Type
{
  NODAL,
  BLOCK,
  NOT_SUPPROTED
};
} // namespace ms::geo

/*










*/

// class declaration
namespace ms::geo
{

class Node_Const_Wrapper
{
public:
  Node_Const_Wrapper(const int dimension, const double* coords_ptr, const int inc = 1);

public:
  std::vector<double> operator-(const Node_Const_Wrapper& other) const;
  double              operator[](const int index) const;
  bool                operator==(const Node_Const_Wrapper& other) const;

public:
  operator std::pair<const double*, int>(void) const;

public:
  int         dimension(void) const;
  void        other_to_this_vector(const Node_Const_Wrapper& other, double* vector_components) const;
  int         size(void) const; // When a point is viewed as tuple, return a size of tuple
  std::string to_string(void) const;

protected:
  int           _dimension             = 0;
  const double* _coordinates_const_ptr = nullptr;
  int           _inc                   = 0;

protected:
  // Acess to the default constructor has been blocked, because it creates objects that do not function properly.
  // The reason why it was not deleted is that the constructor of the Node class uses the default constructor.
  Node_Const_Wrapper(void) = default;
};

/*










*/

class Node_Wrapper : public Node_Const_Wrapper
{
public:
  Node_Wrapper(const int dimension, double* coords_ptr, const int inc = 1)
      : Node_Const_Wrapper(dimension, coords_ptr, inc),
        _coordinates_ptr(coords_ptr){};

public:
  double& operator[](const int position);

public:
  operator std::pair<double*, int>(void) const;

  // resolve base class method hiding problem (overloading across scope problem)
public:
  using Node_Const_Wrapper::operator[];

protected:
  double* _coordinates_ptr = nullptr;

protected:
  // Acess to the default constructor has been blocked, because it creates objects that do not function properly.
  // The reason why it was not deleted is that the constructor of the Node class uses the default constructor.
  Node_Wrapper(void) = default;
};

/*










*/
// class Nodes;

class Node : public Node_Wrapper
{
public:
  Node(const int dim);
  Node(std::vector<double>&& coordinates);
  Node(const Node& other);
  Node(Node&& other) noexcept;

public:
  operator Node_Const_Wrapper(void) const;

private:
  std::vector<double> _coordinates;

  // friend Nodes;
};

/*










*/

class Nodes_Const_Wrapper
{
public:
  Nodes_Const_Wrapper(const Coordinates_Type type, const int num_nodes, const int dimension, const double* coords_ptr);

public:
  Node_Const_Wrapper operator[](const int index) const;

public:
  Node_Const_Wrapper              at(const int index) const;
  void                            copy_coordinates(double* dest) const;
  const double*                   coordinates_ptr(void) const;
  Coordinates_Type                coordinates_type(void) const;
  int                             dimension(void) const;
  int                             num_coordinates(void) const;
  int                             num_nodes(void) const;
  std::vector<Node_Const_Wrapper> nodes_at_indexes(const std::vector<int>& numbers) const;
  void                            nodes_at_indexes(std::vector<Node_Const_Wrapper>& result, const std::vector<int>& numbers) const;

protected:
  bool is_nodal(void) const;

private:
  const double* coordinates_ptr_at(const int index) const;

protected:
  int           _num_nodes             = 0;
  int           _dimension             = 0;
  const double* _coordinates_const_ptr = nullptr;
  int           _inc                   = 0;

protected:
  // Acess to the default constructor has been blocked, because it creates objects that do not function properly.
  // The reason why it was not deleted is that the constructor of the Nodes class uses the default constructor.
  Nodes_Const_Wrapper(void) = default;
};

/*










*/

class Nodes_Wrapper : public Nodes_Const_Wrapper
{
public:
  Nodes_Wrapper(const Coordinates_Type type, const int num_nodes, const int dimension, double* coords_ptr)
      : Nodes_Const_Wrapper(type, num_nodes, dimension, coords_ptr),
        _coordinates_ptr(coords_ptr){};

public:
  Node_Wrapper operator[](const int index);

public:
  Node_Wrapper at(const int index);

private:
  double* coordinates_ptr_at(const int index);

  // resolve base class method hiding problem (overloading across scope problem)
public:
  using Nodes_Const_Wrapper::operator[];
  using Nodes_Const_Wrapper::at;

protected:
  double* _coordinates_ptr = nullptr;

protected:
  // Acess to the default constructor has been blocked, because it creates objects that do not function properly.
  // The reason why it was not deleted is that the constructor of the Node class uses the default constructor.
  Nodes_Wrapper(void) = default;
};

class Nodes : public Nodes_Wrapper
{
public:
  static Nodes Null_Nodes(void);

public:
  Nodes(const Coordinates_Type type, const int num_nodes, const int dimension);
  Nodes(const Coordinates_Type type, const int num_nodes, const int dimension, std::vector<double>&& coordinates);
  Nodes(const Nodes& other);
  Nodes(Nodes&& other) noexcept;

public:
  void operator=(const Nodes& other);
  void operator=(Nodes&& other);

public:
  operator Nodes_Const_Wrapper(void) const;

public:
  void add_node(const Node_Const_Wrapper node);
  void add_nodes(Nodes&& other);
  void to_block_type(void);

private:
  bool is_null_node(void) const;

private:
  std::vector<double> _coordinates;

private:
  Nodes(void) = default;
};

} // namespace ms::geo
