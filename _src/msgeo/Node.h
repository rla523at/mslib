#pragma once
#include <vector>

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
  Node_Const_Wrapper(const int dimension, const double* coords_ptr, const int stride = 1);

public:
  std::vector<double> operator-(const Node_Const_Wrapper& other) const;
  double              operator[](const int index) const;

public:
  operator const double*(void) const;

public:
  int  dimension(void) const;
  void other_to_this_vector(const Node_Const_Wrapper& other, double* vector_components) const;

  // When a point is viewed as tuple, return a size of tuple
  int size(void) const;

private:
  int           _dimension  = 0;
  const double* _coords_ptr = nullptr;
  int           _stride     = 0;
};

class Nodes_Const_Wrapper
{
public:
  Nodes_Const_Wrapper(const Coordinates_Type type, const int num_nodes, const int dimension, const double* coords_ptr);

public:
  Node_Const_Wrapper operator[](const int index) const;

public:
  void                            copy_coordinates(double* dest) const;
  const double*                   coordinates_ptr(void) const;
  int                             dimension(void) const;
  int                             num_coordinates(void) const;
  int                             num_nodes(void) const;
  std::vector<Node_Const_Wrapper> nodes_at_indexes(const std::vector<int>& indexes) const;
  void                            nodes_at_indexes(std::vector<Node_Const_Wrapper>& result, const std::vector<int>& indexes) const;

private:
  const double* coordinate_ptr(const int index) const;
  bool          is_nodal(void) const;

protected:
  int           _num_nodes       = 0;
  int           _dimension       = 0;
  const double* _coordinates_ptr = nullptr;
  int           _stride          = 0;
};

class Nodes : public Nodes_Const_Wrapper
{
public:
  Nodes(const Coordinates_Type type, const int num_nodes, const int dimension, std::vector<double>&& coordinates);

private:
  std::vector<double> _coordinates;
};

} // namespace ms::geo
