#include "Node.h"

#include "msexception/Exception.h"
#include <algorithm>

namespace ms::geo
{

Node_Const_Wrapper::Node_Const_Wrapper(const int n, const double* coords_ptr, const int stride)
    : _dimension(n),
      _coords_ptr(coords_ptr),
      _stride(stride)
{
  REQUIRE(0 < this->_dimension, "dimension should be natural number");
  REQUIRE(this->_coords_ptr != nullptr, "coords ptr should not be null ptr");
  REQUIRE(0 < this->_stride, "dimension should be natural number");
};

double Node_Const_Wrapper::operator[](const int index) const
{
  REQUIRE(index < this->_dimension, "index can not exceed dimension");
  return this->_coords_ptr[index * this->_stride];
}

std::vector<double> Node_Const_Wrapper::operator-(const Node_Const_Wrapper& other) const
{
  std::vector<double> result(this->_dimension);
  for (int i = 0; i < this->_dimension; ++i)
  {
    result[i] = this->_coords_ptr[i] - other._coords_ptr[i];
  }

  return result;
}

Node_Const_Wrapper::operator const double*(void) const
{
  return this->_coords_ptr;
}

int Node_Const_Wrapper::dimension(void) const
{
  return this->_dimension;
}

int Node_Const_Wrapper::size(void) const
{
  return this->_dimension;
}

Nodes_Const_Wrapper::Nodes_Const_Wrapper(const Coordinates_Type type, const int num_nodes, const int dimension, const double* coords_ptr)
    : _num_nodes(num_nodes),
      _dimension(dimension),
      _coordinates_ptr(coords_ptr)
{
  REQUIRE(0 < this->_num_nodes, "number of nodes should be natural number");
  REQUIRE(0 < this->_dimension, "dimension should be natural number");
  REQUIRE(this->_coordinates_ptr != nullptr, "coords ptr should not be null ptr");

  if (type == Coordinates_Type::NODAL)
  {
    this->_stride = 1;
  }
  else if (type == Coordinates_Type::BLOCK)
  {
    this->_stride = this->_num_nodes;
  }
  else
  {
    EXCEPTION("not supported nodes type");
  }
};

Node_Const_Wrapper Nodes_Const_Wrapper::operator[](const int index) const
{
  // convention : index start with 0
  REQUIRE(index < this->_num_nodes, "index can not exceed number of nodes");

  const auto coordinate_ptr = this->coordinate_ptr(index);
  return {this->_dimension, coordinate_ptr, this->_stride};
}

void Nodes_Const_Wrapper::copy_coordinates(double* dest) const
{
  std::copy(this->_coordinates_ptr, this->_coordinates_ptr + this->num_coordinates(), dest);
}

const double* Nodes_Const_Wrapper::coordinates_ptr(void) const
{
  return this->_coordinates_ptr;
}

int Nodes_Const_Wrapper::dimension(void) const
{
  return this->_dimension;
}

int Nodes_Const_Wrapper::num_coordinates(void) const
{
  return this->_dimension * this->_num_nodes;
}

int Nodes_Const_Wrapper::num_nodes(void) const
{
  return this->_num_nodes;
}

std::vector<Node_Const_Wrapper> Nodes_Const_Wrapper::nodes_at_indexes(const std::vector<int>& indexes) const
{
  const auto num_nodes = indexes.size();
  REQUIRE(num_nodes < this->_num_nodes, "number of indexes can not excced number of nodes");

  std::vector<Node_Const_Wrapper> result;
  result.reserve(num_nodes);

  for (const auto index : indexes)
  {
    result.push_back((*this)[index]);
  }

  return result;
}


const double* Nodes_Const_Wrapper::coordinate_ptr(const int index) const
{
  // convention : index start with 0
  REQUIRE(index < this->_num_nodes, "index can not exceed number of nodes");

  auto coordinate_start_index = index;

  if (this->is_nodal())
  {
    coordinate_start_index *= this->_dimension;
  }

  return this->_coordinates_ptr + coordinate_start_index;
}

bool Nodes_Const_Wrapper::is_nodal(void) const
{
  return this->_stride == 1;
}


Nodes::Nodes(const Coordinates_Type type, const int num_nodes, const int dimension, std::vector<double>&& coordinates)
    : Nodes_Const_Wrapper(type, num_nodes, dimension, nullptr),
      _coordinates(std::move(coordinates))
{
  this->_coordinates_ptr = this->_coordinates.data();
}

} // namespace ms::geo
