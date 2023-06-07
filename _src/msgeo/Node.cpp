#include "Node.h"

#include "msexception/Exception.h"
#include <algorithm>
#include <iomanip>
#include <sstream>

namespace ms::geo
{

Node_Const_Wrapper::Node_Const_Wrapper(const int n, const double* coords_ptr, const int stride)
    : _dimension(n),
      _coordinates_const_ptr(coords_ptr),
      _inc(stride)
{
  REQUIRE(0 < this->_dimension, "dimension should be natural number");
  REQUIRE(this->_coordinates_const_ptr != nullptr, "coords ptr should not be null ptr");
  REQUIRE(0 < this->_inc, "dimension should be natural number");
};

double Node_Const_Wrapper::operator[](const int index) const
{
  REQUIRE(index < this->_dimension, "index can not exceed dimension");
  return this->_coordinates_const_ptr[index * this->_inc];
}

bool Node_Const_Wrapper::operator==(const Node_Const_Wrapper& other) const
{
  if (this->_dimension != other._dimension) return false;

  for (int i = 0; i < this->_dimension; ++i)
  {
    if ((*this)[i] != other[i]) return false;
  }

  return true;
}

std::vector<double> Node_Const_Wrapper::operator-(const Node_Const_Wrapper& other) const
{
  std::vector<double> result(this->_dimension);
  for (int i = 0; i < this->_dimension; ++i)
  {
    result[i] = this->_coordinates_const_ptr[i] - other._coordinates_const_ptr[i];
  }

  return result;
}

Node_Const_Wrapper::operator std::pair<const double*, int>(void) const
{
  return {this->_coordinates_const_ptr, this->_inc};
}

int Node_Const_Wrapper::dimension(void) const
{
  return this->_dimension;
}

void Node_Const_Wrapper::other_to_this_vector(const Node_Const_Wrapper& other, double* vector_components) const
{
  for (int i = 0; i < this->_dimension; ++i)
  {
    vector_components[i] = (*this)[i] - other[i];
  }
}

int Node_Const_Wrapper::size(void) const
{
  return this->_dimension;
}

std::string Node_Const_Wrapper::to_string(void) const
{
  std::ostringstream oss;
  oss << std::setprecision(16) << std::showpoint << std::left;
  for (int i = 0; i < this->_dimension; ++i)
  {
    oss << std::setw(25) << this->_coordinates_const_ptr[i];
  }
  return oss.str();
}

/*










*/

double& Node_Wrapper::operator[](const int index)
{
  REQUIRE(index < this->_dimension, "index can not exceed dimension");
  return this->_coordinates_ptr[index * this->_inc];
}

Node_Wrapper::operator std::pair<double*, int>(void) const
{
  return {this->_coordinates_ptr, this->_inc};
}

/*










*/

Node::Node(const int dim)
    : _coordinates(dim)
{
  this->_dimension             = static_cast<int>(this->_coordinates.size());
  this->_coordinates_const_ptr = this->_coordinates.data();
  this->_coordinates_ptr       = this->_coordinates.data();
  this->_inc                   = 1;
  REQUIRE(0 < this->_dimension, "dimension should be natural number");
  REQUIRE(this->_coordinates_const_ptr != nullptr, "coords ptr should not be null ptr");
  REQUIRE(this->_coordinates_ptr != nullptr, "coords ptr should not be null ptr");
}

Node::Node(std::vector<double>&& coordinates)
    : _coordinates(std::move(coordinates))
{
  this->_dimension             = static_cast<int>(this->_coordinates.size());
  this->_coordinates_const_ptr = this->_coordinates.data();
  this->_coordinates_ptr       = this->_coordinates.data();
  this->_inc                   = 1;
  REQUIRE(0 < this->_dimension, "dimension should be natural number");
  REQUIRE(this->_coordinates_const_ptr != nullptr, "coords ptr should not be null ptr");
  REQUIRE(this->_coordinates_ptr != nullptr, "coords ptr should not be null ptr");
}

Node::Node(const Node& other)
    : _coordinates(other._coordinates)
{
  this->_coordinates_const_ptr = this->_coordinates.data();
  this->_coordinates_ptr       = this->_coordinates.data();
  this->_dimension             = other._dimension;
  this->_inc                   = other._inc;
}

Node::Node(Node&& other) noexcept
    : _coordinates(std::move(other._coordinates))
{
  this->_coordinates_const_ptr = this->_coordinates.data();
  this->_coordinates_ptr       = this->_coordinates.data();
  this->_dimension             = other._dimension;
  this->_inc                   = other._inc;

  other._coordinates_const_ptr = nullptr;
  other._coordinates_ptr       = nullptr;
}

Node::operator Node_Const_Wrapper(void) const
{
  return {this->_dimension, this->_coordinates_const_ptr, this->_inc};
}

/*










*/

Nodes_Const_Wrapper::Nodes_Const_Wrapper(const Coordinates_Type type, const int num_nodes, const int dimension, const double* coords_ptr)
    : _num_nodes(num_nodes),
      _dimension(dimension),
      _coordinates_const_ptr(coords_ptr)
{
  REQUIRE(0 < this->_num_nodes, "number of nodes should be natural number");
  REQUIRE(0 < this->_dimension, "dimension should be natural number");
  REQUIRE(this->_coordinates_const_ptr != nullptr, "coords ptr should not be null ptr");

  if (type == Coordinates_Type::NODAL)
  {
    this->_inc = 1;
  }
  else if (type == Coordinates_Type::BLOCK)
  {
    this->_inc = this->_num_nodes;
  }
  else
  {
    EXCEPTION("not supported nodes type");
  }
};

Node_Const_Wrapper Nodes_Const_Wrapper::operator[](const int index) const
{
  return this->at(index);
}

Node_Const_Wrapper Nodes_Const_Wrapper::at(const int index) const
{
  // convention : index start with 0
  REQUIRE(index < this->_num_nodes, "index can not exceed number of nodes");

  const auto coordinate_ptr = this->coordinates_ptr_at(index);
  return {this->_dimension, coordinate_ptr, this->_inc};
}

void Nodes_Const_Wrapper::copy_coordinates(double* dest) const
{
  std::copy(this->_coordinates_const_ptr, this->_coordinates_const_ptr + this->num_coordinates(), dest);
}

const double* Nodes_Const_Wrapper::coordinates_ptr(void) const
{
  return this->_coordinates_const_ptr;
}

Coordinates_Type Nodes_Const_Wrapper::coordinates_type(void) const
{
  if (this->is_nodal()) return Coordinates_Type::NODAL;

  return Coordinates_Type::BLOCK;
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

std::vector<Node_Const_Wrapper> Nodes_Const_Wrapper::nodes_at_indexes(const std::vector<int>& numbers) const
{
  const auto num_nodes = numbers.size();
  REQUIRE(num_nodes < this->_num_nodes, "number of numbers can not excced number of nodes");

  std::vector<Node_Const_Wrapper> result;
  result.reserve(num_nodes);

  for (const auto index : numbers)
  {
    result.push_back((*this)[index]);
  }

  return result;
}

void Nodes_Const_Wrapper::nodes_at_indexes(std::vector<Node_Const_Wrapper>& result, const std::vector<int>& numbers) const
{
  const auto num_nodes = numbers.size();
  REQUIRE(num_nodes < this->_num_nodes, "number of numbers can not excced number of nodes");

  result.clear();
  result.reserve(num_nodes);

  for (const auto index : numbers)
  {
    result.push_back((*this)[index]);
  }
}

const double* Nodes_Const_Wrapper::coordinates_ptr_at(const int index) const
{
  // convention : index start with 0
  REQUIRE(index < this->_num_nodes, "index can not exceed number of nodes");

  auto coordinate_start_index = index;

  if (this->is_nodal())
  {
    coordinate_start_index *= this->_dimension;
  }

  return this->_coordinates_const_ptr + coordinate_start_index;
}

bool Nodes_Const_Wrapper::is_nodal(void) const
{
  return this->_inc == 1;
}

/*










*/

Node_Wrapper Nodes_Wrapper::operator[](const int index)
{
  return this->at(index);
}

Node_Wrapper Nodes_Wrapper::at(const int index)
{
  REQUIRE(this->_coordinates_const_ptr != nullptr, "coordinate ptr should not be nullptr");
  REQUIRE(this->_coordinates_ptr != nullptr, "coordinate ptr should not be nullptr");
  REQUIRE(index < this->_num_nodes, "index can not exceed number of nodes");

  const auto coordinate_ptr = this->coordinates_ptr_at(index);
  return {this->_dimension, coordinate_ptr, this->_inc};
}

double* Nodes_Wrapper::coordinates_ptr_at(const int index)
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

/*










*/

Nodes Nodes::Null_Nodes(void)
{
  return Nodes();
}

Nodes::Nodes(const Coordinates_Type type, const int num_nodes, const int dimension)
{
  this->_num_nodes = num_nodes;
  this->_dimension = dimension;
  REQUIRE(0 < this->_num_nodes, "number of nodes should be natural number");
  REQUIRE(0 < this->_dimension, "dimension should be natural number");

  this->_coordinates.resize(this->_num_nodes * this->_dimension);
  this->_coordinates_const_ptr = this->_coordinates.data();
  this->_coordinates_ptr       = this->_coordinates.data();
  REQUIRE(this->_coordinates_const_ptr != nullptr, "coords ptr should not be null ptr");
  REQUIRE(this->_coordinates_ptr != nullptr, "coords ptr should not be null ptr");

  if (type == Coordinates_Type::NODAL)
  {
    this->_inc = 1;
  }
  else if (type == Coordinates_Type::BLOCK)
  {
    this->_inc = this->_num_nodes;
  }
  else
  {
    EXCEPTION("not supported coordinates type");
  }
}

Nodes::Nodes(const Coordinates_Type type, const int num_nodes, const int dimension, std::vector<double>&& coordinates)
    : _coordinates(std::move(coordinates))
{
  this->_num_nodes             = num_nodes;
  this->_dimension             = dimension;
  this->_coordinates_const_ptr = this->_coordinates.data();
  REQUIRE(0 < this->_num_nodes, "number of nodes should be natural number");
  REQUIRE(0 < this->_dimension, "dimension should be natural number");
  REQUIRE(this->_coordinates_const_ptr != nullptr, "coords ptr should not be null ptr");

  if (type == Coordinates_Type::NODAL)
  {
    this->_inc = 1;
  }
  else if (type == Coordinates_Type::BLOCK)
  {
    this->_inc = this->_num_nodes;
  }
  else
  {
    EXCEPTION("not supported coordinates type");
  }
}

Nodes::Nodes(const Nodes& other)
    : _coordinates(other._coordinates)
{
  this->_coordinates_const_ptr = this->_coordinates.data();
  this->_coordinates_ptr       = this->_coordinates.data();
  this->_num_nodes             = other._num_nodes;
  this->_dimension             = other._dimension;
  this->_inc                   = other._inc;
}

Nodes::Nodes(Nodes&& other) noexcept
    : _coordinates(std::move(other._coordinates))
{
  this->_coordinates_const_ptr = this->_coordinates.data();
  this->_coordinates_ptr       = this->_coordinates.data();
  this->_num_nodes             = other._num_nodes;
  this->_dimension             = other._dimension;
  this->_inc                   = other._inc;

  other._coordinates_const_ptr = nullptr;
  other._coordinates_ptr       = nullptr;
}

void Nodes::operator=(const Nodes& other)
{
  this->_coordinates           = other._coordinates;
  this->_coordinates_const_ptr = this->_coordinates.data();
  this->_coordinates_ptr       = this->_coordinates.data();
  this->_num_nodes             = other._num_nodes;
  this->_dimension             = other._dimension;
  this->_inc                   = other._inc;
}

void Nodes::operator=(Nodes&& other)
{
  this->_coordinates           = std::move(other._coordinates);
  this->_coordinates_const_ptr = this->_coordinates.data();
  this->_coordinates_ptr       = this->_coordinates.data();
  this->_num_nodes             = other._num_nodes;
  this->_dimension             = other._dimension;
  this->_inc                   = other._inc;

  other._coordinates_const_ptr = nullptr;
  other._coordinates_ptr       = nullptr;
}

Nodes::operator Nodes_Const_Wrapper(void) const
{
  if (this->is_nodal())
  {
    return {Coordinates_Type::NODAL, this->_num_nodes, this->_dimension, this->_coordinates_const_ptr};
  }
  else
  {
    return {Coordinates_Type::BLOCK, this->_num_nodes, this->_dimension, this->_coordinates_const_ptr};
  }
}

void Nodes::add_node(const Node_Const_Wrapper node)
{
  const auto dim = node.dimension();

  if (this->is_null_node())
  {
    this->_coordinates.resize(dim);

    for (int i = 0; i < dim; ++i)
    {
      this->_coordinates[i] = node[i];
    }

    this->_coordinates_const_ptr = this->_coordinates.data();
    this->_coordinates_ptr       = this->_coordinates.data();
    this->_dimension             = dim;
    this->_inc                   = 1;
    this->_num_nodes++;
  }
  else
  {
    REQUIRE(this->_dimension == dim, "dimension should be matched for merging");
    const auto c1 = this->_dimension * this->_num_nodes;
    this->_coordinates.resize(c1 + dim);

    if (this->is_nodal())
    {
      for (int i = 0; i < dim; ++i)
      {
        this->_coordinates[c1 + i] = node[i];
      }
    }
    else
    {
      for (int i = 0; i < dim; ++i)
      {
        const auto loc_iter = this->_coordinates.begin() + (i + 1) * this->_num_nodes;

        this->_coordinates.insert(loc_iter, node[i]);
      }
    }

    this->_coordinates_const_ptr = this->_coordinates.data();
    this->_coordinates_ptr       = this->_coordinates.data();
    this->_num_nodes++;
  }
}

void Nodes::add_nodes(Nodes&& other)
{
  if (this->is_null_node())
  {
    this->_coordinates           = std::move(other._coordinates);
    this->_coordinates_const_ptr = this->_coordinates.data();
    this->_coordinates_ptr       = this->_coordinates.data();
    this->_num_nodes             = other._num_nodes;
    this->_dimension             = other._dimension;
    this->_inc                   = other._inc;

    other._coordinates_const_ptr = nullptr;
    other._coordinates_ptr       = nullptr;
  }
  else
  {
    REQUIRE(this->_dimension == other._dimension, "dimension should be matched for adding");
    REQUIRE(this->_inc == other._inc, "inc should be matched for adding");

    if (this->is_nodal())
    {
      this->_coordinates.insert(this->_coordinates.end(), other._coordinates.begin(), other._coordinates.end());
    }
    else
    {
      for (int i = 0; i < this->_dimension; ++i)
      {
        const auto loc_iter   = this->_coordinates.begin() + (i + 1) * this->_num_nodes + i * other._num_nodes;
        const auto start_iter = other._coordinates.begin() + i * other._num_nodes;
        const auto end_iter   = other._coordinates.begin() + (i + 1) * other._num_nodes;

        this->_coordinates.insert(loc_iter, start_iter, end_iter);
      }
    }

    this->_coordinates_const_ptr = this->_coordinates.data();
    this->_coordinates_ptr       = this->_coordinates.data();
    this->_num_nodes += other._num_nodes;

    other._coordinates.clear();
    other._coordinates_const_ptr = nullptr;
    other._coordinates_ptr       = nullptr;
  }
}

void Nodes::to_block_type(void)
{
  if (!this->is_nodal()) return;

  std::vector<double> block_coords(this->_coordinates.size());

  for (int i = 0; i < this->_num_nodes; ++i)
  {
    auto node = this->at(i);

    for (int j = 0; j < this->_dimension; ++j)
    {
      block_coords[j * this->_num_nodes + i] = node[j];
    }
  }
}

bool Nodes::is_null_node(void) const
{
  return this->_coordinates.empty();
}

} // namespace ms::geo
