#pragma once
#include "Data.h"
#include "msgeo/Geometry.h"
#include "msgeo/Node.h"
#include <vector>

// data definition
namespace ms::grid
{

struct Indexed_Node
{
  int                         index;
  ms::geo::Node_Const_Wrapper node;
};

} // namespace ms::grid

/*





*/

// forward declaration
namespace ms::math
{
class Vector_Const_Wrapper;
}

/*





*/

// class declaration
namespace ms::grid
{

class Element
{
public:
  Element(const Element_Type element_type, std::vector<Indexed_Node>&& consisting_nodes, ms::geo::Geometry&& geometry)
      : _type(element_type),
        _consisting_indexed_nodes(std::move(consisting_nodes)),
        _geometry(std::move(geometry)){};

public:
  std::vector<int> find_periodic_matched_node_indexes(const ms::math::Vector_Const_Wrapper& direction_vector, const Element& other) const;

private:
  int dimension(void) const;

private:
  Element_Type              _type;
  std::vector<Indexed_Node> _consisting_indexed_nodes;
  ms::geo::Geometry         _geometry;
};

} // namespace ms::grid