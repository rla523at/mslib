#pragma once
#include "Data.h"
#include "msgeo/Geometry.h"
#include <vector>

namespace ms::grid
{

class Element
{
public:
  Element(const Element_Type element_type, std::vector<int>&& node_indexes, ms::geo::Geometry&& geometry)
      : _type(element_type),
        _consisting_node_indexes(std::move(node_indexes)),
        _geometry(std::move(geometry)){};

public:
  std::vector<int> find_periodic_matched_node_indexes(const Element& other) const;

private:
  bool can_be_periodic_pair(const Element& other) const;

private:
  ms::geo::Geometry _geometry;
  std::vector<int>  _consisting_node_indexes;
  Element_Type      _type;
};

} // namespace ms::grid