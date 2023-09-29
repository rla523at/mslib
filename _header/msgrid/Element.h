#pragma once
#include "Data.h"
#include "msgeo/Geometry.h"
#include "msgeo/Node.h"
#include "msmath/Vector.h"
#include <vector>

// class declaration
namespace ms::grid
{

class Element
{
public:
  // node_numbers는 geometry의 node_views와 일치해야한다.
  Element(const Element_Type type, std::vector<int>&& node_numbers, ms::geo::Geometry&& geometry)
      : _type(type),
        _node_numbers(std::move(node_numbers)),
        _geometry(std::move(geometry)){};

public:
  void reordering_nodes(const std::vector<int>& new_ordered_node_indexes);

public:
  int                           dimension(void) const;
  std::vector<int>              find_periodic_matched_node_numbers(const ms::math::Vector_View direction_vector, const Element& other) const;
  std::vector<std::vector<int>> face_vertex_node_numbers_s(void) const;
  void                          face_vertex_node_numbers_s(std::vector<int>* face_index_to_face_vnode_numbers) const;
  const ms::geo::Geometry&      get_geometry(void) const;
  std::span<const int>          node_numbers(void) const;
  bool                          is_outward_face(const Element& face_element) const;
  Element                       make_face_element(const int face_index) const;
  Element_Type                  type(void) const;
  std::span<const int>          vertex_node_numbers(void) const;

private:
  Element_Type      _type;
  std::vector<int>  _node_numbers;
  ms::geo::Geometry _geometry;
};

} // namespace ms::grid