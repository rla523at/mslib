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
  Element(const Element_Type element_type, std::vector<ms::geo::Numbered_Node_View>&& numbered_node_views, ms::geo::Geometry&& geometry)
      : _type(element_type),
        _numbered_node_views(std::move(numbered_node_views)),
        _geometry(std::move(geometry)){};

public:
  void reordering_nodes(const std::vector<int>& new_ordered_node_indexes);

public:
  void                          accumulate_discrete_node_info(ms::geo::Geometry_Consisting_Nodes_Info& partition_data, const int partition_order) const;
  void                          accumulate_node_info(ms::geo::Geometry_Consisting_Nodes_Info& partition_data, const int partition_order) const;
  int                           dimension(void) const;
  std::vector<int>              find_periodic_matched_node_indexes(const ms::math::Vector_View direction_vector, const Element& other) const;
  std::vector<std::vector<int>> face_index_to_face_vertex_node_numbers(void) const;
  void                          face_index_to_face_vertex_node_numbers(std::vector<int>* face_index_to_face_vnode_numbers) const;
  const ms::geo::Geometry&      get_geometry(void) const;
  bool                          is_outward_face(const Element& face_element) const;
  Element                       make_face_element(const int face_index) const;
  int                           num_nodes(void) const;
  void                          node_numberss(int* node_numberss) const;
  Element_Type                  type(void) const;
  std::vector<int>              vertex_node_numbers(void) const;
  void                          vertex_node_numbers(int* vertex_node_numbers) const;

private:
  Element_Type                             _type;
  std::vector<ms::geo::Numbered_Node_View> _numbered_node_views;
  ms::geo::Geometry                        _geometry;
};

} // namespace ms::grid