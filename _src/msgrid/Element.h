#pragma once
#include "Data.h"
#include "msgeo/Geometry.h"
#include "msgeo/Node.h"
#include "msmath/Vector.h"
#include <vector>

// data definition
namespace ms::grid
{

// Assuming that the nodes that compose the vertices of the Element are listed first.
struct Numbered_Nodes
{
  std::vector<int>                         numbers;
  std::vector<ms::geo::Node_View> nodes;
};

} // namespace ms::grid

/*





*/

// class declaration
namespace ms::grid
{

class Element
{
public:
  Element(const Element_Type element_type, Numbered_Nodes&& numbered_nodes, ms::geo::Geometry&& geometry)
      : _type(element_type),
        _numbered_nodes(std::move(numbered_nodes)),
        _geometry(std::move(geometry)){};

public:
  void reordering_nodes(const std::vector<int>& new_ordered_node_indexes);

public:
  void                          accumulate_discrete_partition_data(ms::geo::Partition_Data& partition_data, const int partition_order) const;
  void                          accumulate_partition_data(ms::geo::Partition_Data& partition_data, const int partition_order) const;
  int                           dimension(void) const;
  std::vector<int>              find_periodic_matched_node_indexes(const ms::math::Vector_Const_Wrapper& direction_vector, const Element& other) const;
  std::vector<std::vector<int>> face_index_to_face_vertex_node_numbers(void) const;
  void                          face_index_to_face_vertex_node_numbers(std::vector<int>* face_index_to_face_vnode_numbers) const;
  const ms::geo::Geometry&      get_geometry(void) const;
  bool                          is_outward_face(const Element& face_element) const;
  Element                       make_face_element(const int face_index) const;
  int                           num_nodes(void) const;
  void                          node_numbers(int* node_numbers) const;
  Element_Type                  type(void) const;
  std::vector<int>              vertex_node_numbers(void) const;
  void                          vertex_node_numbers(int* vertex_node_numbers) const;

private:
  Element_Type      _type;
  Numbered_Nodes    _numbered_nodes;
  ms::geo::Geometry _geometry;
};

} // namespace ms::grid