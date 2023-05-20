#include "Element.h"

#include "msexception/Exception.h"
#include "msmath/Vector.h"
#include <unordered_set>

namespace ms::grid
{

void Element::reordering_nodes(const std::vector<int>& new_ordered_node_indexes)
{
  auto& [numbers, nodes]    = this->_numbered_nodes;
  const auto begin          = new_ordered_node_indexes.begin();
  const auto end            = new_ordered_node_indexes.end();
  const auto original_nodes = nodes;

  const auto num_nodes = nodes.size();

  for (int i = 0; i < num_nodes; ++i)
  {
    const auto index         = numbers[i];
    const auto original_node = original_nodes[i];

    const auto iter = std::find(begin, end, index);
    REQUIRE(iter != end, "For the new face_node_numbers to be valid, they must include all the preceding face_node_numbers");

    const auto new_pos = iter - begin;
    nodes[new_pos]     = original_node;
  }
}

int Element::dimension(void) const
{
  return this->_numbered_nodes.nodes.front().dimension();
}

std::vector<int> Element::find_periodic_matched_node_indexes(const ms::math::Vector_Const_Wrapper& direction_vector, const Element& other) const
{
  // It returns the node indices in "this element" that match the nodes of the "other element".

  REQUIRE(this->_type == Element_Type::PERIODIC && other._type == Element_Type::PERIODIC, "both element should be periodic");

  const auto this_num_nodes  = this->num_nodes();
  const auto other_num_nodes = other.num_nodes();

  if (this_num_nodes != other_num_nodes)
  {
    return {};
  }

  std::unordered_set<int> matched_vnode_index_set;
  matched_vnode_index_set.reserve(other_num_nodes);

  std::vector<int> matched_periodic_node_indexes(this_num_nodes);

  const auto& [this_indexes, this_nodes]   = this->_numbered_nodes;
  const auto& [other_indexes, other_nodes] = other._numbered_nodes;

  const auto       dimension = this->dimension();
  ms::math::Vector node_to_node_vector(dimension);

  for (int i = 0; i < other_num_nodes; ++i)
  {
    const auto  other_node_index = other_indexes[i];
    const auto& other_node       = other_nodes[i];

    bool not_find_pair = true;

    for (int j = 0; j < this_num_nodes; ++j)
    {
      const auto  this_node_index = this_indexes[j];
      const auto& this_node       = this_nodes[j];

      if (matched_vnode_index_set.contains(this_node_index))
      {
        continue;
      }

      this_node.other_to_this_vector(other_node, node_to_node_vector);

      if (direction_vector.is_parallel(node_to_node_vector))
      {
        matched_periodic_node_indexes[i] = this_node_index;
        matched_vnode_index_set.insert(this_node_index);
        not_find_pair = false;
        break;
      }
    }

    if (not_find_pair)
    {
      return {};
    }
  }

  return matched_periodic_node_indexes;
}

std::vector<std::vector<int>> Element::face_index_to_face_vertex_node_numbers(void) const
{
  const auto& numbers = this->_numbered_nodes.numbers;

  const auto face_index_to_face_vnode_indexes = this->_geometry.face_index_to_face_vnode_indexes();
  const auto num_faces                        = face_index_to_face_vnode_indexes.size();

  std::vector<std::vector<int>> face_index_to_face_vnode_numbers(num_faces);

  for (int i = 0; i < num_faces; ++i)
  {
    const auto& face_vnode_indexes = face_index_to_face_vnode_indexes[i];
    auto&       face_vnode_numbers = face_index_to_face_vnode_numbers[i];

    const auto num_face_vnode = face_vnode_indexes.size();
    face_vnode_numbers.resize(num_face_vnode);

    for (int j = 0; j < num_face_vnode; ++j)
    {
      const auto index      = face_vnode_indexes[j];
      face_vnode_numbers[j] = numbers[index];
    }
  }

  return face_index_to_face_vnode_numbers;
}

void Element::face_index_to_face_vertex_node_numbers(std::vector<int>* face_index_to_face_vnode_numbers) const
{
  const auto face_index_to_face_vnode_indexes = this->_geometry.face_index_to_face_vnode_indexes();
  const auto num_faces                        = face_index_to_face_vnode_indexes.size();

  const auto& numbers = this->_numbered_nodes.numbers;

  for (int i = 0; i < num_faces; ++i)
  {
    const auto& face_vnode_indexes = face_index_to_face_vnode_indexes[i];
    auto&       face_vnode_numbers = face_index_to_face_vnode_numbers[i];

    const auto num_face_vnode = face_vnode_indexes.size();
    face_vnode_numbers.resize(num_face_vnode);

    for (int j = 0; j < num_face_vnode; ++j)
    {
      const auto index      = face_vnode_indexes[j];
      face_vnode_numbers[j] = numbers[index];
    }
  }
}

const ms::geo::Geometry& Element::get_geometry(void) const
{
  return this->_geometry;
}

Element Element::make_face_element(const int face_index) const
{
  const auto face_figure       = this->_geometry.face_figure(face_index);
  const auto face_node_indexes = this->_geometry.face_node_indexes(face_index);
  const auto num_face_nodes    = face_node_indexes.size();

  Numbered_Nodes face_numbered_nodes;
  auto& [face_node_numbers, face_nodes] = face_numbered_nodes;
  face_nodes.reserve(num_face_nodes);
  face_node_numbers.reserve(num_face_nodes);

  const auto& [numbers, nodes] = this->_numbered_nodes;
  for (const auto index : face_node_indexes)
  {
    face_node_numbers.push_back(numbers[index]);
    face_nodes.push_back(nodes[index]);
  }

  auto face_geometry = ms::geo::Geometry(face_figure, face_nodes);

  return Element(Element_Type::FACE, std::move(face_numbered_nodes), std::move(face_geometry));
}

int Element::num_nodes(void) const
{
  return static_cast<int>(this->_numbered_nodes.nodes.size());
}

void Element::node_numbers(int* node_numbers) const
{
  const auto  num_nodes = this->num_nodes();
  const auto& numbers   = this->_numbered_nodes.numbers;

  for (int i = 0; i < num_nodes; ++i)
  {
    node_numbers[i] = numbers[i];
  }
}

std::vector<int> Element::vertex_node_numbers(void) const
{
  const auto& numbers      = this->_numbered_nodes.numbers;
  const auto  num_vertices = this->_geometry.num_vertices();

  // Assuming that the nodes that compose the vertices of the Element are listed first.
  return {numbers.begin(), numbers.begin() + num_vertices};
}

void Element::vertex_node_numbers(int* vertex_node_numbers) const
{
  const auto& numbers      = this->_numbered_nodes.numbers;
  const auto  num_vertices = this->_geometry.num_vertices();

  for (int i = 0; i < num_vertices; ++i)
  {
    // Assuming that the nodes that compose the vertices of the Element are listed first.
    vertex_node_numbers[i] = numbers[i];
  }
}

} // namespace ms::grid