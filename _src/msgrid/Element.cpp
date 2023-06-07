#include "Element.h"

#include "msexception/Exception.h"
#include "msmath/Vector.h"
#include <algorithm>
#include <functional>
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

void Element::accumulate_discrete_partition_data(ms::geo::Partition_Data& total_partition_data, const int partition_order) const
{
  REQUIRE(0 <= partition_order, "partition order should not be negative");

  auto&      nodes_in_total_partition_data          = total_partition_data.nodes;
  auto&      connectivities_of_total_connectivities = total_partition_data.connectivities;
  const auto start_node_number                      = nodes_in_total_partition_data.num_nodes();
  const auto num_connectivity                       = connectivities_of_total_connectivities.size();

  auto       partition_data              = this->get_geometry().make_partition_data(partition_order);
  auto&      nodes_in_partition          = partition_data.nodes;
  auto&      connectivities_of_partition = partition_data.connectivities;
  const auto num_new_connectivity        = connectivities_of_partition.size();

  // accumulate nodes in partition
  nodes_in_total_partition_data.add_nodes(std::move(nodes_in_partition));

  // acumulate connectivities of partitions
  connectivities_of_total_connectivities.reserve(num_connectivity + num_new_connectivity);

  std::vector<int> new_connectivity;
  for (int i = 0; i < num_new_connectivity; ++i)
  {
    const auto& ref_connectivity = connectivities_of_partition[i];

    for (const auto node_index : ref_connectivity)
    {
      const auto node_number = start_node_number + node_index;
      new_connectivity.push_back(node_number);
    }
    connectivities_of_total_connectivities.push_back(std::move(new_connectivity));
  }
}

void Element::accumulate_partition_data(ms::geo::Partition_Data& total_partition_data, const int partition_order) const
{
  struct Near
  {
  public:
    bool operator==(const ms::geo::Node_Const_Wrapper other) const
    {
      constexpr auto epsilon = 1.0e-10;

      for (int i = 0; i < node.dimension(); ++i)
      {
        if (epsilon < std::abs(node[i] - other[i])) return false;
      }

      return true;
    }

  public:
    ms::geo::Node_Const_Wrapper node;
  };

  /*





  */

  auto&       nodes_in_total_partition        = total_partition_data.nodes;
  auto&       connectivity_of_total_partition = total_partition_data.connectivities;
  const auto& nodes_in_elem                   = this->_numbered_nodes.nodes;

  auto        partition_data              = this->get_geometry().make_partition_data(partition_order);
  const auto& nodes_in_partition          = partition_data.nodes;
  auto&       connectivities_of_partition = partition_data.connectivities;

  const auto num_nodes_in_partition = nodes_in_partition.num_nodes();

  std::map<int, int> index_to_node_number;

  for (int i = 0; i < num_nodes_in_partition; ++i)
  {
    const auto node_in_partition = nodes_in_partition[i];

    const auto iter = std::find(nodes_in_elem.begin(), nodes_in_elem.end(), Near(node_in_partition));

    if (iter == nodes_in_elem.end())
    {
      nodes_in_total_partition.add_node(node_in_partition);
      const auto node_number = nodes_in_total_partition.num_nodes();

      index_to_node_number.emplace(i, node_number);
    }
    else
    {
      const auto elme_index  = iter - nodes_in_elem.begin();
      const auto node_number = this->_numbered_nodes.numbers[elme_index];

      index_to_node_number.emplace(i, node_number);
    }
  }

  for (auto& connectivity : connectivities_of_partition)
  {
    for (auto& index : connectivity)
    {
      index = index_to_node_number.at(index);
    }
  }

  connectivity_of_total_partition.insert(connectivity_of_total_partition.end(), connectivities_of_partition.begin(), connectivities_of_partition.end());
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

bool Element::is_outward_face(const Element& face_element) const
{
  REQUIRE(this->_type == Element_Type::CELL, "owner cell should be cell type");
  REQUIRE(face_element._type != Element_Type::CELL, "face element should not be cell type");

  const auto& face_geometry = face_element.get_geometry();

  // convention
  if (face_geometry.is_point()) return true;

  const auto bdry_vnode_numbers               = face_element.vertex_node_numbers();
  const auto face_index_to_face_vnode_numbers = this->face_index_to_face_vertex_node_numbers();

  if (face_geometry.is_line())
  {
    for (const auto& face_vnode_numbers : face_index_to_face_vnode_numbers)
    {
      if (!std::is_permutation(face_vnode_numbers.begin(), face_vnode_numbers.end(), bdry_vnode_numbers.begin()))
      {
        continue;
      }

      if (bdry_vnode_numbers == face_vnode_numbers)
      {
        return false;
      }
      else
      {
        return true;
      }
    }
  }
  else
  {
    std::boyer_moore_searcher searcher(bdry_vnode_numbers.begin(), bdry_vnode_numbers.end());

    std::vector<int> temp;
    for (const auto& face_vnode_numbers : face_index_to_face_vnode_numbers)
    {
      if (!std::is_permutation(face_vnode_numbers.begin(), face_vnode_numbers.end(), bdry_vnode_numbers.begin()))
      {
        continue;
      }

      // check circular permutation
      temp.clear();
      temp.insert(temp.end(), face_vnode_numbers.begin(), face_vnode_numbers.end());
      temp.insert(temp.end(), face_vnode_numbers.begin(), face_vnode_numbers.end());

      if (std::search(temp.begin(), temp.end(), searcher) == temp.end())
      {
        return false;
      }
      else
      {
        return true;
      }
    }
  }

  EXCEPTION("This is not my face!");
  return false;
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

Element_Type Element::type(void) const
{
  return this->_type;
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