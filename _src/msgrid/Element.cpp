#include "msgrid/Element.h"

#include "msexception/Exception.h"
#include "msmath/Vector.h"
#include <algorithm>
#include <functional>
#include <unordered_set>

namespace ms::grid
{

void Element::reordering_nodes(const std::vector<int>& new_ordered_node_numbers)
{
  const auto begin = new_ordered_node_numbers.begin();
  const auto end   = new_ordered_node_numbers.end();

  const auto original_numbered_node_views = this->_numbered_node_views;

  const auto num_nodes = this->_numbered_node_views.size();

  for (int i = 0; i < num_nodes; ++i)
  {
    const auto& original_numbered_node_view = original_numbered_node_views[i];
    const auto  number                      = original_numbered_node_view.number;

    const auto iter = std::find(begin, end, number);
    REQUIRE(iter != end, "For the new face_node_numbers to be valid, they must include all the preceding face_node_numbers");

    const auto new_pos = iter - begin;
    nodes[new_pos]     = original_node;
  }
}

void Element::accumulate_discrete_node_info(ms::geo::Geometry_Consisting_Nodes_Info& discrete_node_info, const int partition_order) const
{
  REQUIRE(0 <= partition_order, "partition order should not be negative");

  auto&      nodes             = discrete_node_info.nodes;
  auto&      connectivities    = discrete_node_info.connectivities;
  const auto start_node_number = nodes.num_nodes();
  const auto num_connectivity  = connectivities.size();

  auto       partition_geometry_node_info         = this->get_geometry().make_partitioned_geometry_node_info(partition_order);
  auto&      nodes_in_partition_geometry          = partition_geometry_node_info.nodes;
  auto&      connectivities_of_partition_geometry = partition_geometry_node_info.connectivities;
  const auto num_new_connectivity                 = connectivities_of_partition_geometry.size();

  // accumulate nodes in partition
  nodes.add_nodes(std::move(nodes_in_partition_geometry));

  // acumulate connectivities of partitions
  connectivities.reserve(num_connectivity + num_new_connectivity);

  std::vector<int> new_connectivity;
  for (int i = 0; i < num_new_connectivity; ++i)
  {
    const auto& ref_connectivity = connectivities_of_partition_geometry[i];

    for (const auto node_index : ref_connectivity)
    {
      const auto node_number = start_node_number + node_index;
      new_connectivity.push_back(node_number);
    }
    connectivities.push_back(std::move(new_connectivity));
  }
}

void Element::accumulate_node_info(ms::geo::Geometry_Consisting_Nodes_Info& node_info, const int partition_order) const
{
  // geometry에서 만든 partitioned geometry node info에 알맞은 connectivity를 주는 과정
  auto        pg_node_info         = this->_geometry.make_partitioned_geometry_node_info(partition_order);
  const auto& numbered_nodes_in_pg = pg_node_info.numbered_nodes;
  auto&       connectivities_of_pg = pg_node_info.connectivities;

  int temp_node_number = -1;
  for (const auto& numbered_node_in_pg : numbered_nodes_in_pg)
  {
    const auto& node_in_pg           = numbered_node_in_pg.node;
    const auto  number_of_node_in_pg = numbered_node_in_pg.number;

    for (const auto& numbered_node_in_elem : this->_numbered_node_views)
    {
      constexpr auto epsilon = 1.0e-10;

      const auto& node_in_elem           = numbered_node_in_elem.node_view;
      const auto  number_of_node_in_elem = numbered_node_in_elem.number;

      if (node_in_pg.distance(node_in_elem) <= epsilon)
      {
        for (auto& connectivity : connectivities_of_pg)
        {
          connectivity.change_number(number_of_node_in_pg, number_of_node_in_elem);
        }
      }
      else
      {
        for (auto& connectivity : connectivities_of_pg)
        {
          connectivity.change_number(number_of_node_in_pg, temp_node_number);
        }
        temp_node_number--;
      }
    }
  }

  for (int i = 0; i < num_nodes_in_pg; ++i)
  {
    const auto node_in_pg = numbered_nodes_in_pg[i];
    const auto iter       = std::find(nodes_in_elem.begin(), nodes_in_elem.end(), Near(node_in_pg));

    if (iter == nodes_in_elem.end())
    {
      nodes.add_node(node_in_pg);
      const auto node_number = nodes.num_nodes();

      index_to_node_number.emplace(i, node_number);
    }
    else
    {
      const auto elme_index  = iter - nodes_in_elem.begin();
      const auto node_number = this->_numbered_nodes.numbers[elme_index];

      index_to_node_number.emplace(i, node_number);
    }
  }

  ////

  struct Near
  {
  public:
    bool operator==(const ms::geo::Node_View other) const
    {
      constexpr auto epsilon = 1.0e-10;
      return node.distance(other) <= epsilon;
    }

  public:
    ms::geo::Node_View node;
  };

  /*










  */

  auto&       nodes          = node_info.nodes;
  auto&       connectivities = node_info.connectivities;
  const auto& nodes_in_elem  = this->_numbered_nodes.nodes;

  auto        pg_node_info         = this->_geometry.make_partitioned_geometry_node_info(partition_order);
  const auto& numbered_nodes_in_pg = pg_node_info.numbered_nodes;
  auto&       connectivities_of_pg = pg_node_info.connectivities;

  const auto num_nodes_in_pg = numbered_nodes_in_pg.num_nodes();

  // nodes in pg에 적절한 number를 주는 과정이구나.
  // numbered nodes.
  // connectivity를 numberd nodes와 엮어보자. 의미가 명확해진다!
  // nodeinfo에 nodes가 아니라 numbered nodes와 connectivity로 보자.
  std::map<int, int> index_to_node_number;

  for (int i = 0; i < num_nodes_in_pg; ++i)
  {
    const auto node_in_pg = numbered_nodes_in_pg[i];
    const auto iter       = std::find(nodes_in_elem.begin(), nodes_in_elem.end(), Near(node_in_pg));

    if (iter == nodes_in_elem.end())
    {
      nodes.add_node(node_in_pg);
      const auto node_number = nodes.num_nodes();

      index_to_node_number.emplace(i, node_number);
    }
    else
    {
      const auto elme_index  = iter - nodes_in_elem.begin();
      const auto node_number = this->_numbered_nodes.numbers[elme_index];

      index_to_node_number.emplace(i, node_number);
    }
  }

  for (auto& connectivity : connectivities_of_pg)
  {
    for (auto& index : connectivity)
    {
      index = index_to_node_number.at(index);
    }
  }

  connectivities.insert(connectivities.end(), connectivities_of_pg.begin(), connectivities_of_pg.end());
}

int Element::dimension(void) const
{
  return this->_numbered_nodes.nodes.front().dimension();
}

std::vector<int> Element::find_periodic_matched_node_numbers(const ms::math::Vector_View direction_vector, const Element& other) const
{
  // It returns the node indices in "this element" that match the nodes of the "other element".

  REQUIRE(this->_type == Element_Type::PERIODIC && other._type == Element_Type::PERIODIC, "both element should be periodic");

  const auto this_num_nodes  = this->num_nodes();
  const auto other_num_nodes = other.num_nodes();

  if (this_num_nodes != other_num_nodes)
  {
    return {};
  }

  std::unordered_set<int> matched_vnode_number_set;
  matched_vnode_number_set.reserve(other_num_nodes);

  std::vector<int> matched_node_numbers(this_num_nodes);

  const auto       dimension = this->dimension();
  ms::math::Vector node_to_node_vector(dimension);

  for (int i = 0; i < other_num_nodes; ++i)
  {
    const auto& other_numbered_node_view = other._numbered_node_views[i];
    const auto  other_node_number        = other_numbered_node_view.number;
    const auto& other_node_view          = other_numbered_node_view.node_view;

    bool not_find_pair = true;

    for (int j = 0; j < this_num_nodes; ++j)
    {
      const auto& this_numbered_node_view = this->_numbered_node_views[j];
      const auto  this_node_number        = this_numbered_node_view.number;
      const auto& this_node_view          = this_numbered_node_view.node_view;

      if (matched_vnode_number_set.contains(this_node_number))
      {
        continue;
      }

      this_node_view.other_to_this_vector(other_node_view, node_to_node_vector);

      if (direction_vector.is_parallel(node_to_node_vector))
      {
        matched_node_numbers[i] = this_node_number;
        matched_vnode_number_set.insert(this_node_number);
        not_find_pair = false;
        break;
      }
    }

    if (not_find_pair)
    {
      return {};
    }
  }

  return matched_node_numbers;
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

void Element::node_numberss(int* node_numberss) const
{
  const auto  num_nodes = this->num_nodes();
  const auto& numbers   = this->_numbered_nodes.numbers;

  for (int i = 0; i < num_nodes; ++i)
  {
    node_numberss[i] = numbers[i];
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