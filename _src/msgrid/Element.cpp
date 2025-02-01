#include "msgrid/Element.h"

#include "msexception/Exception.h"
#include "msmath/Vector.h"
#include <algorithm>
#include <functional>
#include <unordered_set>

namespace ms::grid
{

Element::Element(const Element_Type type, std::vector<int>&& node_numbers, ms::geo::Geometry&& geometry)
    : _type(type),
      _node_numbers(std::move(node_numbers)),
      _geometry(std::move(geometry))
{
  // initialize _face_vertex_node_numbers_s
  
  // point element doesn't need to initialize
  if (this->_geometry.is_point()) return;

  const auto& face_vnode_indexes_s = this->_geometry.get_face_vnode_indexes_s();
  const auto  num_faces            = face_vnode_indexes_s.size();

  // index -> node number
  this->_face_vertex_node_numbers_s.resize(num_faces);

  for (int i = 0; i < num_faces; ++i)
  {
    const auto& face_vnode_indexes = face_vnode_indexes_s[i];
    auto&       face_vnode_numbers = this->_face_vertex_node_numbers_s[i];

    const auto num_face_vnode = face_vnode_indexes.size();
    face_vnode_numbers.resize(num_face_vnode);

    for (int j = 0; j < num_face_vnode; ++j)
    {
      const auto index      = face_vnode_indexes[j];
      face_vnode_numbers[j] = this->_node_numbers[index];
    }
  }
}

void Element::reordering_nodes(const std::vector<int>& new_ordered_node_numbers)
{
  const auto begin = new_ordered_node_numbers.begin();
  const auto end   = new_ordered_node_numbers.end();

  std::vector<int> new_orders;
  new_orders.reserve(this->_node_numbers.size());

  for (const auto& node_number : this->_node_numbers)
  {
    const auto iter = std::find(begin, end, node_number);
    REQUIRE(iter != end, "Fail to find node which match with node number");

    const auto new_order = static_cast<int>(iter - begin);
    new_orders.push_back(new_order);
  }

  this->_geometry.reordering_nodes(new_orders);
}

int Element::dimension(void) const
{
  return this->_geometry.dimension();
}

std::vector<int> Element::find_periodic_matched_node_numbers(const ms::math::Vector_View direction_vector, const Element& other) const
{
  // It returns the node indices in "this element" that match the nodes of the "other element".

  REQUIRE(this->_type == Element_Type::PERIODIC && other._type == Element_Type::PERIODIC, "both element should be periodic");

  const auto this_num_nodes  = this->_geometry.num_nodes();
  const auto other_num_nodes = other._geometry.num_nodes();

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
    const auto other_node_number = other._node_numbers[i];
    const auto other_node_view   = other._geometry.node_view(i);

    bool not_find_pair = true;

    for (int j = 0; j < this_num_nodes; ++j)
    {
      const auto this_node_number = this->_node_numbers[j];
      const auto this_node_view   = this->_geometry.node_view(j);

      if (matched_vnode_number_set.contains(this_node_number)) continue;

      ms::geo::A_to_B_vector(node_to_node_vector, other_node_view, this_node_view);

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

std::span<const int> Element::face_vertex_node_numbers(const int face_index) const
{
  REQUIRE(0 <= face_index && face_index < this->_face_vertex_node_numbers_s.size(), "face index is out of range");
  return this->_face_vertex_node_numbers_s[face_index];
}

const ms::geo::Geometry& Element::get_geometry(void) const
{
  return this->_geometry;
}

std::span<const int> Element::node_numbers(void) const
{
  return this->_node_numbers;
}

bool Element::is_outward_face(const Element& face_element) const
{
  REQUIRE(this->_type == Element_Type::CELL, "owner cell should be cell type");
  REQUIRE(face_element._type != Element_Type::CELL, "face element should not be cell type");

  const auto face_index = this->find_face_index(face_element);

  // check direction
  const auto& face_geometry           = face_element.get_geometry();
  const auto& this_face_vnode_numbers = this->_face_vertex_node_numbers_s[face_index];
  const auto  face_vnode_numbers      = face_element.vertex_node_numbers();

  if (face_geometry.is_point())
  {
    // �����ʿ� �ִ� face�� outward face�� �����Ѵ�.
    return face_index == 1;
  }
  else if (face_geometry.is_line())
  {
    // �ݽð� ���� face�� outward face�� �����Ѵ�.
    return std::equal(this_face_vnode_numbers.begin(), this_face_vnode_numbers.end(), face_vnode_numbers.begin());
  }
  else
  {
    // consider circular permutation
    std::vector<int> temp;
    temp.reserve(this_face_vnode_numbers.size() * 2);
    temp.insert(temp.end(), this_face_vnode_numbers.begin(), this_face_vnode_numbers.end());
    temp.insert(temp.end(), this_face_vnode_numbers.begin(), this_face_vnode_numbers.end());

    // �ݽð� ���� face�� outward face�� �����Ѵ�.
    std::boyer_moore_searcher searcher(face_vnode_numbers.begin(), face_vnode_numbers.end());
    return std::search(temp.begin(), temp.end(), searcher) != temp.end();
  }
}

Element Element::make_face_element(const int face_index) const
{
  const auto face_figure       = this->_geometry.face_figure(face_index);
  const auto face_node_indexes = this->_geometry.face_node_indexes(face_index);
  const auto num_face_nodes    = face_node_indexes.size();

  std::vector<ms::geo::Node_View> face_nodes(num_face_nodes);
  std::vector<int>                face_node_numbers(num_face_nodes);

  for (int i = 0; i < num_face_nodes; ++i)
  {
    const auto index     = face_node_indexes[i];
    face_nodes[i]        = this->_geometry.node_view(index);
    face_node_numbers[i] = this->_node_numbers[index];
  }

  auto face_geometry = ms::geo::Geometry(face_figure, std::move(face_nodes));
  return Element(Element_Type::FACE, std::move(face_node_numbers), std::move(face_geometry));
}

std::vector<int> Element::make_sorted_face_vertex_node_numbers(const int face_index) const
{
  const auto vnode_numbers = this->face_vertex_node_numbers(face_index);

  std::vector<int> sorted_vertex_node_numbers(vnode_numbers.size());
  std::copy(vnode_numbers.begin(), vnode_numbers.end(), sorted_vertex_node_numbers.begin());
  std::sort(sorted_vertex_node_numbers.begin(), sorted_vertex_node_numbers.end());

  return sorted_vertex_node_numbers;
}

void Element::make_sorted_face_vertex_node_numbers(const int face_index, std::vector<int>& sorted_vertex_node_numbers) const
{
  const auto vnode_numbers = this->face_vertex_node_numbers(face_index);
  sorted_vertex_node_numbers.resize(vnode_numbers.size());
  std::copy(vnode_numbers.begin(), vnode_numbers.end(), sorted_vertex_node_numbers.begin());
  std::sort(sorted_vertex_node_numbers.begin(), sorted_vertex_node_numbers.end());
}

std::vector<int> Element::make_sorted_vertex_node_numbers(void) const
{
  const auto vnode_numbers = this->vertex_node_numbers();

  std::vector<int> sorted_vertex_node_numbers(vnode_numbers.size());
  std::copy(vnode_numbers.begin(), vnode_numbers.end(), sorted_vertex_node_numbers.begin());
  std::sort(sorted_vertex_node_numbers.begin(), sorted_vertex_node_numbers.end());

  return sorted_vertex_node_numbers;
}

void Element::make_sorted_vertex_node_numbers(std::vector<int>& sorted_vertex_node_numbers) const
{
  const auto vnode_numbers = this->vertex_node_numbers();
  sorted_vertex_node_numbers.resize(vnode_numbers.size());
  std::copy(vnode_numbers.begin(), vnode_numbers.end(), sorted_vertex_node_numbers.begin());
  std::sort(sorted_vertex_node_numbers.begin(), sorted_vertex_node_numbers.end());
}

Element_Type Element::type(void) const
{
  return this->_type;
}

std::span<const int> Element::vertex_node_numbers(void) const
{
  const auto num_vertices = this->_geometry.num_vertices();
  const auto result       = std::span<const int>(this->_node_numbers.begin(), num_vertices); // convention
  return result;
}

int Element::find_face_index(const Element& face_element) const
{
  const auto face_vnode_numbers = face_element.vertex_node_numbers();

  int face_index = -1;
  for (int i = 0; i < this->_face_vertex_node_numbers_s.size(); ++i)
  {
    const auto& this_face_vnode_numbers = _face_vertex_node_numbers_s[i];
    if (!std::is_permutation(this_face_vnode_numbers.begin(), this_face_vnode_numbers.end(), face_vnode_numbers.begin())) continue;

    face_index = i;
    break;
  }

  REQUIRE(face_index != -1, "This is not my face!");
  return face_index;
}

} // namespace ms::grid

//
// std::vector<std::vector<int>> Element::face_vertex_node_numbers_s(void) const
//{
//  const auto& face_vnode_indexes_s = this->_geometry.get_face_vnode_indexes_s();
//  const auto  num_faces            = face_vnode_indexes_s.size();
//
//  // index�� node number�� �ٲٴ� ����
//  std::vector<std::vector<int>> face_vnode_numbers_s(num_faces);
//
//  for (int i = 0; i < num_faces; ++i)
//  {
//    const auto& face_vnode_indexes = face_vnode_indexes_s[i];
//    auto&       face_vnode_numbers = face_vnode_numbers_s[i];
//
//    const auto num_face_vnode = face_vnode_indexes.size();
//    face_vnode_numbers.resize(num_face_vnode);
//
//    for (int j = 0; j < num_face_vnode; ++j)
//    {
//      const auto index      = face_vnode_indexes[j];
//      face_vnode_numbers[j] = this->_node_numbers[index];
//    }
//  }
//
//  return face_vnode_numbers_s;
//}
//
// void Element::face_vertex_node_numbers_s(std::vector<int>* face_vnode_numbers_s) const
//{
//  const auto& face_vnode_indexes_s = this->_geometry.get_face_vnode_indexes_s();
//  const auto  num_faces            = face_vnode_indexes_s.size();
//
//  for (int i = 0; i < num_faces; ++i)
//  {
//    const auto& face_vnode_indexes = face_vnode_indexes_s[i];
//    const auto  num_face_vnode     = face_vnode_indexes.size();
//
//    auto& face_vnode_numbers = face_vnode_numbers_s[i];
//    face_vnode_numbers.resize(num_face_vnode);
//
//    for (int j = 0; j < num_face_vnode; ++j)
//    {
//      const auto index      = face_vnode_indexes[j];
//      face_vnode_numbers[j] = this->_node_numbers[index];
//    }
//  }
//}

// void Element::accumulate_discrete_node_info(ms::geo::Geometry_Consisting_Nodes_Info& discrete_node_info, const int partition_order) const
//{
//   REQUIRE(0 <= partition_order, "partition order should not be negative");
//
//   auto&      nodes             = discrete_node_info.nodes;
//   auto&      connectivities    = discrete_node_info.connectivities;
//   const auto start_node_number = nodes.num_nodes();
//   const auto num_connectivity  = connectivities.size();
//
//   auto       partition_geometry_node_info         = this->get_geometry().make_partitioned_geometry_node_info(partition_order);
//   auto&      nodes_in_partition_geometry          = partition_geometry_node_info.nodes;
//   auto&      connectivities_of_partition_geometry = partition_geometry_node_info.connectivities;
//   const auto num_new_connectivity                 = connectivities_of_partition_geometry.size();
//
//   // accumulate nodes in partition
//   nodes.add_nodes(std::move(nodes_in_partition_geometry));
//
//   // acumulate connectivities of partitions
//   connectivities.reserve(num_connectivity + num_new_connectivity);
//
//   std::vector<int> new_connectivity;
//   for (int i = 0; i < num_new_connectivity; ++i)
//   {
//     const auto& ref_connectivity = connectivities_of_partition_geometry[i];
//
//     for (const auto node_index : ref_connectivity)
//     {
//       const auto node_number = start_node_number + node_index;
//       new_connectivity.push_back(node_number);
//     }
//     connectivities.push_back(std::move(new_connectivity));
//   }
// }
//
// void Element::accumulate_node_info(ms::geo::Geometry_Consisting_Nodes_Info& node_info, const int partition_order) const
//{
//   // geometry���� ���� partitioned geometry node info�� �˸��� connectivity�� �ִ� ����
//   auto        pg_node_info         = this->_geometry.make_partitioned_geometry_node_info(partition_order);
//   const auto& numbered_nodes_in_pg = pg_node_info.numbered_nodes;
//   auto&       connectivities_of_pg = pg_node_info.connectivities;
//
//   int temp_node_number = -1;
//   for (const auto& numbered_node_in_pg : numbered_nodes_in_pg)
//   {
//     const auto& node_in_pg           = numbered_node_in_pg.node;
//     const auto  number_of_node_in_pg = numbered_node_in_pg.number;
//
//     for (const auto& numbered_node_in_elem : this->_numbered_node_views)
//     {
//       constexpr auto epsilon = 1.0e-10;
//
//       const auto& node_in_elem           = numbered_node_in_elem.node_view;
//       const auto  number_of_node_in_elem = numbered_node_in_elem.number;
//
//       if (node_in_pg.distance(node_in_elem) <= epsilon)
//       {
//         for (auto& connectivity : connectivities_of_pg)
//         {
//           connectivity.change_number(number_of_node_in_pg, number_of_node_in_elem);
//         }
//       }
//       else
//       {
//         for (auto& connectivity : connectivities_of_pg)
//         {
//           connectivity.change_number(number_of_node_in_pg, temp_node_number);
//         }
//         temp_node_number--;
//       }
//     }
//   }
//
//   for (int i = 0; i < num_nodes_in_pg; ++i)
//   {
//     const auto node_in_pg = numbered_nodes_in_pg[i];
//     const auto iter       = std::find(nodes_in_elem.begin(), nodes_in_elem.end(), Near(node_in_pg));
//
//     if (iter == nodes_in_elem.end())
//     {
//       nodes.add_node(node_in_pg);
//       const auto node_number = nodes.num_nodes();
//
//       index_to_node_number.emplace(i, node_number);
//     }
//     else
//     {
//       const auto elme_index  = iter - nodes_in_elem.begin();
//       const auto node_number = this->_numbered_nodes.numbers[elme_index];
//
//       index_to_node_number.emplace(i, node_number);
//     }
//   }
//
//   ////
//
//   struct Near
//   {
//   public:
//     bool operator==(const ms::geo::Node_View other) const
//     {
//       constexpr auto epsilon = 1.0e-10;
//       return node.distance(other) <= epsilon;
//     }
//
//   public:
//     ms::geo::Node_View node;
//   };
//
//   /*
//
//
//
//
//
//
//
//
//
//
//   */
//
//   auto&       nodes          = node_info.nodes;
//   auto&       connectivities = node_info.connectivities;
//   const auto& nodes_in_elem  = this->_numbered_nodes.nodes;
//
//   auto        pg_node_info         = this->_geometry.make_partitioned_geometry_node_info(partition_order);
//   const auto& numbered_nodes_in_pg = pg_node_info.numbered_nodes;
//   auto&       connectivities_of_pg = pg_node_info.connectivities;
//
//   const auto num_nodes_in_pg = numbered_nodes_in_pg.num_nodes();
//
//   // nodes in pg�� ������ number�� �ִ� �����̱���.
//   // numbered nodes.
//   // connectivity�� numberd nodes�� �����. �ǹ̰� ��Ȯ������!
//   // nodeinfo�� nodes�� �ƴ϶� numbered nodes�� connectivity�� ����.
//   std::map<int, int> index_to_node_number;
//
//   for (int i = 0; i < num_nodes_in_pg; ++i)
//   {
//     const auto node_in_pg = numbered_nodes_in_pg[i];
//     const auto iter       = std::find(nodes_in_elem.begin(), nodes_in_elem.end(), Near(node_in_pg));
//
//     if (iter == nodes_in_elem.end())
//     {
//       nodes.add_node(node_in_pg);
//       const auto node_number = nodes.num_nodes();
//
//       index_to_node_number.emplace(i, node_number);
//     }
//     else
//     {
//       const auto elme_index  = iter - nodes_in_elem.begin();
//       const auto node_number = this->_numbered_nodes.numbers[elme_index];
//
//       index_to_node_number.emplace(i, node_number);
//     }
//   }
//
//   for (auto& connectivity : connectivities_of_pg)
//   {
//     for (auto& index : connectivity)
//     {
//       index = index_to_node_number.at(index);
//     }
//   }
//
//   connectivities.insert(connectivities.end(), connectivities_of_pg.begin(), connectivities_of_pg.end());
// }


