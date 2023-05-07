#include "Element.h"

#include "msexception/Exception.h"
#include <unordered_set>

namespace ms::grid
{

std::vector<int> Element::find_periodic_matched_node_indexes(const Element& other) const
{
  REQUIRE(this->_type == Element_Type::PERIODIC && other._type == Element_Type::PERIODIC, "both element should be periodic");

  const auto this_num_node  = this->_consisting_node_indexes.size();
  const auto other_num_node = other._consisting_node_indexes.size();

  if (this->_consisting_node_indexes.size() != other._consisting_node_indexes.size())
  {
    return {};
  }

  std::unordered_set<int> matched_vnode_index_set;
  matched_vnode_index_set.reserve(other_num_node);

  std::vector<int> matched_periodic_node_indexes(this_num_node);

  const auto& this_nodes  = this->points_;
  const auto& other_nodes = other.points_;

  for (int i = 0; i < other_num_node; ++i)
  {
    const auto& other_node       = other_nodes[i];
    const auto  other_node_index = other.point_indexes_[i];

    for (int j = 0; j < this_num_node; ++j)
    {
      const auto& this_node       = this_nodes[j];
      const auto  this_node_index = this->point_indexes_[j];

      if (matched_vnode_index_set.contains(other_node_index))
      {
        continue;
      }

      if (this_node.is_axis_translation(other_node))
      {
        matched_periodic_node_indexes[i] = this_node_index;
        matched_vnode_index_set.insert(this_node_index);
        break;
      }
    }

    // when can not find pair
    if (i + 1 != matched_vnode_index_set.size())
      return {};
  }

  return matched_periodic_node_indexes;
}

bool Element::can_be_periodic_pair(const Element& other) const
{
  if (this->_consisting_node_indexes.size() != other._consisting_node_indexes.size())
  {
    return false;
  }

  if (this->_geometry.is_on_same_axis_plane(other._geometry))
  {
  }
  return false;

  return true;
}

} // namespace ms::grid