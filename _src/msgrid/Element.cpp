#include "Element.h"

#include "msexception/Exception.h"
#include "msmath/Vector.h"
#include <unordered_set>

namespace ms::grid
{

std::vector<int> Element::find_periodic_matched_node_indexes(const ms::math::Vector_Const_Wrapper& direction_vector, const Element& other) const
{
  REQUIRE(this->_type == Element_Type::PERIODIC && other._type == Element_Type::PERIODIC, "both element should be periodic");

  const auto this_num_node  = this->_consisting_indexed_nodes.size();
  const auto other_num_node = other._consisting_indexed_nodes.size();

  if (this->_consisting_indexed_nodes.size() != other._consisting_indexed_nodes.size())
  {
    return {};
  }

  std::unordered_set<int> matched_vnode_index_set;
  matched_vnode_index_set.reserve(other_num_node);

  std::vector<int> matched_periodic_node_indexes(this_num_node);

  const auto& this_indexed_nodes  = this->_consisting_indexed_nodes;
  const auto& other_indexed_nodes = other._consisting_indexed_nodes;

  const auto       dimension = this->dimension();
  ms::math::Vector node_to_node_vector(dimension);

  for (int i = 0; i < other_num_node; ++i)
  {
    const auto& [other_node_index, other_node] = other_indexed_nodes[i];

    bool not_find_pair = true;

    for (int j = 0; j < this_num_node; ++j)
    {
      const auto& [this_node_index, this_node] = this_indexed_nodes[i];

      if (matched_vnode_index_set.contains(other_node_index))
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

int Element::dimension(void) const
{
  return this->_consisting_indexed_nodes.front().node.dimension();
}

} // namespace ms::grid