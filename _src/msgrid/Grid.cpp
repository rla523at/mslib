#include "Grid.h"

#include "Data.h"

#include "msexception/Exception.h"
#include "msgeo/Figure.h"
#include "msmath/Vector.h"
#include <map>
#include <unordered_set>

ms::geo::Coordinates_Type convert(const ms::grid::Coordinate_Type type)
{
  switch (type)
  {
  case ms::grid::Coordinate_Type::NODAL:
    return ms::geo::Coordinates_Type::NODAL;
  case ms::grid::Coordinate_Type::BLOCK:
    return ms::geo::Coordinates_Type::BLOCK;
  default:
    EXCEPTION("unsupported coordinate type");
    return ms::geo::Coordinates_Type::NOT_SUPPROTED;
    break;
  }
}

ms::geo::Figure convert(const ms::grid::Figure type)
{
  switch (type)
  {
  case ms::grid::Figure::POINT:
    return ms::geo::Figure::POINT;
  case ms::grid::Figure::LINE:
    return ms::geo::Figure::LINE;
  case ms::grid::Figure::TRIANGLE:
    return ms::geo::Figure::TRIANGLE;
  case ms::grid::Figure::QUADRILATERAL:
    return ms::geo::Figure::QUADRILATERAL;
  case ms::grid::Figure::TETRAHEDRAL:
    return ms::geo::Figure::TETRAHEDRAL;
  case ms::grid::Figure::HEXAHEDRAL:
    return ms::geo::Figure::HEXAHEDRAL;
  case ms::grid::Figure::PRISM:
    return ms::geo::Figure::PRISM;
  case ms::grid::Figure::PYRAMID:
    return ms::geo::Figure::PYRAMID;
  case ms::grid::Figure::NUM_FIGURES:
    return ms::geo::Figure::NUM_FIGURES;
  default:
    EXCEPTION("not supproted Figrue type");
    return ms::geo::Figure::NOT_IN_LIST;
  }
}

namespace ms::grid
{

Grid::Grid(Data&& data)
    : _nodes(convert(data.nodes_data.type), data.nodes_data.num_nodes, data.nodes_data.dimension, std::move(data.nodes_data.coordinates))
{
}

void Grid::make_cell_and_boundary_elements(std::vector<Element_Data>&& element_datas)
{
  const auto num_element_datas = element_datas.size();
  this->_cell_elements.reserve(num_element_datas);
  this->_boundary_elements.reserve(num_element_datas);

  std::vector<Indexed_Node>                consisting_indexed_nodes;
  std::vector<ms::geo::Node_Const_Wrapper> consisting_nodes;

  for (auto& [element_type, figure, node_indexes] : element_datas)
  {
    consisting_nodes.clear();
    consisting_indexed_nodes.clear();

    const auto num_consisting_nodes = node_indexes.size();
    consisting_nodes.reserve(num_consisting_nodes);
    consisting_indexed_nodes.reserve(num_consisting_nodes);

    for (auto i = 0; i < num_consisting_nodes; ++i)
    {
      const auto index = node_indexes[i];
      const auto node  = this->_nodes[index];
      consisting_nodes.push_back(node);
      consisting_indexed_nodes.push_back({index, node});
    }

    const auto geo_figure = convert(figure);
    auto       geometry   = ms::geo::Geometry(geo_figure, std::move(consisting_nodes));

    auto element = Element(element_type, std::move(consisting_indexed_nodes), std::move(geometry));

    if (element_type == Element_Type::CELL)
    {
      this->_cell_elements.push_back(element);
    }
    else
    {
      this->_boundary_elements.push_back(element);
    }
  }

  this->_cell_elements.shrink_to_fit();
  this->_boundary_elements.shrink_to_fit();
}

void Grid::make_periodic_boundary_elements(std::vector<Peridoic_Data>&& periodic_datas)
{
  const auto num_peridoic_datas         = periodic_datas.size();
  const auto num_periodic_element_pairs = num_peridoic_datas / 2;
  this->_periodic_boundary_element_pairs.reserve(num_periodic_element_pairs);

  // sort by direction
  std::map<std::vector<double>, std::vector<Element>> direction_to_periodic_elements;

  std::vector<Indexed_Node>                consisting_indexed_nodes;
  std::vector<ms::geo::Node_Const_Wrapper> consisting_nodes;

  for (auto& [element_data, periodic_direction] : periodic_datas)
  {
    // make element from element data
    auto& [element_type, figure, node_indexes] = element_data;

    consisting_nodes.clear();
    consisting_indexed_nodes.clear();

    const auto num_consisting_nodes = node_indexes.size();
    consisting_nodes.reserve(num_consisting_nodes);
    consisting_indexed_nodes.reserve(num_consisting_nodes);

    for (auto i = 0; i < num_consisting_nodes; ++i)
    {
      const auto index = node_indexes[i];
      const auto node  = this->_nodes[index];
      consisting_nodes.push_back(node);
      consisting_indexed_nodes.push_back({index, node});
    }

    const auto geo_figure = convert(figure);
    auto       geometry   = ms::geo::Geometry(geo_figure, std::move(consisting_nodes));

    auto element = Element(element_type, std::move(consisting_indexed_nodes), std::move(geometry));

    // sorting
    bool       is_new_direction  = true;
    const auto direction_vector1 = ms::math::Vector_Const_Wrapper(periodic_direction);

    for (auto& [periodic_direction2, elements] : direction_to_periodic_elements)
    {
      const auto direction_vector2 = ms::math::Vector_Const_Wrapper(periodic_direction2);

      if (direction_vector1.is_parallel(direction_vector2))
      {
        is_new_direction = false;
        elements.push_back(std::move(element));
        break;
      }
    }

    if (is_new_direction)
    {
      direction_to_periodic_elements.emplace(periodic_direction, std::vector<Element>());
      direction_to_periodic_elements[periodic_direction].push_back(std::move(element));
    }
  }

  // make pair of periodic elements
  for (auto& [direction, periodic_elements] : direction_to_periodic_elements)
  {
    const auto num_periodic_elements = periodic_elements.size();

    std::unordered_set<int> matched_index_set;
    matched_index_set.reserve(num_periodic_elements);

    for (int i = 0; i < num_periodic_elements; ++i)
    {
      if (matched_index_set.contains(i)) continue;

      auto& i_element = periodic_elements[i];

      for (int j = i + 1; j < num_periodic_elements; ++j)
      {
        if (matched_index_set.contains(j)) continue;

        auto& j_element            = periodic_elements[j];
        auto  matched_node_indexes = j_element.find_periodic_matched_node_indexes(direction, j_element);

        if (!matched_node_indexes.empty())
        {
          // reordering 하지 않으면 quadrature points 순서가 두 pair element에서 달라진다.
          // 즉, i element에서 0번째 quadrature point가 j element에서 0번째 quadrature point와 붙어 있지 않게 된다.
          j_element.rearrange_node_indexes(std::move(matched_node_indexes));

          this->_periodic_boundary_element_pairs.push_back(std::make_pair(std::move(i_element), std::move(j_element)));
          matched_index_set.insert(i);
          matched_index_set.insert(j);
          break;
        }
      }
    }
  }
}

} // namespace ms::grid