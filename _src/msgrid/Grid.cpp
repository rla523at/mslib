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

  for (auto& [element_type, figure, figure_order, node_indexes] : data.element_datas)
  {
    const auto geo_figure       = convert(figure);
    auto       consisting_nodes = this->_nodes.nodes_at_indexes(node_indexes);
    auto       geometry         = ms::geo::Geometry(geo_figure, std::move(consisting_nodes));
    auto       element          = Element(element_type, std::move(node_indexes), std::move(geometry));

    if (element_type == Element_Type::CELL)
    {
      this->_cell_elements.push_back(element);
    }
    else
    {
      this->_boundary_elements.push_back(element);
    }
  }

  const auto num_peridoic_datas         = data.periodic_datas.size();
  const auto num_periodic_element_pairs = num_peridoic_datas / 2;
  this->_periodic_boundary_element_pairs.reserve(num_periodic_element_pairs);

  // sorting
  std::map<std::vector<double>, std::vector<Element>> direction_to_periodic_elements;
  for (auto& [element_data, periodic_direction] : data.periodic_datas)
  {
    auto& [element_type, figure, figure_order, node_indexes] = element_data;

    const auto geo_figure       = convert(figure);
    auto       consisting_nodes = this->_nodes.nodes_at_indexes(node_indexes);
    auto       geometry         = ms::geo::Geometry(geo_figure, std::move(consisting_nodes));
    auto       element          = Element(element_type, std::move(node_indexes), std::move(geometry));

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

  for (auto& [direction, periodic_elements] : direction_to_periodic_elements)
  {
    const auto num_periodic_elements = periodic_elements.size();

    std::unordered_set<int> matched_index_set;
    matched_index_set.reserve(num_periodic_elements);

    for (int i = 0; i < num_periodic_elements; ++i)
    {
      if (matched_index_set.contains(i))
      {
        continue;
      }

      auto& i_element = periodic_elements[i];

      for (int j = i + 1; j < num_periodic_elements; ++j)
      {
        if (matched_index_set.contains(j))
        {
          continue;
        }

        auto& j_element = periodic_elements[j];

        auto periodic_matched_node_indexes = j_element.find_periodic_matched_node_indexes(i_element);

        if (!periodic_matched_node_indexes.empty())
        {
          j_element.rearrange_node_indexes(std::move(periodic_matched_node_indexes));

          matched_periodic_element_pairs.push_back(std::make_pair(std::move(i_element), std::move(j_element)));
          matched_index_set.insert(i);
          matched_index_set.insert(j);
          break;
        }
      }
    }
  }
}

} // namespace ms::grid