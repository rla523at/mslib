#include "Grid.h"

#include "Data.h"

#include "msexception/Exception.h"
#include "msgeo/Figure.h"

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

    //if (element_type == Element_Type::CELL)
    //{
    //  this->_cell_elements.push_back(element);
    //}
    //else
    //{
    //  this->_boundary_elements.push_back(element);
    //}
  }
}

} // namespace ms::grid