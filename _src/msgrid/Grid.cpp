#include "Grid.h"

#include "Data.h"

#include "msexception/Exception.h"
#include "msgeo/Figure.h"
#include "msmath/Vector.h"
#include <algorithm>
#include <map>
#include <set>
#include <unordered_set>

// anonymous function definition
namespace
{
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
    return ms::geo::Figure::NOT_FIGURE;
  }
}
} // namespace

/*




*/

// class definition
namespace ms::grid
{

Grid::Grid(Grid_Data&& data)
    : _grid_nodes(convert(data.nodes_data.type), data.nodes_data.num_nodes, data.nodes_data.dimension, std::move(data.nodes_data.coordinates))
{
  this->make_cell_and_boundary_elements(std::move(data.element_datas));
  this->make_periodic_boundary_elements(std::move(data.periodic_datas));
  this->make_inter_cell_face_elements();
}

void Grid::cell_center(double* cell_center, const int cell_number) const
{
  const auto& geometry = this->_cell_elements[cell_number].get_geometry();
  geometry.center(cell_center);
}

void Grid::cell_projected_volume(double* cell_projected_volume, const int cell_number) const
{
  const auto& geometry = this->_cell_elements[cell_number].get_geometry();
  geometry.cal_projected_volumes(cell_projected_volume);
}

double Grid::cell_volume(const int cell_number) const
{
  const auto& geometry = this->_cell_elements[cell_number].get_geometry();
  return geometry.cal_volume();
}

int Grid::dimension(void) const
{
  return this->_cell_elements.front().dimension();
}

bool Grid::is_line_cell(const int cell_number) const
{
  return this->_cell_elements[cell_number].get_geometry().is_line();
}

int Grid::num_cells(void) const
{
  return static_cast<int>(this->_cell_elements.size());
}

void Grid::make_cell_and_boundary_elements(std::vector<Grid_Element_Data>&& element_datas)
{
  const auto num_element_datas = element_datas.size();
  this->_cell_elements.reserve(num_element_datas);
  this->_boundary_elements.reserve(num_element_datas);

  Numbered_Nodes element_numbered_nodes;
  auto& [numbers, element_nodes] = element_numbered_nodes;

  for (auto& [element_type, figure, node_numbers] : element_datas)
  {
    numbers = std::move(node_numbers);
    this->_grid_nodes.nodes_at_indexes(element_nodes, numbers);

    const auto geo_figure = convert(figure);
    auto       geometry   = ms::geo::Geometry(geo_figure, element_nodes);

    auto element = Element(element_type, std::move(element_numbered_nodes), std::move(geometry));

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

void Grid::make_periodic_boundary_elements(std::vector<Grid_Peridoic_Data>&& periodic_datas)
{
  const auto num_peridoic_datas         = periodic_datas.size();
  const auto num_periodic_element_pairs = num_peridoic_datas / 2;
  this->_periodic_boundary_element_pairs.reserve(num_periodic_element_pairs);

  // sort by direction
  std::map<std::vector<double>, std::vector<Element>> direction_to_periodic_elements;

  Numbered_Nodes consisting_indexed_nodes;
  auto& [numbers, element_consisting_nodes] = consisting_indexed_nodes;

  for (auto& [element_data, periodic_direction] : periodic_datas)
  {
    // make element from element data
    auto& [element_type, figure, node_numbers] = element_data;

    numbers = std::move(node_numbers);
    this->_grid_nodes.nodes_at_indexes(element_consisting_nodes, numbers);

    const auto geo_figure = convert(figure);
    auto       geometry   = ms::geo::Geometry(geo_figure, element_consisting_nodes);

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
        auto  matched_node_indexes = j_element.find_periodic_matched_node_indexes(direction, i_element);

        if (matched_node_indexes.empty()) continue;

        // Reordering is necessary for the quadrature points of the two elements to match
        j_element.reordering_nodes(matched_node_indexes);

        this->_periodic_boundary_element_pairs.push_back(std::make_pair(std::move(i_element), std::move(j_element)));
        matched_index_set.insert(i);
        matched_index_set.insert(j);
        break;
      }
    }
  }
}

void Grid::make_inter_cell_face_elements(void)
{
  // check constructed face vnode numbers
  std::set<std::vector<int>> constructed_face_vnode_numbers_set;

  std::vector<int> vnode_numbers;
  for (const auto& boundray_element : this->_boundary_elements)
  {
    const auto& geometry = boundray_element.get_geometry();

    if (geometry.is_point())
    {
      const auto num_nodes = boundray_element.num_nodes();

      vnode_numbers.resize(num_nodes);
      boundray_element.node_numbers(vnode_numbers.data());
    }
    else
    {
      const auto num_vertices = geometry.num_vertices();

      vnode_numbers.resize(num_vertices);
      boundray_element.vertex_node_numbers(vnode_numbers.data());
      std::sort(vnode_numbers.begin(), vnode_numbers.end()); // to ignore index order
    }

    constructed_face_vnode_numbers_set.insert(std::move(vnode_numbers));
  }

  for (const auto& [pbdry_elem1, pbdry_elem2] : this->_periodic_boundary_element_pairs)
  {
    const auto& geometry = pbdry_elem1.get_geometry();

    if (geometry.is_point())
    {
      const auto num_nodes = pbdry_elem1.num_nodes();

      vnode_numbers.resize(num_nodes);
      pbdry_elem1.node_numbers(vnode_numbers.data());
      constructed_face_vnode_numbers_set.insert(std::move(vnode_numbers));

      vnode_numbers.resize(num_nodes);
      pbdry_elem2.node_numbers(vnode_numbers.data());
      constructed_face_vnode_numbers_set.insert(std::move(vnode_numbers));
    }
    else
    {
      const auto num_vertices = geometry.num_vertices();

      vnode_numbers.resize(num_vertices);
      pbdry_elem1.vertex_node_numbers(vnode_numbers.data());
      std::sort(vnode_numbers.begin(), vnode_numbers.end()); // to ignore index order
      constructed_face_vnode_numbers_set.insert(std::move(vnode_numbers));

      vnode_numbers.resize(num_vertices);
      pbdry_elem2.vertex_node_numbers(vnode_numbers.data());
      std::sort(vnode_numbers.begin(), vnode_numbers.end()); // to ignore index order
      constructed_face_vnode_numbers_set.insert(std::move(vnode_numbers));
    }
  }

  // construct intercell face

  std::vector<std::vector<int>> face_index_to_face_vnode_numbers;
  for (const auto& cell_element : this->_cell_elements)
  {
    const auto& geometry = cell_element.get_geometry();

    const auto num_faces = geometry.num_faces();
    face_index_to_face_vnode_numbers.resize(num_faces);
    cell_element.face_index_to_face_vertex_node_numbers(face_index_to_face_vnode_numbers.data());

    for (int i = 0; i < num_faces; ++i)
    {
      auto& face_vnode_numbers = face_index_to_face_vnode_numbers[i];
      std::sort(face_vnode_numbers.begin(), face_vnode_numbers.end()); // to ignore index order

      if (constructed_face_vnode_numbers_set.contains(face_vnode_numbers))
      {
        continue;
      }
      else
      {
        this->_inter_cell_face_elements.push_back(cell_element.make_face_element(i));
        constructed_face_vnode_numbers_set.insert(std::move(face_vnode_numbers));
      }
    }
  }
}

} // namespace ms::grid