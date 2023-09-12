#include "Grid.h"

#include "Data.h"

#include "msexception/Exception.h"
#include "msgeo/Figure.h"
#include "msmath/Vector.h"
#include <algorithm>
#include <functional>
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

template <typename Container1, typename Container2, std::enable_if_t<std::is_same_v<typename Container1::value_type, typename Container2::value_type>, bool> = true>
std::vector<typename Container1::value_type> set_intersection(const Container1& container1, const Container2& container2)
{
  std::vector<typename Container1::value_type> intersection;
  std::set_intersection(container1.begin(), container1.end(), container2.begin(), container2.end(), std::back_inserter(intersection));

  return intersection;
}

template <typename Container>
std::vector<typename Container::value_type> set_intersection(const std::vector<Container>& containers)
{
  const auto num_container = containers.size();
  REQUIRE(2 <= num_container, "number of container should be greater than 2");

  const auto& set0         = containers.at(0);
  const auto& set1         = containers.at(1);
  auto        intersection = set_intersection(set0, set1);

  if (2 < num_container)
  {
    for (int i = 2; i < num_container; ++i)
    {
      const auto& set_i = containers.at(i);
      auto        temp  = set_intersection(intersection, set_i);

      std::swap(intersection, temp);
    }
  }

  return intersection;
}

template <typename Container>
std::vector<typename Container::value_type> set_intersection(const std::vector<Container*>& containers)
{
  const auto num_container = containers.size();
  REQUIRE(2 <= num_container, "number of container should be greater than 2");

  const auto set0         = containers.at(0);
  const auto set1         = containers.at(1);
  auto       intersection = set_intersection(*set0, *set1);

  if (2 < num_container)
  {
    for (int i = 2; i < num_container; ++i)
    {
      const auto set_i = containers.at(i);
      auto       temp  = set_intersection(intersection, *set_i);

      std::swap(intersection, temp);
    }
  }

  return intersection;
}

template <typename Container1, typename Container2, std::enable_if_t<std::is_same_v<typename Container1::value_type, typename Container2::value_type>, bool> = true>
std::vector<typename Container1::value_type> set_difference(const Container1& container1, const Container2& container2)
{
  std::vector<typename Container1::value_type> difference;
  std::set_difference(container1.begin(), container1.end(), container2.begin(), container2.end(), std::back_inserter(difference));

  return difference;
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

  // precalculate _vnode_number_to_share_cell_number_set_ignore_pbdry
  const auto num_cell = this->_cell_elements.size();

  for (int i = 0; i < num_cell; ++i)
  {
    const auto vnode_numbers = this->_cell_elements[i].vertex_node_numbers();
    for (const auto vnode_number : vnode_numbers)
    {
      if (!this->_vnode_number_to_share_cell_number_set_ignore_pbdry.contains(vnode_number))
      {
        this->_vnode_number_to_share_cell_number_set_ignore_pbdry.emplace(vnode_number, std::set<int>());
      }

      this->_vnode_number_to_share_cell_number_set_ignore_pbdry.at(vnode_number).insert(i);
    }
  }

  // precalculate _vnode_number_to_share_cell_number_set_consider_pbdry
  this->_vnode_number_to_share_cell_number_set_consider_pbdry = this->_vnode_number_to_share_cell_number_set_ignore_pbdry;

  const auto pbdry_vnode_number_to_matched_vnode_number_set = this->peridoic_boundary_vertex_node_number_to_matched_vertex_node_number_set();

  for (const auto& [pbdry_vnode_number, matched_vnode_number_set] : pbdry_vnode_number_to_matched_vnode_number_set)
  {
    for (const auto matched_vnode_number : matched_vnode_number_set)
    {
      auto&       i_set = this->_vnode_number_to_share_cell_number_set_consider_pbdry.at(pbdry_vnode_number);
      const auto& j_set = this->_vnode_number_to_share_cell_number_set_consider_pbdry.at(matched_vnode_number);

      const auto difference = set_difference(j_set, i_set);

      i_set.insert(difference.begin(), difference.end());
    }
  }
}

int Grid::boundary_owner_cell_number(const int boundary_number) const
{
  const auto vnode_numbers = this->_boundary_elements[boundary_number].vertex_node_numbers();
  const auto cell_numbers  = this->find_cell_numbers_have_these_vertex_nodes_ignore_pbdry(vnode_numbers);
  REQUIRE(cell_numbers.size() == 1, "boundary should have unique owner cell");

  return cell_numbers.front();
}

void Grid::boundary_outward_unit_normal_at_center(double* normal, const int boundary_number, const int owner_cell_number) const
{
  const auto& bdry_element  = this->_boundary_elements[boundary_number];
  const auto& oc_element    = this->_cell_elements[owner_cell_number];
  const auto& bdry_geometry = bdry_element.get_geometry();

  const auto bdry_center = bdry_geometry.center();
  bdry_geometry.cal_normal(normal, bdry_center);

  ms::math::Vector_Wrap normal_v(normal, this->dimension());
  normal_v.normalize();

  if (!oc_element.is_outward_face(bdry_element))
  {
    normal_v *= -1.0;
  }
}

Element_Type Grid::boundary_type(const int boundary_number) const
{
  return this->_boundary_elements[boundary_number].type();
}

double Grid::boundary_volume(const int boundary_number) const
{
  const auto& bdry_geo = this->_boundary_elements[boundary_number].get_geometry();
  return bdry_geo.cal_volume();
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

std::pair<int, int> Grid::inter_cell_face_owner_neighbor_cell_number_pair(const int inter_cell_face_number) const
{
  const auto num_inter_cell_face = static_cast<int>(this->_inter_cell_face_elements.size());

  if (inter_cell_face_number < num_inter_cell_face)
  {
    const auto& element       = this->_inter_cell_face_elements.at(inter_cell_face_number);
    const auto  vnode_numbers = element.vertex_node_numbers();
    const auto  cell_numbers  = this->find_cell_numbers_have_these_vertex_nodes_ignore_pbdry(vnode_numbers);
    REQUIRE(cell_numbers.size() == 2, "inter cell face should have an unique owner neighbor cell pair");

    // set first number as owner cell number
    const auto oc_number = cell_numbers[0];
    const auto nc_number = cell_numbers[1];
    return {oc_number, nc_number};
  }
  else
  {
    const auto pbdry_pair_number = inter_cell_face_number - num_inter_cell_face;
    // set first element as owner cell side element
    const auto& [oc_side_element, nc_side_element] = this->_periodic_boundary_element_pairs.at(pbdry_pair_number);

    const auto oc_numbers = this->find_cell_numbers_have_these_vertex_nodes_ignore_pbdry(oc_side_element.vertex_node_numbers());
    const auto nc_numbers = this->find_cell_numbers_have_these_vertex_nodes_ignore_pbdry(nc_side_element.vertex_node_numbers());
    REQUIRE(oc_numbers.size() == 1, "periodic boundary should have unique owner cell");
    REQUIRE(nc_numbers.size() == 1, "periodic boundary should have unique neighbor cell");

    const auto oc_number = oc_numbers.front();
    const auto nc_number = nc_numbers.front();
    return {oc_number, nc_number};
  }
}

void Grid::inter_cell_face_outward_unit_normal_at_center(double* normal, const int inter_cell_face_number, const int owner_cell_number) const
{
  const auto num_inter_cell_face = this->_inter_cell_face_elements.size();

  ms::math::Vector_Wrap normal_v(normal, this->dimension());

  if (inter_cell_face_number < num_inter_cell_face)
  {
    const auto& infc_element  = this->_inter_cell_face_elements[inter_cell_face_number];
    const auto& infc_geometry = infc_element.get_geometry();
    const auto& oc_element    = this->_cell_elements[owner_cell_number];

    const auto center = infc_geometry.center();
    infc_geometry.cal_normal(normal, center);
    normal_v.normalize();

    if (!oc_element.is_outward_face(infc_element))
    {
      normal_v *= -1.0;
    }
  }
  else
  {
    const auto pbdry_pair_number                   = inter_cell_face_number - num_inter_cell_face;
    const auto& [oc_side_element, nc_side_element] = this->_periodic_boundary_element_pairs.at(pbdry_pair_number);
    const auto& oc_side_geometry                   = oc_side_element.get_geometry();
    const auto& oc_element                         = this->_cell_elements[owner_cell_number];

    const auto center = oc_side_geometry.center();
    oc_side_geometry.cal_normal(normal, center);
    normal_v.normalize();

    if (!oc_element.is_outward_face(oc_side_element))
    {
      normal_v *= -1.0;
    }
  }
}

double Grid::inter_cell_face_volume(const int inter_cell_face_number) const
{
  const auto num_inter_cell_face = this->_inter_cell_face_elements.size();

  if (inter_cell_face_number < num_inter_cell_face)
  {
    const auto& element  = this->_inter_cell_face_elements.at(inter_cell_face_number);
    const auto& geometry = element.get_geometry();

    return geometry.cal_volume();
  }
  else
  {
    const auto pbdry_pair_number                   = inter_cell_face_number - num_inter_cell_face;
    const auto& [oc_side_element, nc_side_element] = this->_periodic_boundary_element_pairs.at(pbdry_pair_number);
    const auto& ocs_geometry                       = oc_side_element.get_geometry();

    return ocs_geometry.cal_volume();
  }
}

int Grid::num_boundaries(void) const
{
  return static_cast<int>(this->_boundary_elements.size());
}

int Grid::num_inter_cell_faces(void) const
{
  return static_cast<int>(this->_inter_cell_face_elements.size() + this->_periodic_boundary_element_pairs.size());
}

int Grid::num_cells(void) const
{
  return static_cast<int>(this->_cell_elements.size());
}

ms::geo::Partition_Data Grid::make_discrete_partition_data(const int partition_order) const
{
  ms::geo::Partition_Data discrete_partition_data;

  for (const auto& cell_elem : _cell_elements)
  {
    cell_elem.accumulate_discrete_partition_data(discrete_partition_data, partition_order);
  }

  return discrete_partition_data;
}

ms::geo::Partition_Data Grid::make_partition_data(const int partition_order) const
{
  ms::geo::Partition_Data partition_data;
  partition_data.nodes = this->_grid_nodes;

  for (const auto& cell_elem : this->_cell_elements)
  {
    cell_elem.accumulate_partition_data(partition_data, partition_order);
  }

  return partition_data;
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

  Numbered_Nodes consisting_numbered_nodes;
  auto& [numbers, element_consisting_nodes] = consisting_numbered_nodes;

  for (auto& [element_data, periodic_direction] : periodic_datas)
  {
    // make element from element data
    auto& [element_type, figure, node_numbers] = element_data;

    numbers = std::move(node_numbers);
    this->_grid_nodes.nodes_at_indexes(element_consisting_nodes, numbers);

    const auto geo_figure = convert(figure);
    auto       geometry   = ms::geo::Geometry(geo_figure, element_consisting_nodes);

    auto element = Element(element_type, std::move(consisting_numbered_nodes), std::move(geometry));

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

std::vector<int> Grid::find_cell_numbers_have_these_vertex_nodes_ignore_pbdry(const std::vector<int>& vnode_numbers) const
{
  const auto num_vnode = vnode_numbers.size();

  std::vector<const std::set<int>*> share_cell_index_set_ptrs;
  share_cell_index_set_ptrs.reserve(num_vnode);

  for (int i = 0; i < num_vnode; ++i)
  {
    share_cell_index_set_ptrs.push_back(&this->_vnode_number_to_share_cell_number_set_ignore_pbdry.at(vnode_numbers[i]));
  }

  return set_intersection(share_cell_index_set_ptrs);
}

std::unordered_map<int, std::set<int>> Grid::peridoic_boundary_vertex_node_number_to_matched_vertex_node_number_set(void) const
{
  std::unordered_map<int, std::set<int>> pbdry_vnode_number_to_matched_vnode_number_set;

  for (const auto& [oc_side_element, nc_side_element] : this->_periodic_boundary_element_pairs)
  {
    const auto oc_side_vnode_numbers = oc_side_element.vertex_node_numbers();
    const auto nc_side_vnode_numbers = nc_side_element.vertex_node_numbers();

    const auto num_vnode = oc_side_vnode_numbers.size();
    for (int i = 0; i < num_vnode; ++i)
    {
      const auto i_vnode_index = oc_side_vnode_numbers[i];
      const auto j_vnode_index = nc_side_vnode_numbers[i];

      if (!pbdry_vnode_number_to_matched_vnode_number_set.contains(i_vnode_index))
      {
        pbdry_vnode_number_to_matched_vnode_number_set.emplace(i_vnode_index, std::set<int>());
      }

      if (!pbdry_vnode_number_to_matched_vnode_number_set.contains(j_vnode_index))
      {
        pbdry_vnode_number_to_matched_vnode_number_set.emplace(j_vnode_index, std::set<int>());
      }

      pbdry_vnode_number_to_matched_vnode_number_set.at(i_vnode_index).insert(j_vnode_index);
      pbdry_vnode_number_to_matched_vnode_number_set.at(j_vnode_index).insert(i_vnode_index);
    }
  }

  // consider pbdry conner
  const auto dim = this->dimension();
  for (int i = 0; i < dim - 1; ++i)
  {
    for (auto& [pbdry_vnode_number, matched_vnode_number_set] : pbdry_vnode_number_to_matched_vnode_number_set)
    {
      if (matched_vnode_number_set.size() == 1)
      {
        continue;
      }

      for (const auto matched_vnode_index : matched_vnode_number_set)
      {
        const auto& other_matched_vnode_index_set = pbdry_vnode_number_to_matched_vnode_number_set.at(matched_vnode_index);

        auto&       i_set = matched_vnode_number_set;
        const auto& j_set = other_matched_vnode_index_set;

        const auto difference = set_difference(j_set, i_set);

        if (difference.empty())
        {
          continue;
        }

        i_set.insert(difference.begin(), difference.end());
        i_set.erase(pbdry_vnode_number);
      }
    }
  }

  return pbdry_vnode_number_to_matched_vnode_number_set;
}

} // namespace ms::grid