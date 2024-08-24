#include "msgrid/Grid.h"

#include "msexception/Exception.h"
#include "msgeo/Figure.h"
#include "msgrid/Grid_Data.h"
#include "msmath/Vector.h"
#include <algorithm>
#include <functional>
#include <map>
#include <set>
#include <unordered_set>

// anonymous function definition
namespace
{

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

template <typename Container1, typename Container2, std::enable_if_t<std::is_same_v<typename Container1::value_type, typename Container2::value_type>, bool> = true>
void set_intersection(std::vector<typename Container1::value_type>& intersection, const Container1& container1, const Container2& container2)
{
  std::set_intersection(container1.begin(), container1.end(), container2.begin(), container2.end(), std::back_inserter(intersection));
}

template <typename Container1, typename Container2, std::enable_if_t<std::is_same_v<typename Container1::value_type, typename Container2::value_type>, bool> = true>
void set_intersection(std::set<typename Container1::value_type>& intersection, const Container1& container1, const Container2& container2)
{
  std::set_intersection(container1.begin(), container1.end(), container2.begin(), container2.end(), std::inserter(intersection, intersection.begin()));
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
  if (num_container == 1)
  {
    const auto& container = containers[0];

    std::vector<typename Container::value_type> result;
    result.insert(result.end(), container->begin(), container->end());

    return result;
  }

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
{
  this->make_grid_nodes(std::move(data.nodes_data));
  this->make_cell_and_boundary_elements(std::move(data.element_datas));
  this->make_periodic_boundary_elements(std::move(data.periodic_datas));

  // precalculate _node_number_to_share_cell_index_set_ignore_pbdry
  const auto num_cell = this->_cell_elements.size();

  for (int i = 0; i < num_cell; ++i)
  {
    const auto vnode_numbers = this->_cell_elements[i].vertex_node_numbers();
    for (const auto vnode_number : vnode_numbers)
    {
      if (!this->_vnode_number_to_share_cell_index_set_ignore_pbdry.contains(vnode_number))
      {
        this->_vnode_number_to_share_cell_index_set_ignore_pbdry.emplace(vnode_number, std::set<int>());
      }

      this->_vnode_number_to_share_cell_index_set_ignore_pbdry.at(vnode_number).insert(i);
    }
  }

  // precalculate _vnode_number_to_share_cell_number_set_consider_pbdry
  this->_vnode_number_to_share_cell_index_set_consider_pbdry = this->_vnode_number_to_share_cell_index_set_ignore_pbdry;

  const auto pbdry_vnode_number_to_matched_vnode_number_set = this->peridoic_boundary_vertex_node_number_to_matched_vertex_node_number_set();

  for (const auto& [pbdry_vnode_number, matched_vnode_number_set] : pbdry_vnode_number_to_matched_vnode_number_set)
  {
    for (const auto matched_vnode_number : matched_vnode_number_set)
    {
      auto&       i_set = this->_vnode_number_to_share_cell_index_set_consider_pbdry.at(pbdry_vnode_number);
      const auto& j_set = this->_vnode_number_to_share_cell_index_set_consider_pbdry.at(matched_vnode_number);

      const auto difference = set_difference(j_set, i_set);

      i_set.insert(difference.begin(), difference.end());
    }
  }

  this->make_inter_cell_face_elements();
}

int Grid::boundary_owner_cell_index(const int boundary_index) const
{
  const auto vnode_numbers = this->_boundary_elements[boundary_index].vertex_node_numbers();
  const auto cell_indexes  = this->find_cell_indexes_have_these_nodes_ignore_pbdry(vnode_numbers);
  REQUIRE(cell_indexes.size() == 1, "boundary should have unique owner cell");

  return cell_indexes.front();
}

void Grid::boundary_outward_unit_normal_at_center(double* normal, const int boundary_index, const int owner_cell_index) const
{
  const auto& bdry_element  = this->_boundary_elements[boundary_index];
  const auto& oc_element    = this->_cell_elements[owner_cell_index];
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

std::pair<int, int> Grid::inter_cell_face_owner_neighbor_cell_index_pair(const int inter_cell_face_index) const
{
  const auto num_inter_cell_face = static_cast<int>(this->_inter_cell_face_elements.size());

  if (inter_cell_face_index < num_inter_cell_face)
  {
    const auto& element       = this->_inter_cell_face_elements.at(inter_cell_face_index);
    const auto  vnode_numbers = element.vertex_node_numbers();
    const auto  cell_indexes  = this->find_cell_indexes_have_these_nodes_ignore_pbdry(vnode_numbers);
    REQUIRE(cell_indexes.size() == 2, "inter cell face should have an unique owner neighbor cell pair");

    // set first number as owner cell number
    const auto oc_index = cell_indexes[0];
    const auto nc_index = cell_indexes[1];
    return {oc_index, nc_index};
  }
  else
  {
    const auto pbdry_pair_number = inter_cell_face_index - num_inter_cell_face;
    // set first element as owner cell side element
    const auto& [oc_side_element, nc_side_element] = this->_periodic_boundary_element_pairs.at(pbdry_pair_number);

    const auto oc_indexes = this->find_cell_indexes_have_these_nodes_ignore_pbdry(oc_side_element.vertex_node_numbers());
    const auto nc_indexes = this->find_cell_indexes_have_these_nodes_ignore_pbdry(nc_side_element.vertex_node_numbers());
    REQUIRE(oc_indexes.size() == 1, "periodic boundary should have unique owner cell");
    REQUIRE(nc_indexes.size() == 1, "periodic boundary should have unique neighbor cell");

    const auto oc_index = oc_indexes.front();
    const auto nc_index = nc_indexes.front();
    return {oc_index, nc_index};
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

ms::geo::Geometry_Consisting_Nodes_Info Grid::make_discrete_partition_grid_nodes_info(const int partition_order) const
{
  auto  discrete_partition_grid_nodes_info = ms::geo::Geometry_Consisting_Nodes_Info();
  auto& pgrid_nodes                        = discrete_partition_grid_nodes_info.nodes;
  auto& pgrid_connectivites                = discrete_partition_grid_nodes_info.connectivities;

  auto connectivity_start_number = 0;
  for (const auto& cell_elem : _cell_elements)
  {
    const auto& geometry        = cell_elem.get_geometry();
    auto        pgeo_nodes_info = geometry.make_partitioned_geometry_node_info(partition_order);

    auto& pgeo_nodes         = pgeo_nodes_info.nodes;
    auto& pgeo_connectivites = pgeo_nodes_info.connectivities;

    pgrid_nodes.add_nodes(std::move(pgeo_nodes));

    for (auto& connectivity : pgeo_connectivites)
    {
      connectivity.add_number(connectivity_start_number);
    }
    pgrid_connectivites.insert(pgrid_connectivites.end(), pgeo_connectivites.begin(), pgeo_connectivites.end());

    connectivity_start_number += pgrid_nodes.num_nodes();
  }

  return discrete_partition_grid_nodes_info;
}

ms::geo::Geometry_Consisting_Nodes_Info Grid::make_partition_grid_nodes_info(const int partition_order) const
{
  using Node_2_Index = std::map<ms::geo::Node_View, int, ms::geo::Node_Compare>;

  //

  auto  partition_grid_nodes_info = ms::geo::Geometry_Consisting_Nodes_Info();
  auto& pgrid_nodes               = partition_grid_nodes_info.nodes;
  auto& pgrid_connectivites       = partition_grid_nodes_info.connectivities;

  const auto                num_cells = this->num_cells();
  std::vector<Node_2_Index> pgrid_node_to_index_s(num_cells);

  for (int i = 0; i < num_cells; ++i)
  {
    const auto& cell_elem          = this->_cell_elements[i];
    const auto& cell_geo           = cell_elem.get_geometry();
    auto        pgeo_nodes_info    = cell_geo.make_partitioned_geometry_node_info(partition_order);
    auto&       pgeo_connectivites = pgeo_nodes_info.connectivities;
    const auto& pgeo_nodes         = pgeo_nodes_info.nodes;
    const auto  num_pgeo_nodes     = pgeo_nodes.num_nodes();

    std::vector<int> pgrid_indexes(num_pgeo_nodes);

    const auto cell_index_set_to_check = this->find_face_sharing_cell_index_set_ignore_pbdry_and_less_then(i);

    for (int j = 0; j < num_pgeo_nodes; ++j)
    {
      const auto pgeo_node_view = pgeo_nodes[j]; 

      // check pgeo node view already in pgrid nodes
      // if it is in pgrid nodes, pgrid_indexes set as a pgrid index
      // if not, pgrid_indexes set as a new number and pgeo node view added in pgrid nodes.
      bool is_in_pgrid_nodes = false;
      for (const auto cell_index_to_check : cell_index_set_to_check)
      {
        const auto& pgrid_node_to_index_for_check = pgrid_node_to_index_s[cell_index_to_check];

        const auto iter = pgrid_node_to_index_for_check.find(pgeo_node_view);
        if (iter != pgrid_node_to_index_for_check.end())
        {
          is_in_pgrid_nodes = true;
          pgrid_indexes[j]  = iter->second;
          break;
        }
      }

      if (!is_in_pgrid_nodes)
      {
        pgrid_indexes[j] = pgrid_nodes.num_nodes();
        pgrid_nodes.add_node(pgeo_node_view);
      }

      // update pgrid_node_to_pgrid_index
      const auto pgrid_node_index = pgrid_indexes[j];
      const auto pgrid_node_view  = pgrid_nodes[pgrid_node_index];

      auto& pgrid_node_to_index = pgrid_node_to_index_s[i];
      pgrid_node_to_index.insert({pgrid_node_view, pgrid_node_index});
    }

    // pgeo_connectivities -> pgrid_connectivities
    for (auto& connectivity : pgeo_connectivites)
    {
      for (auto& old_index : connectivity)
      {
        old_index = pgrid_indexes[old_index];
      }
    }

    pgrid_connectivites.insert(pgrid_connectivites.end(), pgeo_connectivites.begin(), pgeo_connectivites.end());
  }

  return partition_grid_nodes_info;
}

void Grid::make_grid_nodes(Grid_Nodes_Data&& nodes_data)
{
  this->_grid_nodes = std::move(nodes_data.nodes);

  const auto num_nodes = this->_grid_nodes.num_nodes();
  this->_node_number_to_index.reserve(num_nodes);

  for (int i = 0; i < num_nodes; ++i)
  {
    const auto node_number = nodes_data.numbers[i];
    this->_node_number_to_index.emplace(node_number, i);
  }
}

void Grid::make_cell_and_boundary_elements(std::vector<Grid_Element_Data>&& element_datas)
{
  const auto num_element_datas = element_datas.size();
  this->_cell_elements.reserve(num_element_datas);
  this->_boundary_elements.reserve(num_element_datas);

  for (auto& element_data : element_datas)
  {
    auto element = this->make_element(element_data);

    if (element_data.type == Element_Type::CELL)
    {
      this->_cell_elements.push_back(std::move(element));
      // this->_cell_number_to_index.emplace(number, this->_cell_elements.size() - 1);
    }
    else
    {
      this->_boundary_elements.push_back(std::move(element));
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

  std::map<std::vector<double>, std::vector<Element>> direction_to_periodic_elements;

  // make element and sort by direction
  for (auto& [element_data, periodic_direction] : periodic_datas)
  {
    auto element = this->make_element(element_data);

    bool       is_new_direction  = true;
    const auto direction_vector1 = ms::math::Vector_View(periodic_direction);

    for (auto& [periodic_direction2, elements] : direction_to_periodic_elements)
    {
      const auto direction_vector2 = ms::math::Vector_View(periodic_direction2);

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
      direction_to_periodic_elements.at(periodic_direction).push_back(std::move(element));
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
        auto  matched_node_numbers = j_element.find_periodic_matched_node_numbers(direction, i_element);

        if (matched_node_numbers.empty()) continue;

        // Reordering is necessary for the quadrature points of the two elements to match
        j_element.reordering_nodes(matched_node_numbers);

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
  // list up boundary vertex node numbers
  std::set<std::vector<int>> bdry_vnode_numbers_set;
  for (const auto& boundary_element : this->_boundary_elements)
  {
    bdry_vnode_numbers_set.insert(boundary_element.make_sorted_vertex_node_numbers());
  }

  // list up periodic boundary vertex node numbers
  std::set<std::vector<int>> pbdry_vnode_numbers_set;
  for (const auto& [pbdry_elem1, pbdry_elem2] : this->_periodic_boundary_element_pairs)
  {
    pbdry_vnode_numbers_set.insert(pbdry_elem1.make_sorted_vertex_node_numbers());
    pbdry_vnode_numbers_set.insert(pbdry_elem2.make_sorted_vertex_node_numbers());
  }

  const auto                              num_cells = this->num_cells();
  std::vector<std::set<std::vector<int>>> constructed_face_vnode_numbers_sets(num_cells);

  for (int i = 0; i < num_cells; ++i)
  {
    const auto& cell_element = this->_cell_elements[i];
    const auto& geometry     = cell_element.get_geometry();
    const auto  num_face     = geometry.num_faces();

    const auto cell_index_set_to_check = this->find_face_sharing_cell_index_set_ignore_pbdry_and_less_then(i);

    for (int j = 0; j < num_face; ++j)
    {
      auto face_sorted_vnode_numbers = cell_element.make_sorted_face_vertex_node_numbers(j);

      if (bdry_vnode_numbers_set.contains(face_sorted_vnode_numbers)) continue;
      if (pbdry_vnode_numbers_set.contains(face_sorted_vnode_numbers)) continue;

      bool is_constructed_face = false;

      for (const auto cell_index : cell_index_set_to_check)
      {
        if (constructed_face_vnode_numbers_sets[cell_index].contains(face_sorted_vnode_numbers))
        {
          is_constructed_face = true;
          break;
        }
      }

      if (is_constructed_face) continue;

      // construct intercell face
      this->_inter_cell_face_elements.push_back(cell_element.make_face_element(j));
      constructed_face_vnode_numbers_sets[i].insert(std::move(face_sorted_vnode_numbers));
    }
  }
}

Element Grid::make_element(Grid_Element_Data& element_data) const
{
  auto& [number, type, figure, node_numbers] = element_data;

  const auto geo_figure = convert(figure);

  const auto num_nodes = node_numbers.size();

  std::vector<ms::geo::Node_View> node_views(num_nodes);
  for (int i = 0; i < num_nodes; ++i)
  {
    auto& node_view = node_views[i];

    const auto node_number = node_numbers[i];
    const auto index       = this->_node_number_to_index.at(node_number);
    node_view              = this->_grid_nodes[index];
  }

  auto geometry = ms::geo::Geometry(geo_figure, std::move(node_views));
  auto element  = Element(type, std::move(node_numbers), std::move(geometry));
  return element;
}

std::vector<int> Grid::find_cell_indexes_have_these_nodes_ignore_pbdry(const std::span<const int> vnode_numbers) const
{
  const auto num_vnode = vnode_numbers.size();

  std::vector<const std::set<int>*> share_cell_index_set_ptrs;
  share_cell_index_set_ptrs.reserve(num_vnode);

  for (int i = 0; i < num_vnode; ++i)
  {
    share_cell_index_set_ptrs.push_back(&this->_vnode_number_to_share_cell_index_set_ignore_pbdry.at(vnode_numbers[i]));
  }

  return set_intersection(share_cell_index_set_ptrs);
}

std::unordered_map<int, std::set<int>> Grid::peridoic_boundary_vertex_node_number_to_matched_vertex_node_number_set(void) const
{
  std::unordered_map<int, std::set<int>> pbdry_vnode_number_to_matched_vnode_number_set;

  // Grid::make_periodic_boundary_elements에서 periodic boundary element pair를 만들 때, vertex node끼리 정렬되게 만들었음을 가정한다.
  for (const auto& [oc_side_element, nc_side_element] : this->_periodic_boundary_element_pairs)
  {
    const auto oc_side_vnode_numbers = oc_side_element.vertex_node_numbers();
    const auto nc_side_vnode_numbers = nc_side_element.vertex_node_numbers();

    const auto num_vnode = oc_side_vnode_numbers.size();
    for (int i = 0; i < num_vnode; ++i)
    {
      const auto i_vnode_number = oc_side_vnode_numbers[i];
      const auto j_vnode_number = nc_side_vnode_numbers[i];

      if (!pbdry_vnode_number_to_matched_vnode_number_set.contains(i_vnode_number))
      {
        pbdry_vnode_number_to_matched_vnode_number_set.emplace(i_vnode_number, std::set<int>());
      }

      if (!pbdry_vnode_number_to_matched_vnode_number_set.contains(j_vnode_number))
      {
        pbdry_vnode_number_to_matched_vnode_number_set.emplace(j_vnode_number, std::set<int>());
      }

      pbdry_vnode_number_to_matched_vnode_number_set.at(i_vnode_number).insert(j_vnode_number);
      pbdry_vnode_number_to_matched_vnode_number_set.at(j_vnode_number).insert(i_vnode_number);
    }
  }

  // consider pbdry conner
  const auto dim = this->dimension();
  for (int i = 0; i < dim - 1; ++i)
  {
    for (auto& [pbdry_vnode_number, matched_vnode_number_set] : pbdry_vnode_number_to_matched_vnode_number_set)
    {
      if (matched_vnode_number_set.size() == 1) continue;

      for (const auto matched_vnode_number : matched_vnode_number_set)
      {
        const auto& other_matched_vnode_number_set = pbdry_vnode_number_to_matched_vnode_number_set.at(matched_vnode_number);

        auto&       i_set      = matched_vnode_number_set;
        const auto& j_set      = other_matched_vnode_number_set;
        const auto  difference = set_difference(j_set, i_set);

        if (difference.empty()) continue;

        i_set.insert(difference.begin(), difference.end());
        i_set.erase(pbdry_vnode_number);
      }
    }
  }

  return pbdry_vnode_number_to_matched_vnode_number_set;
}

std::set<int> Grid::find_face_sharing_cell_index_set_ignore_pbdry(const int cell_index) const
{
  const auto& cell_element = this->_cell_elements[cell_index];
  const auto& geometry     = cell_element.get_geometry();
  const auto  num_faces    = geometry.num_faces();

  std::set<int> face_sharing_cell_index_set;

  for (int i = 0; i < num_faces; ++i)
  {
    const auto face_vnode_numbers = cell_element.face_vertex_node_numbers(i);
    auto       cell_indexes       = this->find_cell_indexes_have_these_nodes_ignore_pbdry(face_vnode_numbers);

    auto iter = std::find(cell_indexes.begin(), cell_indexes.end(), cell_index);
    REQUIRE(iter != cell_indexes.end(), "face_sharing_cell_index should be in cell_indexes");
    cell_indexes.erase(iter);

    for (const auto face_sharing_cell_index : cell_indexes)
    {
      face_sharing_cell_index_set.insert(face_sharing_cell_index);
    }
  }

  return face_sharing_cell_index_set;
}

std::set<int> Grid::find_face_sharing_cell_index_set_ignore_pbdry_and_less_then(const int cell_index) const
{
  const auto& cell_element = this->_cell_elements[cell_index];
  const auto& geometry     = cell_element.get_geometry();
  const auto  num_faces    = geometry.num_faces();

  std::set<int> face_sharing_cell_index_set;

  for (int i = 0; i < num_faces; ++i)
  {
    const auto& face_vnode_numbers = cell_element.face_vertex_node_numbers(i);
    auto        cell_indexes       = this->find_cell_indexes_have_these_nodes_ignore_pbdry(face_vnode_numbers);

    auto iter = std::find(cell_indexes.begin(), cell_indexes.end(), cell_index);
    REQUIRE(iter != cell_indexes.end(), "face_sharing_cell_index should be in cell_indexes");
    cell_indexes.erase(iter);

    for (const auto face_sharing_cell_index : cell_indexes)
    {
      if (face_sharing_cell_index < cell_index)
      {
        face_sharing_cell_index_set.insert(face_sharing_cell_index);
      }
    }
  }

  return face_sharing_cell_index_set;
}

} // namespace ms::grid