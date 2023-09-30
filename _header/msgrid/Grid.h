#pragma once
#include "Data.h"
#include "Element.h"

#include "msgeo/Geometry.h"
#include "msgeo/Node.h"
#include <set>
#include <unordered_map>

// class declaration
namespace ms::grid
{

class Grid
{
public:
  Grid(Grid_Data&& data);

public:
  int                                     boundary_owner_cell_number(const int boundary_number) const;
  void                                    boundary_outward_unit_normal_at_center(double* normal, const int boundary_number, const int owner_cell_number) const;
  Element_Type                            boundary_type(const int boundary_number) const;
  double                                  boundary_volume(const int boundary_number) const;
  void                                    cell_center(double* cell_center, const int cell_number) const;
  void                                    cell_projected_volume(double* cell_projected_volume, const int cell_number) const;
  double                                  cell_volume(const int cell_number) const;
  int                                     dimension(void) const;
  bool                                    is_line_cell(const int cell_number) const;
  std::pair<int, int>                     inter_cell_face_owner_neighbor_cell_number_pair(const int inter_cell_face_number) const;
  void                                    inter_cell_face_outward_unit_normal_at_center(double* normal, const int inter_cell_face_number, const int owner_cell_number) const;
  double                                  inter_cell_face_volume(const int inter_cell_face_number) const;
  int                                     num_boundaries(void) const;
  int                                     num_inter_cell_faces(void) const;
  int                                     num_cells(void) const;
  ms::geo::Geometry_Consisting_Nodes_Info make_discrete_partition_grid_nodes_info(const int partition_order) const;
  ms::geo::Geometry_Consisting_Nodes_Info make_partition_grid_nodes_info(const int partition_order) const;

private:
  void    make_grid_nodes(Grid_Nodes_Data&& nodes_data);
  void    make_cell_and_boundary_elements(std::vector<Grid_Element_Data>&& element_datas);
  void    make_periodic_boundary_elements(std::vector<Grid_Peridoic_Data>&& periodic_datas);
  void    make_inter_cell_face_elements(void);
  Element make_element(Grid_Element_Data& element_data) const;

private:
  std::vector<int>                       find_cell_indexes_have_these_nodes_ignore_pbdry(const std::span<const int> vnode_numbers) const;
  std::unordered_map<int, std::set<int>> peridoic_boundary_vertex_node_number_to_matched_vertex_node_number_set(void) const;
  std::set<int>                          find_face_sharing_cell_index_set_ignore_pbdry(const int cell_index) const;

private:
  ms::geo::Nodes                           _grid_nodes;
  std::vector<Element>                     _cell_elements;
  std::vector<Element>                     _boundary_elements;
  std::vector<Element>                     _inter_cell_face_elements;
  std::vector<std::pair<Element, Element>> _periodic_boundary_element_pairs;
  std::unordered_map<int, int>             _node_number_to_index;

  // 이것들이 필요한지 잘 모르겠네..
  // std::unordered_map<int, int>             _cell_number_to_index;
  // std::unordered_map<int, int>             _boundary_number_to_index;
  // std::unordered_map<int, int>             _inter_cell_face_number_to_index;
  // std::unordered_map<int, int>             _periodic_boundary_number_to_index;

  std::unordered_map<int, std::set<int>> _node_number_to_share_cell_index_set_ignore_pbdry;
  std::unordered_map<int, std::set<int>> _node_number_to_share_cell_index_set_consider_pbdry;
  // std::unordered_map<int, std::set<int>> _cell_index_to_share_cell_index_set_ignore_pbdry;
};

} // namespace ms::grid