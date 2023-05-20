#pragma once
#include "Element.h"
#include "msgeo/Node.h"

// forward declaration
namespace ms::grid
{
struct Grid_Data;
struct Grid_Element_Data;
struct Grid_Peridoic_Data;
} // namespace ms::grid

/*





*/

// class declaration
namespace ms::grid
{

class Grid
{
public:
  Grid(Grid_Data&& data);

public:
  void   cell_center(double* cell_center, const int cell_number) const;
  void   cell_projected_volume(double* cell_projected_volume, const int cell_number) const;
  double cell_volume(const int cell_number) const;
  int    dimension(void) const;
  bool   is_line_cell(const int cell_number) const;
  int    num_cells(void) const;

private:
  void make_cell_and_boundary_elements(std::vector<Grid_Element_Data>&& element_datas);
  void make_periodic_boundary_elements(std::vector<Grid_Peridoic_Data>&& periodic_datas);
  void make_inter_cell_face_elements(void);

private:
  ms::geo::Nodes                           _grid_nodes;
  std::vector<Element>                     _cell_elements;
  std::vector<Element>                     _boundary_elements;
  std::vector<Element>                     _inter_cell_face_elements;
  std::vector<std::pair<Element, Element>> _periodic_boundary_element_pairs;
};

} // namespace ms::grid