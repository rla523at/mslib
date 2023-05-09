#pragma once
#include "Element.h"
#include "msgeo/Node.h"

// forward declaration
namespace ms::grid
{
struct Data;
struct Element_Data;
struct Peridoic_Data;
} // namespace ms::grid

/*





*/

// class declaration
namespace ms::grid
{

class Grid
{
public:
  Grid(Data&& data);

private:
  void make_cell_and_boundary_elements(std::vector<Element_Data>&& element_datas);
  void make_periodic_boundary_elements(std::vector<Peridoic_Data>&& periodic_datas);

private:
  ms::geo::Nodes                           _nodes;
  std::vector<Element>                     _cell_elements;
  std::vector<Element>                     _boundary_elements;
  std::vector<Element>                     _inter_cell_face_elements;
  std::vector<std::pair<Element, Element>> _periodic_boundary_element_pairs;
};

} // namespace ms::grid