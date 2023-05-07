#pragma once
#include "Element.h"
#include "msgeo/Node.h"

namespace ms::grid
{
struct Data;
}

namespace ms::grid
{

class Grid
{
public:
  Grid(Data&& data);

private:
  ms::geo::Nodes                           _nodes;
  std::vector<Element>                     _cell_elements;
  std::vector<Element>                     _boundary_elements;
  std::vector<Element>                     _inter_cell_face_elements;
  std::vector<std::pair<Element, Element>> _periodic_boundary_element_pairs;
};

} // namespace ms::grid