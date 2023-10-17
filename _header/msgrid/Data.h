#pragma once
#include <vector>

#include "msgeo/Node.h"

namespace ms::grid
{

enum class Element_Type
{
  CELL,
  FACE,
  SLIP_WALL,
  SUPERSONIC_INLET,
  SUPERSONIC_OUTLET,
  INTIAL_CONSTANT,
  PERIODIC,
  NOT_IN_LIST
};

enum class Figure
{
  POINT         = 0,
  LINE          = 1,
  TRIANGLE      = 2,
  QUADRILATERAL = 3,
  TETRAHEDRAL   = 4,
  HEXAHEDRAL    = 5,
  PRISM         = 6,
  PYRAMID       = 7,
  NUM_FIGURES,
  NOT_IN_LIST = -1
};

struct Grid_Nodes_Data
{
  ms::geo::Nodes   nodes;
  std::vector<int> numbers;
};

// 2.node numbers���� vertex node number�� ���� ���´�.
//  2-1. vertex node number�� vertex�� �ݽð� �������� ���� ������� ���´�.
struct Grid_Element_Data
{
  int              number;
  Element_Type     type;
  Figure           figure;
  std::vector<int> node_numbers;
};

struct Grid_Peridoic_Data
{
  Grid_Element_Data   element_data;
  std::vector<double> periodic_direction;
};

struct Grid_Data
{
  Grid_Nodes_Data                 nodes_data;
  std::vector<Grid_Element_Data>  element_datas;
  std::vector<Grid_Peridoic_Data> periodic_datas;
};

} // namespace ms::grid