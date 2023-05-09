#pragma once
#include <vector>

namespace ms::grid
{

enum class Coordinate_Type
{
  NODAL,
  BLOCK
};

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

struct Nodes_Data
{
public:
  const double* coordinates_ptr(const int node_index) const;
  int           stride(void) const;

public:
  Coordinate_Type     type;
  int                 dimension;
  int                 num_nodes;
  std::vector<double> coordinates;
};

struct Element_Data
{
  Element_Type     element_type;
  Figure           figure;
  std::vector<int> node_indexes;
};

struct Peridoic_Data
{
  Element_Data        element_data;
  std::vector<double> periodic_direction;
};

struct Data
{
  Nodes_Data                 nodes_data;
  std::vector<Element_Data>  element_datas;
  std::vector<Peridoic_Data> periodic_datas;
};

} // namespace ms::grid