#pragma once
#include <string>
#include <vector>

namespace ms::tecplot
{

enum class Variable_Location
{
  NODAL,
  CELLCENTERED
};

// 1. coordinates_ptr format :  dimension X num nodes  (x1,...,xn,y1,...,yn,z1,...,zn)
// 2. connectivities start from 1
struct Grid_Data
{
public:
  Grid_Data(int dimension, int num_nodes, int num_elements, const double* coordinates_ptr, const std::vector<std::vector<int>>& connectivities);

public:
  int                                  dimension;
  int                                  num_nodes;
  int                                  num_elements;
  const double*                        coordinates_ptr;
  const std::vector<std::vector<int>>& connectivities;
};

// 1. values_ptr format
// num variables X num_nodes ((v1)1, ..., (v1)nn , ... , (vn)1,...,(vn)nn) (variable location = NODAL)
// num variables X num_cells ((v1)1, ..., (v1)nc , ... , (vn)1,...,(vn)nc) (variable location = CELLCENTERD)
struct Solution_Data
{
public:
  std::string var_location_str(void) const;
  int         num_variables(void) const;

public:
  const double*            values_ptr;
  double                   solution_time;
  std::vector<std::string> variable_strs;
  const Variable_Location* var_locations_ptr;
};

} // namespace ms::tecplot