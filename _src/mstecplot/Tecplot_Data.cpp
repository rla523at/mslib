#include "mstecplot/Tecplot_Data.h"

#include "msexception/Exception.h"

namespace ms::tecplot
{

Grid_Data::Grid_Data(int dimension, int num_nodes, int num_elements, const double* coordinates_ptr, const std::vector<std::vector<int>>& connectivities)
    : dimension(dimension), num_nodes(num_nodes), num_elements(num_elements), coordinates_ptr(coordinates_ptr), connectivities(connectivities)
{
  REQUIRE(this->connectivities.size() == this->num_elements, "connectivity should be givne for all elements");
}

int Solution_Data::num_variables(void) const
{
  return static_cast<int>(this->variable_strs.size());
}

std::string Solution_Data::var_location_str(void) const
{
  const auto num_variables = this->num_variables();

  std::vector<int> cell_centered_var_indexes;
  cell_centered_var_indexes.reserve(num_variables);

  for (int i = 0; i < num_variables; ++i)
  {
    if (var_locations_ptr[i] == Variable_Location::CELL_CENTER)
    {
      const auto var_index = i + 1;
      cell_centered_var_indexes.push_back(var_index);
    }
  }

  if (cell_centered_var_indexes.empty())
  {
    return "NODE";
  }
  else
  {
    std::string var_location_str = "([";

    for (auto var_index : cell_centered_var_indexes)
    {
      var_location_str += std::to_string(var_index) + ",";
    }

    var_location_str.pop_back();
    var_location_str += "] = CELL_CENTER)";

    return var_location_str;
  }
}

} // namespace ms::tecplot


