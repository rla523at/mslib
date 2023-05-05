#include "mstecplot/Writer.h"
#include "mstecplot/Data.h"
#include "gtest/gtest.h"

TEST(msTecplot, ASCII_Grid_Test1)
{
  ms::tecplot::Configuration configuration;
  configuration.file_fomrat = ms::tecplot::File_Format::ASCII;
  configuration.folder_path = "TEST/";
  configuration.title       = "ASCII_GRID_TEST1";

  ms::tecplot::Writer::set_configuration(std::move(configuration));

  constexpr auto                dimension      = 2;
  constexpr auto                num_nodes      = 4;
  constexpr auto                num_elements   = 2;
  std::vector<double>           coordinates    = {1, 2, 1, 2, 1, 1, 2, 2};
  std::vector<std::vector<int>> connectivities = {{1, 2, 4}, {1, 4, 3}};

  ms::tecplot::Grid_Data grid_data(dimension, num_nodes, num_elements, coordinates.data(), connectivities);
  ms::tecplot::Writer::write_grid_file(grid_data);
}

TEST(msTecplot, ASCII_Solution_Test1)
{
  ms::tecplot::Configuration configuration;
  configuration.file_fomrat = ms::tecplot::File_Format::ASCII;
  configuration.folder_path = "TEST/";
  configuration.title       = "ASCII_SOLUTION_TEST1";

  ms::tecplot::Writer::set_configuration(std::move(configuration));

  constexpr auto                dimension      = 2;
  constexpr auto                num_nodes      = 4;
  constexpr auto                num_elements   = 2;
  std::vector<double>           coordinates    = {1, 2, 1, 2, 1, 1, 2, 2};
  std::vector<std::vector<int>> connectivities = {{1, 2, 4}, {1, 4, 3}};

  ms::tecplot::Grid_Data grid_file_data(dimension, num_nodes, num_elements, coordinates.data(), connectivities);
  ms::tecplot::Writer::write_grid_file(grid_file_data);

  std::vector<double>                         variables      = {1, 2, 3, 4, -1, -2, 11, 12, 13, 14, 20, 30};
  constexpr auto                              solution_times = 1.0;
  std::vector<std::string>                    variable_strs  = {"A", "B", "C", "D"};
  std::vector<ms::tecplot::Variable_Location> var_locations  = {ms::tecplot::Variable_Location::NODAL, ms::tecplot::Variable_Location::CELLCENTERED, ms::tecplot::Variable_Location::NODAL, ms::tecplot::Variable_Location::CELLCENTERED};

  ms::tecplot::Solution_Data solution_data;
  solution_data.values_ptr        = variables.data();
  solution_data.solution_time     = solution_times;
  solution_data.variable_strs     = variable_strs;
  solution_data.var_locations_ptr = var_locations.data();

  ms::tecplot::Writer::write_solution_file(solution_data);
}

TEST(msTecplot, Binary_Grid_Test1)
{
  ms::tecplot::Configuration configuration;
  configuration.file_fomrat = ms::tecplot::File_Format::Binary;
  configuration.folder_path = "TEST/";
  configuration.title       = "Binary_GRID_TEST1";

  ms::tecplot::Writer::set_configuration(std::move(configuration));

  constexpr auto                dimension      = 2;
  constexpr auto                num_nodes      = 4;
  constexpr auto                num_elements   = 2;
  std::vector<double>           coordinates    = {1, 2, 1, 2, 1, 1, 2, 2};
  std::vector<std::vector<int>> connectivities = {{1, 2, 4}, {1, 4, 3}};

  ms::tecplot::Grid_Data grid_data(dimension, num_nodes, num_elements, coordinates.data(), connectivities);
  ms::tecplot::Writer::write_grid_file(grid_data);
}

TEST(msTecplot, Binary_Solution_Test1)
{
  ms::tecplot::Configuration configuration;
  configuration.file_fomrat = ms::tecplot::File_Format::Binary;
  configuration.folder_path = "TEST/";
  configuration.title       = "Binary_Solution_TEST1";

  ms::tecplot::Writer::set_configuration(std::move(configuration));

  // grid
  constexpr auto                dimension      = 2;
  constexpr auto                num_nodes      = 4;
  constexpr auto                num_elements   = 2;
  std::vector<double>           coordinates    = {1, 2, 1, 2, 1, 1, 2, 2};
  std::vector<std::vector<int>> connectivities = {{1, 2, 4}, {1, 4, 3}};

  ms::tecplot::Grid_Data grid_file_info(dimension, num_nodes, num_elements, coordinates.data(), connectivities);
  ms::tecplot::Writer::write_grid_file(grid_file_info);

  // solution
  std::vector<double>                         variables      = {1, 2, 3, 4, -1, -2, 11, 12, 13, 14, 20, 30};
  constexpr auto                              solution_times = 1.0;
  std::vector<std::string>                    variable_strs  = {"A", "B", "C", "D"};
  std::vector<ms::tecplot::Variable_Location> var_locations  = {ms::tecplot::Variable_Location::NODAL, ms::tecplot::Variable_Location::CELLCENTERED, ms::tecplot::Variable_Location::NODAL, ms::tecplot::Variable_Location::CELLCENTERED};

  ms::tecplot::Solution_Data solution_data;
  solution_data.values_ptr        = variables.data();
  solution_data.solution_time     = solution_times;
  solution_data.variable_strs     = variable_strs;
  solution_data.var_locations_ptr = var_locations.data();

  ms::tecplot::Writer::write_solution_file(solution_data);
}