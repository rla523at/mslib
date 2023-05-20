#include "msgrid/Grid.h"
#include "msgrid/Gmsh_Reader.h"
#include "gtest/gtest.h"

namespace ms::grid
{
TEST(Grid, construct)
{
  constexpr auto dimension = 1;
  Gmsh_Reader    grid_file_reader(dimension);

  constexpr auto file_path      = "TEST/1D.msh";
  auto           grid_file_data = grid_file_reader.read(file_path);

  Grid grid(std::move(grid_file_data));
}

} // namespace ms::grid