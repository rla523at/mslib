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

TEST(Grid, make_discrete_partition_data1)
{
  constexpr auto dimension = 1;
  Gmsh_Reader    grid_file_reader(dimension);

  constexpr auto file_path      = "TEST/1D.msh";
  auto           grid_file_data = grid_file_reader.read(file_path);

  Grid       grid(std::move(grid_file_data));
  const auto discrete_partition_data = grid.make_discrete_partition_data(0);

  EXPECT_EQ(discrete_partition_data.nodes.num_nodes(), 200);
  EXPECT_EQ(discrete_partition_data.connectivities.size(), 100);
}
TEST(Grid, make_discrete_partition_data2)
{
  constexpr auto dimension = 1;
  Gmsh_Reader    grid_file_reader(dimension);

  constexpr auto file_path      = "TEST/1D.msh";
  auto           grid_file_data = grid_file_reader.read(file_path);

  Grid       grid(std::move(grid_file_data));
  const auto discrete_partition_data = grid.make_discrete_partition_data(1);

  EXPECT_EQ(discrete_partition_data.nodes.num_nodes(), 300);
  EXPECT_EQ(discrete_partition_data.connectivities.size(), 200);
}

TEST(Grid, make_partition_data1)
{
  constexpr auto dimension = 1;
  Gmsh_Reader    grid_file_reader(dimension);

  constexpr auto file_path      = "TEST/1D.msh";
  auto           grid_file_data = grid_file_reader.read(file_path);

  Grid       grid(std::move(grid_file_data));
  const auto partition_data = grid.make_partition_data(0);

  EXPECT_EQ(partition_data.nodes.num_nodes(), 101);
  EXPECT_EQ(partition_data.connectivities.size(), 100);
}
TEST(Grid, make_partition_data2)
{
  constexpr auto dimension = 1;
  Gmsh_Reader    grid_file_reader(dimension);

  constexpr auto file_path      = "TEST/1D.msh";
  auto           grid_file_data = grid_file_reader.read(file_path);

  Grid       grid(std::move(grid_file_data));
  const auto partition_data = grid.make_partition_data(1);

  EXPECT_EQ(partition_data.nodes.num_nodes(), 201);
  EXPECT_EQ(partition_data.connectivities.size(), 200);
}

} // namespace ms::grid