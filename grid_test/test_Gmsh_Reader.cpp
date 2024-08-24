#include "msgrid/Grid_Data.h"
#include "msgrid/Gmsh_Reader.h"
#include "gtest/gtest.h"

namespace GFR_TEST_API
{
void is_equal(const ms::grid::Grid_Nodes_Data& data1, const ms::grid::Grid_Nodes_Data& data2)
{
  EXPECT_EQ(data1.nodes, data2.nodes);
  EXPECT_EQ(data1.numbers, data2.numbers);
}

void is_equal(const ms::grid::Grid_Element_Data& data1, const ms::grid::Grid_Element_Data& data2)
{
  EXPECT_EQ(data1.node_numbers, data2.node_numbers);
  EXPECT_EQ(data1.figure, data2.figure);
  EXPECT_EQ(data1.type, data2.type);
}
} // namespace GFR_TEST_API

namespace ms::grid
{

TEST(Gmsh_Reader, read1)
{
  constexpr auto dimension = 2;
  Gmsh_Reader    grid_file_reader(dimension);

  constexpr auto file_path      = "TEST/gmsh_test1.msh";
  const auto     grid_file_data = grid_file_reader.read(file_path);

  const auto num_nodes = grid_file_data.nodes_data.nodes.num_nodes();

  constexpr auto ref_num_nodes = 505;
  EXPECT_EQ(num_nodes, ref_num_nodes);
}

TEST(Gmsh_Reader, read2)
{
  constexpr auto dimension = 2;
  Gmsh_Reader    grid_file_reader(dimension);

  constexpr auto file_path      = "TEST/gmsh_test1.msh";
  const auto     grid_file_data = grid_file_reader.read(file_path);

  const auto     num_elements     = grid_file_data.element_datas.size() + grid_file_data.periodic_datas.size();
  constexpr auto ref_num_elements = 608;
  EXPECT_EQ(num_elements, ref_num_elements);

  const auto& elem_data = grid_file_data.element_datas[42];

  Grid_Element_Data ref_elem_data;
  ref_elem_data.figure       = Figure::QUADRILATERAL;
  ref_elem_data.type         = Element_Type::CELL;
  ref_elem_data.node_numbers = {43, 44, 145, 144};

  GFR_TEST_API::is_equal(ref_elem_data, elem_data);
}

} // namespace ms::grid
