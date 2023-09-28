#include "msgrid/Gmsh_Reader.h"
#include "msgrid/Data.h"
#include "gtest/gtest.h"

namespace GFR_TEST_API
{
void is_equal(const ms::grid::Grid_Nodes_Data& data1, const ms::grid::Grid_Nodes_Data& data2)
{
  EXPECT_EQ(data1.coordinates, data2.coordinates);
  EXPECT_EQ(data1.dimension, data2.dimension);
  EXPECT_EQ(data1.num_nodes, data2.num_nodes);
  EXPECT_EQ(data1.type, data2.type);
}

void is_equal(const ms::grid::Grid_Element_Data& data1, const ms::grid::Grid_Element_Data& data2)
{
  EXPECT_EQ(data1.node_numberss, data2.node_numberss);
  EXPECT_EQ(data1.figure, data2.figure);
  EXPECT_EQ(data1.element_type, data2.element_type);
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

  const auto num_nodes = grid_file_data.nodes_data.num_nodes;

  constexpr auto ref_num_nodes = 505;
  EXPECT_EQ(num_nodes, ref_num_nodes);

  // const auto&         coordinates      = grid_file_data.nodes_data.coordinates;
  // std::vector<double> ref_coordinates = {0.4299999999999999, 0.01, 0};

  // EXPECT_EQ(coordinates, ref_coordinates);
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
  ref_elem_data.element_type = Element_Type::CELL;
  ref_elem_data.node_numberss = {42, 43, 144, 143}; // 1¾¿ »©¼­ ÀÐÀ½

  GFR_TEST_API::is_equal(ref_elem_data, elem_data);
}

} // namespace ms::grid
