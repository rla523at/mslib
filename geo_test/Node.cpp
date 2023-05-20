#include "msgeo/Node.h"
#include "msmath/Vector.h"
#include "mssym/Polynomial.h"
#include "gtest/gtest.h"

TEST(Nodes_Const_Wrapper, Block_Nodes)
{
  const auto          type        = ms::geo::Coordinates_Type::BLOCK;
  const auto          num_nodes   = 3;
  const auto          dimension   = 2;
  std::vector<double> coordinates = {1, 2, 3, 3, 2, 1};

  ms::geo::Nodes_Const_Wrapper nodes(type, num_nodes, dimension, coordinates.data());

  const auto& node1 = nodes[1];
  EXPECT_DOUBLE_EQ(node1[0], 2);
  EXPECT_DOUBLE_EQ(node1[1], 2);
}
TEST(Nodes_Const_Wrapper, other_node_to_this_node)
{
  constexpr auto      dimension   = 3;
  std::vector<double> coordinates = {1, 2, 3, 3, 2, 1};

  ms::geo::Node_Const_Wrapper node1(dimension, coordinates.data());
  ms::geo::Node_Const_Wrapper node2(dimension, coordinates.data() + dimension);

  ms::math::Vector<3> v;
  node1.other_to_this_vector(node2, v);

  ms::math::Vector<3> ref = {-2, 0, 2};
  EXPECT_EQ(v, ref);
}
TEST(Nodes_Const_Wrapper, conversion1)
{
  constexpr auto               dimension   = 3;
  std::vector<double>          coordinates = {1, 2, 3, 4, 5, 6};
  ms::geo::Nodes_Const_Wrapper nodes(ms::geo::Coordinates_Type::BLOCK, 2, dimension, coordinates.data());

  ms::sym::Polynomial x("x0");
  ms::sym::Polynomial y("x1");
  ms::sym::Polynomial z("x2");

  auto p = x * x + y + y * z;

  const auto     val1 = p(nodes[0]);
  constexpr auto ref1 = 1 * 1 + 3 + 3 * 5;
  EXPECT_EQ(val1, ref1);

  const auto     val2 = p(nodes[1]);
  constexpr auto ref2 = 2 * 2 + 4 + 4 * 6;
  EXPECT_EQ(val2, ref2);
}