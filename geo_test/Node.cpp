#include "gtest/gtest.h"
#include "msgeo/Node.h"
#include "msmath/Vector.h"

TEST(Nodes_Const_Wrapper, Block_Nodes)
{
  const auto type = ms::geo::Coordinates_Type::BLOCK;
  const auto num_nodes = 3;
  const auto dimension = 2;
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
  ms::geo::Node_Const_Wrapper node2(dimension, coordinates.data()+dimension);

  ms::math::Vector<3> v;
  node1.other_to_this_vector(node2, v);

  ms::math::Vector<3> ref = {-2, 0, 2};
  EXPECT_EQ(v, ref);
}