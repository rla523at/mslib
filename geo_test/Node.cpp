#include "gtest/gtest.h"
#include "msgeo/Node.h"

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