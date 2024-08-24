#include "msgeo/Node.h"
#include "msmath/Vector.h"
#include "mssym/Polynomial.h"
#include "gtest/gtest.h"

std::ostream& operator<<(std::ostream& os, const ms::geo::Node_View& node)
{
  return os << node.to_string();
}
std::ostream& operator<<(std::ostream& os, const ms::geo::Node& node)
{
  return os << node.to_string();
}
std::ostream& operator<<(std::ostream& os, const ms::geo::Node_Wrap& node)
{
  return os << node.to_string();
}

TEST(Node_Compare, compare1)
{
  ms::geo::Node node1({1, 2});
  ms::geo::Node node2({1, 2});

  ms::geo::Node_Compare compare;
  EXPECT_FALSE(compare(node1, node2));
}
TEST(Node_Compare, compare2)
{
  ms::geo::Node node1({1, 2});
  ms::geo::Node node2({1.1, 2});

  ms::geo::Node_Compare compare;
  EXPECT_TRUE(compare(node1, node2));
}
TEST(Node_Compare, compare3)
{
  constexpr auto epsilon = 1.0e-13;
  ms::geo::Node node1({1, 2});
  ms::geo::Node node2({1+epsilon, 2-epsilon});

  ms::geo::Node_Compare compare;
  EXPECT_FALSE(compare(node1, node2));
}
TEST(Node_Compare, compare4)
{
  ms::geo::Node node1({1, 2});
  ms::geo::Node node2({0.9, 2});

  ms::geo::Node_Compare compare;
  EXPECT_FALSE(compare(node1, node2));
}

TEST(Nodes_View, construct)
{
  const auto          num_nodes   = 3;
  const auto          dimension   = 2;
  std::vector<double> coordinates = {1, 2, 3, 3, 2, 1};

  ms::geo::Nodes_View nodes(num_nodes, dimension, coordinates);

  const auto& node1 = nodes[1];
  EXPECT_DOUBLE_EQ(node1[0], 3);
  EXPECT_DOUBLE_EQ(node1[1], 3);
}
TEST(Nodes_View, other_node_to_this_node)
{
  constexpr auto      dimension    = 3;
  std::vector<double> coordinates  = {1, 2, 3};
  std::vector<double> coordinates2 = {3, 2, 1};

  ms::geo::Node_View node1(coordinates);
  ms::geo::Node_View node2(coordinates2);

  ms::math::Vector<3> v;
  node1.other_to_this_vector(node2, v);

  ms::math::Vector<3> ref = {-2, 0, 2};
  EXPECT_EQ(v, ref);
}
TEST(Nodes_View, conversion1)
{
  constexpr auto      dimension   = 3;
  std::vector<double> coordinates = {1, 3, 5, 2, 4, 6};
  ms::geo::Nodes_View nodes(2, dimension, coordinates);

  ms::sym::Polynomial x("x0");
  ms::sym::Polynomial y("x1");
  ms::sym::Polynomial z("x2");

  auto p = x * x + y + y * z;

  const auto     val1 = p(nodes[0].to_vector_view());
  constexpr auto ref1 = 1 * 1 + 3 + 3 * 5;
  EXPECT_EQ(val1, ref1);

  const auto     val2 = p(nodes[1].to_vector_view());
  constexpr auto ref2 = 2 * 2 + 4 + 4 * 6;
  EXPECT_EQ(val2, ref2);
}
TEST(Nodes_View, axis_coordinates1)
{
  std::vector<double> coordinates = {1, 3, 5, 2, 4, 6, 1, 4, 3};
  constexpr auto      dimension   = 3;
  constexpr auto      num_nodes   = 3;
  ms::geo::Nodes_View nodes(num_nodes, dimension, coordinates);

  const auto result0 = nodes.axis_coordinates_vector_view(0);
  const auto result1 = nodes.axis_coordinates_vector_view(1);
  const auto result2 = nodes.axis_coordinates_vector_view(2);

  std::vector<double> ref0 = {1, 2, 1};
  std::vector<double> ref1 = {3, 4, 4};
  std::vector<double> ref2 = {5, 6, 3};
  EXPECT_EQ(result0, ref0);
  EXPECT_EQ(result1, ref1);
  EXPECT_EQ(result2, ref2);
}

TEST(Nodes, add_node1)
{
  ms::geo::Nodes nodes;

  ms::geo::Node node1 = {1, 2};
  ms::geo::Node node2({3, 4});
  ms::geo::Node node3({5, 6});

  nodes.add_node(node1);
  nodes.add_node(node2);
  nodes.add_node(node3);

  EXPECT_EQ(nodes[0], node1);
  EXPECT_EQ(nodes[1], node2);
  EXPECT_EQ(nodes[2], node3);
}
TEST(Nodes, add_node2)
{
  ms::geo::Nodes nodes(1, 2);
  ms::geo::Node  node1({1, 2});
  ms::geo::Node  node2({3, 4});
  ms::geo::Node  node3({5, 6});

  nodes.add_node(node1);
  nodes.add_node(node2);
  nodes.add_node(node3);

  EXPECT_EQ(nodes[0], ms::geo::Node(2));
  EXPECT_EQ(nodes[1], node1);
  EXPECT_EQ(nodes[2], node2);
  EXPECT_EQ(nodes[3], node3);
}
TEST(Nodes, copy_test1)
{
  constexpr auto num_nodes = 4;
  constexpr auto dimension = 2;
  ms::geo::Nodes nodes(num_nodes, dimension);

  for (int i = 0; i < num_nodes; ++i)
  {
    auto node_wrap = nodes[i];

    node_wrap[0] = i + 1;
    node_wrap[1] = i + 1;
  }

  auto other_nodes = nodes;

  auto node_wrap = nodes[2];
  node_wrap[0]   = -1;

  auto other_node_wrap = other_nodes[2];
  EXPECT_EQ(other_node_wrap[0], 3);
}

#ifdef _DEBUG
TEST(Nodes, move_test1)
{
  constexpr auto num_nodes = 4;
  constexpr auto dimension = 2;
  ms::geo::Nodes nodes(num_nodes, dimension);

  for (int i = 0; i < num_nodes; ++i)
  {
    auto node_wrap = nodes[i];

    node_wrap[0] = i + 1;
    node_wrap[1] = i + 1;
  }

  auto other_nodes = std::move(nodes);

  // EXPECT_ANY_THROW(nodes[2]); // access moved value is undefined behavior
}
#endif