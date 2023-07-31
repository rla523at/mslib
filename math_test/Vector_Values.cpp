#include "gtest\gtest.h"
#include "msmath/Vector.h"

TEST(Vector_Values_Const_Iterator, is_input_iterator)
{
  EXPECT_TRUE(std::input_iterator<ms::math::Vector_Values_Const_Iterator>);
}

TEST(Vector_Values_Const_Wrapper, const_span)
{
  EXPECT_TRUE(ms::math::const_span<std::span<double>>);
}

template <std::input_iterator T, ms::math::input_iterator_or_int U>
void test(T t, U u)
{
  std::span<const double> s(t, u);
};

TEST(Vector_Values_Const_Wrapper, is_const_span_constructable)
{
  const std::vector<double> vec = {1, 2, 3};

  test(vec.begin(), vec.end());
}

TEST(Vector_Values_Const_Wrapper, construct)
{
  std::vector<double>                   v = {1, 2, 3, 4, 5, 6};
  ms::math::Vector_Values_Const_Wrapper values(v, 5);
}
TEST(Vector_Values_Const_Wrapper, iterator1)
{
  std::vector<double>                   v = {1, 2, 3, 4, 5, 6};
  ms::math::Vector_Values_Const_Wrapper values(v, 2);

  auto                iter = values.begin();
  std::vector<double> ref  = {1, 3, 5};
  for (int i = 0; i < 3; ++i)
  {
    EXPECT_EQ(ref[i], *iter++);
  }
}
TEST(Vector_Values_Const_Wrapper, iterator2)
{
  std::vector<double>                   v = {1, 2, 3, 4, 5, 6, 7};
  ms::math::Vector_Values_Const_Wrapper values(v, 3);

  auto                iter = values.begin();
  std::vector<double> ref  = {1, 4, 7};
  for (int i = 0; i < 3; ++i)
  {
    EXPECT_EQ(ref[i], *iter++);
  }
}

TEST(Vector_Values_Const_Wrapper, part1)
{
  std::vector<double>                   v = {1, 2, 3, 4, 5, 6, 7};
  ms::math::Vector_Values_Const_Wrapper values(v);

  auto part = values.part(1, 3);

  std::vector<double>                   ref_v = {2, 3};
  ms::math::Vector_Values_Const_Wrapper ref_values(ref_v);
  EXPECT_EQ(ref_values, part);
}
TEST(Vector_Values_Const_Wrapper, part2)
{
  std::vector<double>                   v = {1, 2, 3, 4, 5, 6, 7};
  ms::math::Vector_Values_Const_Wrapper values(v, 3);

  auto part = values.part(1, 3);

  std::vector<double>                   ref_v = {4, 7};
  ms::math::Vector_Values_Const_Wrapper ref_values(ref_v);
  EXPECT_EQ(ref_values, part);
}

TEST(Vector_Values_Wrapper, construct1)
{
  std::vector<double> vec = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  constexpr auto                  inc = 3;
  ms::math::Vector_Values_Wrapper v   = {vec, inc};

  std::vector<double>             ref_v = {1, 4, 7};
  ms::math::Vector_Values_Wrapper ref   = ref_v;

  EXPECT_EQ(v, ref);
}