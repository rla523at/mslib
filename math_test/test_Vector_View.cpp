#include "gtest\gtest.h"
#include "msmath/Vector.h"

namespace ms::math
{

// for google test cout message
std::ostream& operator<<(std::ostream& os, const Vector_View cvw)
{
  return os << cvw.to_string();
}

template <size_t dim>
std::ostream& operator<<(std::ostream& os, const Vector<dim>& vec)
{
  return os << vec.to_string();
}
// for google test cout message

TEST(Vector_Const_Iterator, is_input_iterator)
{
  EXPECT_TRUE(std::input_iterator<ms::math::Vector_Const_Iterator>);
}
TEST(Vector_Const_Iterator, contiguous_iterator)
{
  EXPECT_FALSE(std::contiguous_iterator<ms::math::Vector_Const_Iterator>);
}
TEST(Vector_View, const_span)
{
  EXPECT_TRUE(ms::math::const_span<std::span<double>>);
}

template <std::input_iterator T, ms::math::contiguous_iterator_or_int U>
void test(T t, U u)
{
  std::span<const double> s(t, u);
};

TEST(Vector_View, is_const_span_constructable)
{
  const std::vector<double> vec = {1, 2, 3};

  test(vec.begin(), vec.end());
}

TEST(Vector_View, construct_1)
{
  const std::vector<double> vec1 = {1, 2, 3};
  Vector_View               cvw1 = vec1;
  const auto                cvw2 = cvw1;

  EXPECT_EQ(cvw1, cvw2);
}
TEST(Vector_View, construct_2)
{
  const std::vector<double> vals = {1, 2, 3};
  Vector_View               vec  = vals;

  for (int i = 0; i < 3; ++i)
  {
    EXPECT_EQ(vals[i], vec[i]);
  }
}
TEST(Vector_View, construct_3)
{
  const std::vector<double> vals = {1, 2, 3};
  Vector_View               vec  = {vals, 2};

  const std::vector<double> ref = {1, 3};
  for (int i = 0; i < 2; ++i)
  {
    EXPECT_EQ(ref[i], vec[i]);
  }
}
TEST(Vector_View, construct_4)
{
  std::vector<double> v = {1, 2, 3, 4, 5, 6};
  Vector_View         vec(v.begin(), v.end() - 1);

  const std::vector<double> ref = {1, 2, 3, 4, 5};
  for (int i = 0; i < 5; ++i)
  {
    EXPECT_EQ(ref[i], vec[i]);
  }
}

TEST(Vector_View, iterator1)
{
  std::vector<double>   v = {1, 2, 3, 4, 5, 6};
  ms::math::Vector_View values(v, 2);

  auto                iter = values.begin();
  std::vector<double> ref  = {1, 3, 5};
  for (int i = 0; i < 3; ++i)
  {
    EXPECT_EQ(ref[i], *iter++);
  }
}
TEST(Vector_View, iterator2)
{
  std::vector<double>   v = {1, 2, 3, 4, 5, 6, 7};
  ms::math::Vector_View values(v, 3);

  auto                iter = values.begin();
  std::vector<double> ref  = {1, 4, 7};
  for (int i = 0; i < 3; ++i)
  {
    EXPECT_EQ(ref[i], *iter++);
  }
}
TEST(Vector_View, iterator3)
{
  std::vector<double>   v = {1, 2, 3, 4, 5, 6, 7};
  ms::math::Vector_View values1(v, 3);
  ms::math::Vector_View values2 = values1;

  auto iter1 = values1.begin();
  auto iter2 = values2.begin();

  std::vector<double> ref = {1, 4, 7};
  for (int i = 0; i < 3; ++i)
  {
    EXPECT_EQ(*iter1++, *iter2++);
  }
}

TEST(Vector_View, part_view1)
{
  std::vector<double>   v = {1, 2, 3, 4, 5, 6, 7};
  ms::math::Vector_View values(v);

  auto part = values.sub_view(1, 3);

  std::vector<double>   ref_v = {2, 3};
  ms::math::Vector_View ref_values(ref_v);
  EXPECT_EQ(ref_values, part);
}
TEST(Vector_View, part_view2)
{
  std::vector<double>   v = {1, 2, 3, 4, 5, 6, 7};
  ms::math::Vector_View values(v, 3);

  auto part = values.sub_view(1, 3);

  std::vector<double>   ref_v = {4, 7};
  ms::math::Vector_View ref_values(ref_v);
  EXPECT_EQ(ref_values, part);
}
TEST(Vector_View, part_view3)
{
  std::vector<double>   v = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  ms::math::Vector_View values(v);

  constexpr auto start_position = 0;
  constexpr auto inc            = 3;
  constexpr auto num_values     = 3;
  auto           sub_view      = values.sub_view(start_position, inc, num_values);

  std::vector<double> ref_v = {1, 4, 7};
  EXPECT_EQ(sub_view, ref_v);
}
TEST(Vector_View, part_view4)
{
  std::vector<double>   v = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  ms::math::Vector_View values(v, 3); // 1 4 7 10

  constexpr auto start_position = 0;
  constexpr auto inc            = 3;
  constexpr auto num_values     = 2;
  const auto     sub_view      = values.sub_view(start_position, inc, num_values);

  std::vector<double> ref_v = {1, 10};
  EXPECT_EQ(sub_view, ref_v);
}
TEST(Vector_View, part_view5)
{
  std::vector<double>   v = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
  ms::math::Vector_View values(v, 3); // 1 4 7 10 13

  constexpr auto start_position = 1;
  constexpr auto inc            = 3;
  constexpr auto num_values     = 2;
  const auto     sub_view      = values.sub_view(start_position, inc, num_values);

  std::vector<double> ref_v = {4, 13};
  EXPECT_EQ(sub_view, ref_v);
}

TEST(Vector_View, operator_addition_1)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {4, 5, 7};

  const Vector_View cvw1   = vec1;
  const auto        result = cvw1 + vec2;

  const Vector ref = {5, 7, 10};
  EXPECT_EQ(result, ref);
}
TEST(Vector_View, operator_addition_2)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {4, 5, 6};

  const Vector_View cvw1   = vec1;
  const Vector_View cvw2   = vec2;
  const auto        result = cvw1 + cvw2;

  const Vector ref = {5, 7, 9};
  EXPECT_EQ(result, ref);
}
TEST(Vector_View, operator_addition_3)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {4, 5, 6, 7, 8, 9};

  const Vector_View cvw1 = vec1;

  const auto        value_ptr = vec2.data();
  constexpr auto    num_value = 3;
  const Vector_View cvw2      = {value_ptr, num_value};

  const auto result = cvw1 + cvw2;

  const Vector ref = {5, 7, 9};
  EXPECT_EQ(result, ref);
}
TEST(Vector_View, operator_addition_4)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {4, 5, 6, 7, 8, 9};

  const Vector_View cvw1 = vec1;

  const auto        value_ptr = vec2.data() + 2;
  constexpr auto    num_value = 3;
  const Vector_View cvw2      = {value_ptr, num_value};

  const auto result = cvw1 + cvw2;

  const Vector ref = {7, 9, 11};
  EXPECT_EQ(result, ref);
}
TEST(Vector_View, operator_addition_5)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {4, 5, 6, 7, 8, 9};

  const Vector_View cvw1 = vec1;
  const Vector_View cvw2 = vec2;

  constexpr auto start_index  = 2;
  constexpr auto end_index    = 5;
  const auto     part_of_cvw2 = cvw2.sub_view(start_index, end_index);

  const auto result = cvw1 + part_of_cvw2;

  const Vector ref = {7, 9, 11};
  EXPECT_EQ(result, ref);
}
TEST(Vector_View, operator_addition_6)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {4, 5, 6, 7, 8, 9};

  const Vector_View cvw1 = vec1;

  constexpr auto    inc  = 2;
  const Vector_View cvw2 = {vec2, inc};

  const auto result = cvw1 + cvw2;

  const Vector ref = {5, 8, 11};
  EXPECT_EQ(result, ref);
}
TEST(Vector_View, operator_addition_7)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {4, 5, 6, 7, 8, 9};

  const Vector_View cvw1 = vec1;

  std::span<const double> s(vec2.begin() + 1, vec2.end());
  constexpr auto          inc  = 2;
  const Vector_View       cvw2 = {s, inc};

  const auto result = cvw1 + cvw2;

  const Vector ref = {6, 9, 12};
  EXPECT_EQ(result, ref);
}
TEST(Vector_View, operator_addition_8)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

  const Vector_View cvw1 = vec1;

  std::span<const double> s(vec2.begin(), 10);
  constexpr auto          inc  = 3;
  const Vector_View       cvw2 = {s, inc};

  constexpr auto start_pos = 0;
  constexpr auto end_pos   = 3;
  const auto     pcvw2     = cvw2.sub_view(start_pos, end_pos);

  const auto result = cvw1 + pcvw2;

  const Vector ref = {2, 6, 10};
  EXPECT_EQ(result, ref);
}
TEST(Vector_View, operator_scalar_multiplication_1)
{
  const std::vector<double> vec = {1, 2, 3};

  const Vector_View cvw1   = vec;
  const auto        result = 2 * cvw1;

  const Vector ref = {2, 4, 6};
  EXPECT_EQ(result, ref);
}
TEST(Vector_View, operator_scalar_multiplication_2)
{
  const std::vector<double> vec = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  std::span<const double> s      = {vec.begin() + 2, vec.end()};
  constexpr auto          inc    = 3;
  const Vector_View       cvw    = {s, inc};
  const auto              result = 2 * cvw;

  const Vector ref = {6, 12, 18};
  EXPECT_EQ(result, ref);
}
TEST(Vector_View, operator_equal_1)
{
  const std::vector<double> vec1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  const std::vector<double> vec2 = {3, 2, 1, 6, 5, 4, 9, 8, 7};

  std::span<const double> s1  = {vec1.begin() + 2, vec1.end()};
  std::span<const double> s2  = {vec2.begin(), vec2.end()};
  constexpr auto          inc = 3;

  Vector_View cvw1 = {s1, inc};
  Vector_View cvw2 = {s2, inc};

  EXPECT_TRUE(cvw1 == cvw2);
}
TEST(Vector_View, cosine_1)
{
  std::vector<double> values1 = {1, 2};
  std::vector<double> values2 = {-1, -2};

  Vector_View cvw1   = values1;
  Vector_View cvw2   = values2;
  const auto  result = cvw1.cosine(cvw2);

  const auto ref = -1.0;
  EXPECT_DOUBLE_EQ(ref, result);
}
TEST(Vector_View, cosine_2)
{
  std::vector<double> values1 = {1.123123, 2.92378479823};
  std::vector<double> values2 = {-1.123123, -2.92378479823};

  Vector_View cvw1   = values1;
  Vector_View cvw2   = values2;
  const auto  result = cvw1.cosine(cvw2);

  const auto ref = -1.0;
  EXPECT_DOUBLE_EQ(ref, result);
}
TEST(Vector_View, cross_product_1)
{
  std::vector<double> values = {1, 2};

  Vector_View cvw1   = values;
  Vector_View cvw2   = values;
  const auto  result = cvw1.cross_product<2>(cvw2);

  Vector ref = {0, 0, 0};
  EXPECT_EQ(ref, result);
}
TEST(Vector_View, cross_product_2)
{
  std::vector<double> values1 = {1, 2, 3};
  std::vector<double> values2 = {3, 2, 1};

  Vector_View cvw1   = values1;
  Vector_View cvw2   = values2;
  const auto  result = cvw1.cross_product<2>(cvw2);

  Vector ref = {0, 0, -4};
  EXPECT_EQ(ref, result);
}
TEST(Vector_View, cross_product_3)
{
  std::vector<double> values1 = {1, 2, 3};
  std::vector<double> values2 = {3, 2, 1};

  Vector_View cvw1   = values1;
  Vector_View cvw2   = values2;
  const auto  result = cvw1.cross_product<3>(cvw2);

  Vector ref = {-4, 8, -4};
  EXPECT_EQ(ref, result);
}
TEST(Vector_View, cross_product_4)
{
  std::vector<double> vec1 = {7, 8, 1, 2, 3, 2, 1, 1, 3};
  std::vector<double> vec2 = {1, 3, 2, 5, 2, 2, 3, 1, 1};

  constexpr auto inc  = 3;
  Vector_View    cvw1 = {vec1.begin() + 2, vec1.end(), inc};
  Vector_View    cvw2 = {vec2.begin() + 1, vec2.end() - 1, inc};

  const auto result = cvw1.cross_product<3>(cvw2);

  Vector ref = {-4, 8, -4};
  EXPECT_EQ(ref, result);
}
TEST(Vector_View, inner_product_1)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {4, 5, 6};

  const Vector_View cv1    = vec1;
  const auto        result = cv1.inner_product(vec2);

  const double ref = 32;
  EXPECT_EQ(result, ref);
}
TEST(Vector_View, inner_product_2)
{
  std::vector<double> vec = {3, 4};

  const Vector_View cv1    = vec;
  const auto        result = cv1.inner_product(cv1);

  const auto ref = 25;
  EXPECT_EQ(result, ref);
}
TEST(Vector_View, is_parallel_1)
{
  std::vector<double> values1 = {1, 2};
  std::vector<double> values2 = {-1, -2};

  Vector_View cvw1 = values1;
  Vector_View cvw2 = values2;
  EXPECT_TRUE(cvw1.is_parallel(cvw2));
}
TEST(Vector_View, is_parallel_2)
{
  std::vector<double> values1 = {1.123123, 2.92378479823};
  std::vector<double> values2 = {-1.123123, -2.92378479823};

  Vector_View cvw1 = values1;
  Vector_View cvw2 = values2;
  EXPECT_TRUE(cvw1.is_parallel(cvw2));
}
TEST(Vector_View, is_parallel_3)
{
  std::vector<double> values1 = {1.123123, 2.92378479823};

  Vector_View cvw1 = values1;
  auto        cvw2 = cvw1 * 2.1421323;
  EXPECT_TRUE(cvw1.is_parallel(cvw2));
}
TEST(Vector_View, L1_norm_1)
{
  std::vector<double> vec = {-3, 4};

  const Vector_View cv1    = vec;
  const auto        result = cv1.L1_norm();

  const auto ref = 7;
  EXPECT_EQ(result, ref);
}
TEST(Vector_View, L2_norm_1)
{
  std::vector<double> vec = {3, 4};

  const Vector_View cv1    = vec;
  const auto        result = cv1.L2_norm();

  const auto ref = 5;
  EXPECT_EQ(result, ref);
}

// Caution!
#ifdef _DEBUG
TEST(Vector_View, different_size_1)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {4, 5, 6, 7, 8};

  const Vector_View cvw1 = vec1;
  EXPECT_ANY_THROW(cvw1 + vec2);
}

#endif

TEST(Vector_View, changed_vector_1)
{
  std::vector<double> vec1 = {1, 2, 3};

  const Vector_View cvw1 = vec1;
  vec1                   = {4, 3, 2};

  EXPECT_EQ(cvw1, vec1);
}

TEST(Vector_View, moved_vector_1)
{
  std::vector<double> vec = {1, 2, 3};

  Vector_View v1 = vec;

  std::vector<double> new_vec = {4, 5, 6, 7, 8, 9};
  vec                         = std::move(new_vec); // move may invoke changing data ptr

  const std::vector<double> ref = {4, 5, 6, 7, 8, 9};
  EXPECT_NE(v1, ref);
}
} // namespace ms::math
