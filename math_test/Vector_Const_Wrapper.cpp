#include "gtest\gtest.h"
#include "msmath/Vector.h"

namespace ms::math
{

// for google test cout message
std::ostream& operator<<(std::ostream& os, const Vector_Const_Wrapper& cvw)
{
  return os << cvw.to_string();
}

template <size_t dim>
std::ostream& operator<<(std::ostream& os, const Vector<dim>& vec)
{
  return os << vec.to_string();
}
// for google test cout message

TEST(Vector_Const_Wrapper, construct_1)
{
  const std::vector<double>  vec1 = {1, 2, 3};
  const Vector_Const_Wrapper cvw1 = vec1;
  const auto                 cvw2 = cvw1;

  EXPECT_EQ(cvw1, cvw2);
}

TEST(Vector_Const_Wrapper, operator_addition_1)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {4, 5, 7};

  const Vector_Const_Wrapper cvw1   = vec1;
  const auto                 result = cvw1 + vec2;

  const Vector ref = {5, 7, 10};
  EXPECT_EQ(result, ref);
}
TEST(Vector_Const_Wrapper, operator_addition_2)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {4, 5, 6};

  const Vector_Const_Wrapper cvw1   = vec1;
  const Vector_Const_Wrapper cvw2   = vec2;
  const auto                 result = cvw1 + cvw2;

  const Vector ref = {5, 7, 9};
  EXPECT_EQ(result, ref);
}
TEST(Vector_Const_Wrapper, operator_addition_3)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {4, 5, 6, 7, 8, 9};

  const Vector_Const_Wrapper cvw1 = vec1;

  const auto                 value_ptr = vec2.data();
  constexpr auto             num_value = 3;
  const Vector_Const_Wrapper cvw2      = {value_ptr, num_value};

  const auto result = cvw1 + cvw2;

  const Vector ref = {5, 7, 9};
  EXPECT_EQ(result, ref);
}
TEST(Vector_Const_Wrapper, operator_addition_4)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {4, 5, 6, 7, 8, 9};

  const Vector_Const_Wrapper cvw1 = vec1;

  const auto                 value_ptr = vec2.data() + 2;
  constexpr auto             num_value = 3;
  const Vector_Const_Wrapper cvw2      = {value_ptr, num_value};

  const auto result = cvw1 + cvw2;

  const Vector ref = {7, 9, 11};
  EXPECT_EQ(result, ref);
}
TEST(Vector_Const_Wrapper, operator_addition_5)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {4, 5, 6, 7, 8, 9};

  const Vector_Const_Wrapper cvw1 = vec1;
  const Vector_Const_Wrapper cvw2 = vec2;

  constexpr auto start_index  = 2;
  constexpr auto end_index    = 5;
  const auto     part_of_cvw2 = cvw2.part(start_index, end_index);

  const auto result = cvw1 + part_of_cvw2;

  const Vector ref = {7, 9, 11};
  EXPECT_EQ(result, ref);
}
TEST(Vector_Const_Wrapper, operator_addition_6)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {4, 5, 6, 7, 8, 9};

  const Vector_Const_Wrapper cvw1 = vec1;

  constexpr auto             inc  = 2;
  const Vector_Const_Wrapper cvw2 = {vec2, inc};

  const auto result = cvw1 + cvw2;

  const Vector ref = {5, 8, 11};
  EXPECT_EQ(result, ref);
}
TEST(Vector_Const_Wrapper, operator_addition_7)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {4, 5, 6, 7, 8, 9};

  const Vector_Const_Wrapper cvw1 = vec1;

  const auto                 start_ptr = vec2.data() + 1;
  constexpr auto             dim       = 3;
  constexpr auto             inc       = 2;
  const Vector_Const_Wrapper cvw2      = {start_ptr, dim, inc};

  const auto result = cvw1 + cvw2;

  const Vector ref = {6, 9, 12};
  EXPECT_EQ(result, ref);
}
TEST(Vector_Const_Wrapper, operator_addition_8)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

  const Vector_Const_Wrapper cvw1 = vec1;

  const auto                 start_ptr = vec2.data();
  constexpr auto             dim       = 4;
  constexpr auto             inc       = 3;
  const Vector_Const_Wrapper cvw2      = {start_ptr, dim, inc};

  constexpr auto start_pos = 0;
  constexpr auto end_pos   = 3;
  const auto     pcvw2     = cvw2.part(start_pos, end_pos);

  const auto result = cvw1 + pcvw2;

  const Vector ref = {2, 6, 10};
  EXPECT_EQ(result, ref);
}
TEST(Vector_Const_Wrapper, operator_scalar_multiplication_1)
{
  const std::vector<double> vec = {1, 2, 3};

  const Vector_Const_Wrapper cvw1   = vec;
  const auto                 result = 2 * cvw1;

  const Vector ref = {2, 4, 6};
  EXPECT_EQ(result, ref);
}
TEST(Vector_Const_Wrapper, operator_scalar_multiplication_2)
{
  const std::vector<double> vec = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  const auto     start_ptr = vec.data() + 2;
  constexpr auto n         = 3;
  constexpr auto inc       = 3;

  const Vector_Const_Wrapper cvw    = {start_ptr, n, inc};
  const auto                 result = 2 * cvw;

  const Vector ref = {6, 12, 18};
  EXPECT_EQ(result, ref);
}
TEST(Vector_Const_Wrapper, operator_equal_1)
{
  const std::vector<double> vec1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  const std::vector<double> vec2 = {3, 2, 1, 6, 5, 4, 9, 8, 7};

  const auto     start_ptr1 = vec1.data() + 2;
  const auto     start_ptr2 = vec2.data();
  constexpr auto n          = 3;
  constexpr auto inc        = 3;

  Vector_Const_Wrapper cvw1 = {start_ptr1, n, inc};
  Vector_Const_Wrapper cvw2 = {start_ptr2, n, inc};

  EXPECT_TRUE(cvw1 == cvw2);
}
TEST(Vector_Const_Wrapper, cosine_1)
{
  std::vector<double> values1 = {1, 2};
  std::vector<double> values2 = {-1, -2};

  Vector_Const_Wrapper cvw1   = values1;
  Vector_Const_Wrapper cvw2   = values2;
  const auto           result = cvw1.cosine(cvw2);

  const auto ref = -1.0;
  EXPECT_DOUBLE_EQ(ref, result);
}
TEST(Vector_Const_Wrapper, cosine_2)
{
  std::vector<double> values1 = {1.123123, 2.92378479823};
  std::vector<double> values2 = {-1.123123, -2.92378479823};

  Vector_Const_Wrapper cvw1   = values1;
  Vector_Const_Wrapper cvw2   = values2;
  const auto           result = cvw1.cosine(cvw2);

  const auto ref = -1.0;
  EXPECT_DOUBLE_EQ(ref, result);
}
TEST(Vector_Const_Wrapper, cross_product_1)
{
  std::vector<double> values = {1, 2};

  Vector_Const_Wrapper cvw1   = values;
  Vector_Const_Wrapper cvw2   = values;
  const auto           result = cvw1.cross_product<2>(cvw2);

  Vector ref = {0, 0, 0};
  EXPECT_EQ(ref, result);
}
TEST(Vector_Const_Wrapper, cross_product_2)
{
  std::vector<double> values1 = {1, 2, 3};
  std::vector<double> values2 = {3, 2, 1};

  Vector_Const_Wrapper cvw1   = values1;
  Vector_Const_Wrapper cvw2   = values2;
  const auto           result = cvw1.cross_product<2>(cvw2);

  Vector ref = {0, 0, -4};
  EXPECT_EQ(ref, result);
}
TEST(Vector_Const_Wrapper, cross_product_3)
{
  std::vector<double> values1 = {1, 2, 3};
  std::vector<double> values2 = {3, 2, 1};

  Vector_Const_Wrapper cvw1   = values1;
  Vector_Const_Wrapper cvw2   = values2;
  const auto           result = cvw1.cross_product<3>(cvw2);

  Vector ref = {-4, 8, -4};
  EXPECT_EQ(ref, result);
}
TEST(Vector_Const_Wrapper, cross_product_4)
{
  std::vector<double> values1 = {7, 8, 1, 2, 3, 2, 1, 1, 3};
  std::vector<double> values2 = {1, 3, 2, 5, 2, 2, 3, 1, 1};

  const auto           start_data_ptr1 = values1.data() + 2;
  constexpr auto       n               = 3;
  constexpr auto       inc             = 3;
  Vector_Const_Wrapper cvw1            = {start_data_ptr1, n, inc};

  const auto           start_data_ptr2 = values2.data() + 1;
  Vector_Const_Wrapper cvw2            = {start_data_ptr2, n, inc};
  ;

  const auto result = cvw1.cross_product<3>(cvw2);

  Vector ref = {-4, 8, -4};
  EXPECT_EQ(ref, result);
}
TEST(Vector_Const_Wrapper, inner_product_1)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {4, 5, 6};

  const Vector_Const_Wrapper cv1    = vec1;
  const auto                 result = cv1.inner_product(vec2);

  const double ref = 32;
  EXPECT_EQ(result, ref);
}
TEST(Vector_Const_Wrapper, inner_product_2)
{
  std::vector<double> vec = {3, 4};

  const Vector_Const_Wrapper cv1    = vec;
  const auto                 result = cv1.inner_product(cv1);

  const auto ref = 25;
  EXPECT_EQ(result, ref);
}
TEST(Vector_Const_Wrapper, is_parallel_1)
{
  std::vector<double> values1 = {1, 2};
  std::vector<double> values2 = {-1, -2};

  Vector_Const_Wrapper cvw1 = values1;
  Vector_Const_Wrapper cvw2 = values2;
  EXPECT_TRUE(cvw1.is_parallel(cvw2));
}
TEST(Vector_Const_Wrapper, is_parallel_2)
{
  std::vector<double> values1 = {1.123123, 2.92378479823};
  std::vector<double> values2 = {-1.123123, -2.92378479823};

  Vector_Const_Wrapper cvw1 = values1;
  Vector_Const_Wrapper cvw2 = values2;
  EXPECT_TRUE(cvw1.is_parallel(cvw2));
}
TEST(Vector_Const_Wrapper, is_parallel_3)
{
  std::vector<double> values1 = {1.123123, 2.92378479823};

  Vector_Const_Wrapper cvw1 = values1;
  auto                 cvw2 = cvw1 * 2.1421323;
  EXPECT_TRUE(cvw1.is_parallel(cvw2));
}
TEST(Vector_Const_Wrapper, L1_norm_1)
{
  std::vector<double> vec = {-3, 4};

  const Vector_Const_Wrapper cv1    = vec;
  const auto                 result = cv1.L1_norm();

  const auto ref = 7;
  EXPECT_EQ(result, ref);
}
TEST(Vector_Const_Wrapper, L2_norm_1)
{
  std::vector<double> vec = {3, 4};

  const Vector_Const_Wrapper cv1    = vec;
  const auto                 result = cv1.L2_norm();

  const auto ref = 5;
  EXPECT_EQ(result, ref);
}

// Caution!
#ifdef _DEBUG
TEST(Vector_Const_Wrapper, different_size_1)
{
  const std::vector<double> vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {4, 5, 6, 7, 8};

  const Vector_Const_Wrapper cvw1 = vec1;
  EXPECT_ANY_THROW(cvw1 + vec2);
}
#endif

TEST(Vector_Const_Wrapper, changed_vector_1)
{
  std::vector<double> vec1 = {1, 2, 3};

  const Vector_Const_Wrapper cvw1 = vec1;
  vec1                            = {4, 3, 2};

  EXPECT_EQ(cvw1, vec1);
}

TEST(Vector_Const_Wrapper, moved_vector_1)
{
  std::vector<double> vec = {1, 2, 3};

  Vector_Const_Wrapper v1 = vec;

  std::vector<double> new_vec = {4, 5, 6, 7, 8, 9};
  vec                         = std::move(new_vec); // move may invoke changing data ptr

  const std::vector<double> ref = {4, 5, 6, 7, 8, 9};
  EXPECT_NE(v1, ref);
}
} // namespace ms::math
