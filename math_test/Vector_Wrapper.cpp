#include "gtest\gtest.h"
#include "msmath/Vector.h"

namespace ms::math
{

TEST(Vector_Wrap, operator_assign_1)
{
  std::vector<double>       vec1 = {0, 0, 0};
  const std::vector<double> vec2 = {1, 2, 3};
  const std::vector<double> vec3 = {6, 5, 4};

  Vector_Wrap       vw1  = vec1;
  const Vector_View cvw2 = vec2;

  vw1.change_value(cvw2 + vec3);

  const std::vector<double> ref = {7, 7, 7};
  EXPECT_EQ(vec1, ref);
}
TEST(Vector_Wrap, operator_addition_assign_1)
{
  std::vector<double>       vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {6, 5, 4};

  Vector_Wrap vw1 = vec1;
  vw1 += vec2;

  const std::vector<double> ref = {7, 7, 7};
  EXPECT_EQ(vec1, ref);
}
TEST(Vector_Wrap, operator_addition_assign_2)
{
  std::vector<double>       vec1 = {1, 2, 3};
  const std::vector<double> vec2 = {6, 5, 4, 3, 2, 1};

  Vector_Wrap v1 = vec1;
  Vector_View cv2(vec2.data(), 3);
  v1 += cv2;

  const Vector ref = {7, 7, 7};
  EXPECT_EQ(v1, ref);
}
TEST(Vector_Wrap, operator_scalar_multiplication_assign_1)
{
  std::vector<double> vec = {1, 2, 3};

  Vector_Wrap v = vec;
  v *= 2;

  const std::vector<double> ref = {2, 4, 6};
  EXPECT_EQ(vec, ref);
}
TEST(Vector_Wrap, operator_scalar_multiplication_assign_2)
{
  std::vector<double> vec = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  constexpr auto      inc = 3;
  Vector_Wrap         v   = {vec, inc};
  v *= 2;

  const std::vector<double> ref = {2, 2, 3, 8, 5, 6, 14, 8, 9};
  EXPECT_EQ(vec, ref);
}

// Caution!
#ifdef _DEBUG
TEST(Vector_Wrap, different_size_1)
{
  std::vector<double>       vec1 = {1, 2, 3, 4};
  const std::vector<double> vec2 = {1, 2};

  Vector_Wrap vw1 = vec1;
  EXPECT_ANY_THROW(vw1 += vec2);
}
#endif

} // namespace ms::math