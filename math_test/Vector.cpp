#include "msmath/Vector.h"
#include "gtest\gtest.h"

namespace ms::math
{

// for google test cout message
template <size_t dim>
std::ostream& operator<<(std::ostream& os, const Vector<dim>& vec)
{
  return os << vec.to_string();
}
// for google test cout message

TEST(Vector, construct_1)
{
  std::vector<double> values = {1, 2, 3};
  Vector<3>           v      = {1, 2, 3};

  for (size_t i = 0; i < values.size(); i++)
  {
    EXPECT_EQ(v[i], values[i]);
  }
}
TEST(Vector, construct_2)
{
  std::vector<double> values = {1, 2, 3};
  Vector              v      = {1, 2, 3};

  for (size_t i = 0; i < values.size(); i++)
  {
    EXPECT_EQ(v[i], values[i]);
  }
}
TEST(Vector, construct_3)
{
  std::vector<double> values = {1, 2};
  Vector              v      = {1, 2};

  for (size_t i = 0; i < values.size(); i++)
  {
    EXPECT_EQ(v[i], values[i]);
  }
}
TEST(Vector, construct_4)
{
  std::vector<double> values = {1, 2, 3};
  Vector              v(values.begin(), values.end());

  for (size_t i = 0; i < values.size(); i++)
  {
    EXPECT_EQ(v[i], values[i]);
  }
}
TEST(Vector, construct_5)
{
  std::vector<double> values(5);
  Vector              v(5);

  for (size_t i = 0; i < values.size(); i++)
  {
    EXPECT_EQ(v[i], values[i]);
  }
}
TEST(Vector, construct_6)
{
  Vector<0> v = {1, 2, 3};

  std::vector<double> ref = {1, 2, 3};

  for (size_t i = 0; i < ref.size(); i++)
  {
    EXPECT_EQ(v[i], ref[i]);
  }
}
TEST(Vector, construct_7)
{
  Vector<0> v = {1, 2, 3};
  v           = {4, 5, 6, 7};

  std::vector<double> ref = {4, 5, 6, 7};

  for (size_t i = 0; i < ref.size(); i++)
  {
    EXPECT_EQ(v[i], ref[i]);
  }
}
TEST(Vector, operator_addition_1)
{
  const Vector v1     = {1, 2, 3};
  const Vector v2     = {4, 5, 6};
  const auto   result = v1 + v2;

  const Vector ref = {5, 7, 9};
  EXPECT_EQ(result, ref);
}
TEST(Vector, operator_addition_2)
{
  const Vector v1     = {1, 2, 3, 4};
  const Vector v2     = {4, 5, 6, 7};
  const auto   result = v1 + v2;

  const Vector ref = {5, 7, 9, 11};
  EXPECT_EQ(result, ref);
}
TEST(Vector, operator_addition_3)
{
  const Vector v1     = {1, 2, 3, 4};
  const Vector v2     = v1;
  const Vector v3     = {4, 5, 6, 7};
  const auto   result = v1 + v2 + v3;

  const Vector ref = {6, 9, 12, 15};
  EXPECT_EQ(result, ref);
}
TEST(Vector, operator_addition_4)
{
  Vector       v1 = {1, 2, 3, 4};
  const Vector v2 = v1;
  const Vector v3 = {4, 5, 6, 7};

  v1 *= 0;
  const auto result = v2 + v3;

  const Vector ref = {5, 7, 9, 11};
  EXPECT_EQ(result, ref);
}
TEST(Vector, operator_addition_5)
{
  Vector       v1 = {1, 2, 3, 4};
  const Vector v2 = std::move(v1);
  const Vector v3 = {4, 5, 6, 7};

  v1                = {2, 3, 4, 5};
  const auto result = v2 + v3;

  const Vector ref = {5, 7, 9, 11};
  EXPECT_EQ(result, ref);
}
TEST(Vector, operator_addition_6)
{
  Vector<3> v1 = {1, 2, 3};
  Vector<3> v2 = {1, 2};

  const auto result = v1 + v2;

  const Vector ref = {2, 4, 3};
  EXPECT_EQ(ref, result);
}
TEST(Vector, operator_addition_7)
{
  Vector<0> v1 = {1, 2, 3};
  Vector<3> v2 = {1, 2};

  const auto result = v1 + v2;

  const Vector ref = {2, 4, 3};
  EXPECT_EQ(ref, result);
}
TEST(Vector, operator_addition_assign_1)
{
  Vector v1 = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  v1 += v1;

  Vector ref = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
  EXPECT_EQ(v1, ref);
}
TEST(Vector, operator_addition_assign_2)
{
  Vector       v1 = {1, 2, 3, 4};
  const Vector v2 = {4, 5, 6, 7};
  v1 += v2;

  const Vector ref = {5, 7, 9, 11};
  EXPECT_EQ(v1, ref);
}
TEST(Vector, operator_addition_assign_3)
{
  Vector<0>    v1 = {1, 2, 3, 4};
  const Vector v2 = {4, 5, 6, 7};
  v1 += v2;

  const Vector ref = {5, 7, 9, 11};
  EXPECT_EQ(v1, ref);
}
TEST(Vector, operator_substraction_1)
{
  const Vector v1     = {1, 2, 3};
  const Vector v2     = {4, 5, 6};
  const auto   result = v1 - v2;

  const Vector ref = {-3, -3, -3};
  EXPECT_EQ(result, ref);
}
TEST(Vector, operator_scalar_multiplication_1)
{
  const Vector v1     = {1, 2, 3};
  const auto   result = v1 * 2;

  const Vector ref = {2, 4, 6};
  EXPECT_EQ(result, ref);
}
TEST(Vector, operator_scalar_multiplication_2)
{
  const Vector v1     = {1, 2, 3};
  const auto   result = 2 * v1;

  const Vector ref = {2, 4, 6};
  EXPECT_EQ(result, ref);
}
TEST(Vector, cross_product1)
{
  Vector v1 = {1, 2};
  Vector v2 = {1, 3};

  const auto result = v1.cross_product<2>(v2);

  const Vector ref = {0, 0, 1};
  EXPECT_EQ(ref, result);
}
TEST(Vector, cross_product2)
{
  Vector v1 = {1, 2, 3};
  Vector v2 = {1, 3, 3};

  const auto result = v1.cross_product<3>(v2);

  const Vector ref = {-3, 0, 1};
  EXPECT_EQ(ref, result);
}
TEST(Vector, cross_product3)
{
  Vector<0> v1 = {1, 2, 3};
  Vector    v2 = {1, 3, 3};

  const auto result = v1.cross_product<3>(v2);

  const Vector ref = {-3, 0, 1};
  EXPECT_EQ(ref, result);
}
TEST(Vector, inc_1)
{
  Vector     v(5);
  const auto result = v.inc();

  constexpr auto ref = 1;
  EXPECT_EQ(ref, result);
}
TEST(Vector, inner_product_1)
{
  const Vector v1     = {1, 2, 3};
  const Vector v2     = {4, 5, 6};
  const auto   result = v1.inner_product(v2);

  const double ref = 32;
  EXPECT_EQ(result, ref);
}
TEST(Vector, inner_product_2)
{
  Vector     v1     = {3, 4};
  const auto result = v1.inner_product(v1);

  const auto ref = 25;
  EXPECT_EQ(result, ref);
}
TEST(Vector, L2_norm_1)
{
  const Vector v1     = {3, 4};
  const auto   result = v1.L2_norm();

  const auto ref = 5;
  EXPECT_EQ(result, ref);
}
TEST(Vector, normalize_1)
{
  Vector v1 = {3, 4};
  v1.normalize();

  Vector ref = {3.0 / 5.0, 4.0 / 5.0};

  for (int i = 0; i < v1.dimension(); ++i)
    EXPECT_DOUBLE_EQ(v1[i], ref[i]);
}

// Caution!
#ifdef _DEBUG
TEST(Vector, different_size_1)
{
  Vector v1 = {0.1, 0.0, 0.6};
  Vector v2 = {0.1, 0.4};

  Vector ref = {0.2, 0.4, 0.6};
  EXPECT_ANY_THROW(v1 += v2);
}
#endif

} // namespace ms::math

//
// TEST(ms, theta_1)
//{
//  std::vector<double> values1 = { 1,0 };
//  std::vector<double> values2 = { 3,0 };
//
//  const auto result = ms::theta(values1, values2);
//
//  constexpr auto ref = 0;
//  EXPECT_EQ(ref, result);
//}
// TEST(ms, theta_2)
//{
//  std::vector<double> values1 = { 1,0 };
//  std::vector<double> values2 = { 0.1,0 };
//
//  const auto result = ms::theta(values1, values2);
//
//  constexpr auto ref = 0;
//  EXPECT_EQ(ref, result);
//}
// TEST(ms, theta_3)
//{
//  std::vector<double> values1 = { 1,0 };
//  std::vector<double> values2 = { -1,0 };
//
//  const auto result = ms::theta(values1, values2);
//
//  constexpr auto ref = 180;
//  EXPECT_EQ(ref, result);
//}
// TEST(ms, theta_4)
//{
//  std::vector<double> values1 = { 1,0 };
//  std::vector<double> values2 = { 0,-1 };
//
//  const auto result = ms::theta(values1, values2);
//
//  constexpr auto ref = 90;
//  EXPECT_EQ(ref, result);
//}
// TEST(ms, theta_5)
//{
//  std::vector<double> values1 = { 1,0,0 };
//  std::vector<double> values2 = { 1,1,std::sqrt(2) };
//
//  const auto result = ms::theta(values1, values2);
//
//  constexpr auto ref = 60;
//  EXPECT_DOUBLE_EQ(ref, result); //consider floating point epsilon
//}
// TEST(ms, theta_6)
//{
//  std::vector<double> values1 = { 1,0,0 };
//  std::vector<double> values2 = { 1,1,0 };
//
//  const auto result = ms::theta(values1, values2);
//
//  constexpr auto ref = 45;
//  EXPECT_DOUBLE_EQ(ref, result); //consider floating point epsilon
//}
// TEST(ms, singed_theta_1)
//{
//  std::vector<double> values1 = { 1,0 };
//  std::vector<double> values2 = { 0,-1 };
//
//  const auto result = ms::signed_theta<2>(values1, values2);
//
//  constexpr auto ref = -90;
//  EXPECT_EQ(ref, result);
//}
// TEST(ms, singed_theta_2)
//{
//  std::vector<double> values1 = { 1,0 };
//  std::vector<double> values2 = { 1,-1 };
//
//  const auto result = ms::signed_theta<2>(values1, values2);
//
//  constexpr auto ref = -45;
//  EXPECT_DOUBLE_EQ(ref, result); //consider floating point epsilon
//}
// TEST(ms, singed_theta_3)
//{
//  std::vector<double> values1 = { 1,0 };
//  std::vector<double> values2 = { -1,-1 };
//
//  const auto result = ms::signed_theta<2>(values1, values2);
//
//  constexpr auto ref = -135;
//  EXPECT_EQ(ref, result);
//}
// TEST(ms, singed_theta_4)
//{
//  std::vector<double> values1 = { 1,0 };
//  std::vector<double> values2 = { 0,1 };
//
//  const auto result = ms::signed_theta<2>(values1, values2);
//
//  constexpr auto ref = 90;
//  EXPECT_EQ(ref, result);
//}

////TEST(Dynamic_Real_Vector, be_absolute_1)
////{
////	Dynamic_Real_Vector v1 = { -1, -2 ,3 ,-4 ,5 };
////	v1.be_absolute();
////
////	const Dynamic_Real_Vector ref = { 1,2,3,4,5 };
////	EXPECT_EQ(v1, ref);
////}
////
////TEST(ms, min_value_gathering_vector_1)
////{
////	Real_Vector v1 = { 1,2,3 };
////	Real_Vector v2 = { 2,3,1 };
////	Real_Vector v3 = { 3,1,2 };
////
////	std::vector<Real_Vector<3>> vec = { v1,v2,v3 };
////
////	const auto result = ms::gather_min_value(vec);
////	const Real_Vector ref = { 1,1,1 };
////	EXPECT_EQ(result, ref);
////}
////
////TEST(ms, max_value_gathering_vector_1)
////{
////	Real_Vector v1 = { 1,2,3 };
////	Real_Vector v2 = { 2,3,1 };
////	Real_Vector v3 = { 3,2,1 };
////
////	std::vector<Real_Vector<3>> vec = { v1,v2,v3 };
////
////	const auto result = ms::gather_max_value(vec);
////	const Real_Vector ref = { 3,3,3 };
////	EXPECT_EQ(result, ref);
////}

//
// TEST(Vector_Function, constructor)
//{
//	constexpr ushort domain_dimension = 2;
//
//	Polynomial x("x0");
//	Polynomial y("x1");
//
//	constexpr ushort range_dimension = 3;
//	Vector_Function<Polynomial> vf = { 1,x,y };
//
//	EXPECT_EQ(vf.at(0), 1);
//	EXPECT_EQ(vf.at(1), x);
//	EXPECT_EQ(vf.at(2), y);
//}
// TEST(Vector_Function, at_1)
//{
//	Polynomial x("x0");
//	Polynomial y("x1");
//
//	Vector_Function<Polynomial> vf = { 0 , 2 * x + y };
//
//	const auto ref = 2 * x + y;
//	EXPECT_EQ(ref, vf.at(1));
//}
// TEST(Vector_Function, operator_call_1)
//{
//	Polynomial x("x0");
//	Polynomial y("x1");
//
//	Vector_Function<Polynomial> vf = { x + y , 2 * x + y };
//	std::vector<double> v = { 1,1 };
//	const auto result = vf(v);
//
//	std::vector<double> ref = { 2,3 };
//	EXPECT_EQ(ref, result);
//}
// TEST(Vector_Function, operator_call_2)
//{
//	Polynomial x("x0");
//	Polynomial y("x1");
//	Polynomial z("x2");
//
//	auto p1 = (x ^ 2) + 3 * (x ^ 2) * y + (y ^ 3) + (z ^ 2) - 6;
//	auto p2 = x + y + z - 3;
//	auto p3 = (y ^ 2) * z + x * z - 2;
//
//	Vector_Function<Polynomial> f = { p1,p2,p3 };
//	std::vector<double> node = { 1,1,1 };
//	const auto result = f(node);
//
//	std::vector<double> ref = { 0,0,0 };
//	EXPECT_EQ(result, ref);
//}
// TEST(Vector_Function, get_differentiate_1)
//{
//	Polynomial x("x0");
//	Polynomial y("x1");
//	Polynomial z("x2");
//
//	auto p1 = x + 1;
//	auto p2 = x + y + z - 3;
//	auto p3 = y + z;
//	Vector_Function<Polynomial> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 0;
//	const auto result = f.get_differentiate(variable_index);
//
//	Vector_Function<Polynomial> ref = { 1,1,0 };
//	EXPECT_EQ(result, ref);
//}
// TEST(Vector_Function, get_differentiate_2)
//{
//	Polynomial x("x0");
//	Polynomial y("x1");
//	Polynomial z("x2");
//
//	auto p1 = x * y + 1;
//	auto p2 = x + y * z + z - 3;
//	auto p3 = y + z;
//	Vector_Function<Polynomial> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 0;
//	const auto result = f.get_differentiate(variable_index);
//
//	Vector_Function<Polynomial> ref = { y,1,0 };
//	EXPECT_EQ(result, ref);
//}
// TEST(Vector_Function, get_differentiate_3)
//{
//	Polynomial x("x0");
//	Polynomial y("x1");
//	Polynomial z("x2");
//
//	auto p1 = x * y + 1;
//	auto p2 = x + y * z + z - 3;
//	auto p3 = y + z;
//	Vector_Function<Polynomial> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 1;
//	const auto result = f.get_differentiate(variable_index);
//
//	Vector_Function<Polynomial> ref = { x,z,1 };
//	EXPECT_EQ(result, ref);
//}
// TEST(Vector_Function, get_differentiate_4)
//{
//	Polynomial x("x0");
//	Polynomial y("x1");
//	Polynomial z("x2");
//
//	auto p1 = 1.5 * x + 0.5 * y + 3;
//	auto p2 = y + 3;
//	auto p3 = 0;
//	Vector_Function<Polynomial> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 0;
//	const auto result = f.get_differentiate(variable_index);
//
//	Vector_Function<Polynomial> ref = { 1.5,0,0 };
//	EXPECT_EQ(result, ref);
//}
// TEST(Vector_Function, get_differentiate_5)
//{
//	Polynomial x("x0");
//	Polynomial y("x1");
//	Polynomial z("x2");
//
//	auto p1 = 1.5 * x + 0.5 * y + 3;
//	auto p2 = y + 3;
//	auto p3 = 0;
//	Vector_Function<Polynomial> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 1;
//	const auto result = f.get_differentiate(variable_index);
//
//	Vector_Function<Polynomial> ref = { 0.5,1,0 };
//	EXPECT_EQ(result, ref);
//}
// TEST(Vector_Function, get_differentiate_6)
//{
//	Polynomial x("x0");
//	Polynomial y("x1");
//
//	Vector_Function<Polynomial> vf = { x + y , 2 * x + y };
//	const auto result = vf.get_differentiate(0);
//
//	Vector_Function<Polynomial> ref = { 1,2 };
//	EXPECT_EQ(ref, result);
//}
// TEST(Vector_Function, get_differentiate_7)
//{
//	Polynomial x("x0");
//	Polynomial y("x1");
//
//	Vector_Function<Polynomial> vf = { 0.25 * x * y + 1.25 * x + 0.25 * y + 2.25, -0.25 * x * y - 0.75 * x + 0.25 * y + 1.75 };
//	const auto result = vf.get_differentiate(0);
//
//	Vector_Function<Polynomial> ref = { 0.25 * y + 1.25,-0.25 * y - 0.75 };
//	EXPECT_EQ(ref, result);
//}
// TEST(Vector_Function, cross_product_1)
//{
//	Vector_Function<Polynomial> vf1 = { 1.5,0,0 };
//	Vector_Function<Polynomial> vf2 = { 0.5,1,0 };
//	const auto result = vf1.cross_product(vf2);
//
//	Vector_Function<Polynomial> ref = { 0,0,1.5 };
//	EXPECT_EQ(result, ref);
//}
// TEST(Vector_Function, cross_product_2)
//{
//	Polynomial x("x0");
//	Polynomial y("x1");
//
//	Vector_Function<Polynomial> vf1 = { 1 + x ,		2 * x + y,		3 };
//	Vector_Function<Polynomial> vf2 = { 2 * x - y ,	1 + (-1 * x),	2 };
//	const auto result = vf1.cross_product(vf2);
//
//	Vector_Function<Polynomial> ref = { 7 * x + 2 * y - 3,4 * x - 3 * y - 2, 1 - 5 * (x ^ 2) + (y ^ 2) };
//	EXPECT_EQ(ref, result);
//}
// TEST(Vector_Function, cross_product_3)
//{
//	Polynomial x("x0");
//	Polynomial y("x1");
//
//	Vector_Function<Polynomial> vf1 = { 1 + x ,		2 * x + y };
//	Vector_Function<Polynomial> vf2 = { 2 * x - y ,	1 + (-1 * x) };
//	const auto result = vf1.cross_product(vf2);
//
//	Vector_Function<Polynomial> ref = { 0, 0, 1 - 5 * (x ^ 2) + (y ^ 2) };
//	EXPECT_EQ(ref, result);
//}
////TEST(Vector_Function, L2_NORM1)
////{
////	Vector_Function<Polynomial> vf = { 0,0,1.5 };
////	const auto result = vf.L2_norm();
////
////	Polynomial ref = 1.5;
////	EXPECT_EQ(result, ref);
////}
//// TEST(Vector_Function, L2_norm_1) {
////constexpr size_t domain_dimension = 2;
////Polynomial x("x0");
////Polynomial y("x1");
////
////constexpr ushort range_dimension = 3;
////
////Vector_Function<Polynomial, range_dimension> vf1 = { 1, x, y };
////const auto result = vf1.L2_norm();
////
////Irrational_Function ref(1 + (x ^ 2) + (y ^ 2), 0.5);
////EXPECT_EQ(ref, result);
////}
////TEST(ms, scalar_triple_product_1) {
////	constexpr ushort domain_dimension = 3;
////
////	Polynomial x("x0");
////	Polynomial y("x1");
////	Polynomial z("x1");
////
////	Vector_Function<Polynomial, domain_dimension> a = { x + y, x, y };
////	Vector_Function<Polynomial, domain_dimension> b = { x, x, 0 };
////	Vector_Function<Polynomial, domain_dimension> c = { -1 * y, y, 1 };
////	const auto result = ms::scalar_triple_product(a, b, c);
////
////	Polynomial ref = x * y + 2 * x * (y ^ 2);
////	EXPECT_EQ(result, ref);
////}