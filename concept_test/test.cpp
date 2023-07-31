#include "msconcept/concept.h"
#include "gtest/gtest.h"
#include <array>

template <typename... Args>
bool g(const Args...)
{
  return false;
}

template <ms::const_span_constructable T>
bool g(const T& t)
{
  std::span<const double> s(t);
  return true;
}

template <typename T, typename U>
requires ms::is_const_span_constructable<T, U> 
bool g(const T& t, const U& u)
{
  std::span<const double> s(t, u);
  return true;
}

TEST(concept, test5)
{
  std::vector<double>   v;
  std::array<double, 3> ar = {0};
  double                rar[3];

  auto& v_ref   = v;
  auto& ar_ref  = ar;
  auto& rar_ref = rar;

  const auto& v_cref   = v;
  const auto& ar_cref  = ar;
  const auto& rar_cref = rar;

  EXPECT_TRUE(g(v));
  EXPECT_TRUE(g(ar));
  EXPECT_TRUE(g(rar));
  EXPECT_TRUE(g(v.begin(), v.end()));
  EXPECT_TRUE(g(ar.begin(), ar.end()));

  EXPECT_TRUE(g(v_ref));
  EXPECT_TRUE(g(ar_ref));
  EXPECT_TRUE(g(rar_ref));
  EXPECT_TRUE(g(v_ref.begin(), v_ref.end()));
  EXPECT_TRUE(g(ar_ref.begin(), ar_ref.end()));

  EXPECT_TRUE(g(v_cref));
  EXPECT_TRUE(g(ar_cref));
  EXPECT_TRUE(g(rar_cref));
  EXPECT_TRUE(g(v_cref.begin(), v_cref.end()));
  EXPECT_TRUE(g(ar_cref.begin(), ar_cref.end()));
}