#pragma once
#include "mssym/Symbol.h"
#include "mssym/Polynomial.h"
#include "gtest/gtest.h"

namespace ms::sym
{

Symbol cal_curvature(const Polynomials& parametric_curve);

TEST(Symbol, Complicate_Operation1)
{
  Polynomial x("x0");
  auto       p1 = x + 1;
  auto       p2 = x + 3;
  auto       p3 = x * x + 2;
  auto       p4 = x + 4;

  Symbol sym1 = p1.to_sym_base();
  sym1.pow(-1.5);

  Symbol sym2 = p2.to_sym_base();

  Symbol sym3 = p3.to_sym_base();
  sym3.pow(-0.5);

  Symbol sym4 = p4.to_sym_base();

  auto sym = sym1 * sym2 + sym3 * sym4;

  double     d      = 1.0;
  auto       result = sym(&d);
  const auto ref    = std::sqrt(2.0) + (5 * std::sqrt(3) / 3);
  EXPECT_DOUBLE_EQ(result, ref);
}
TEST(Symbol, Complicate_Operation2)
{
  Polynomial x("x0");

  Polynomials pfs(3);
  pfs[0] = 1.5 * x + 2.5;
  pfs[1] = 2;
  pfs[2] = 2 * x + 5;

  const auto     curvature = cal_curvature(pfs);
  const double   d         = 1.0;
  const auto     result    = curvature(&d);
  constexpr auto ref       = 0.0;
  EXPECT_DOUBLE_EQ(result, ref);
}

TEST(Symbol, divide)
{
  Polynomial x("x0");
  auto       p1 = x * x + 2 * x + 3;
  auto       p2 = 2 * x + 3;

  Symbol sym1 = p1.to_sym_base();
  Symbol sym2 = p2.to_sym_base();

  auto sym = sym1 / sym2;

  double         d      = 1.0;
  const auto     result = sym(&d);
  constexpr auto ref    = 1.2;
  EXPECT_DOUBLE_EQ(result, ref);
}

TEST(Symbol, root)
{
  Polynomial x("x0");
  auto       p1 = x * x + 2 * x + 6;

  Symbol sym = p1.to_sym_base();
  sym.pow(0.5);

  double         d      = 1.0;
  const auto     result = sym(&d);
  constexpr auto ref    = 3;
  EXPECT_DOUBLE_EQ(result, ref);
}

TEST(Symbol, differentiate)
{
  Polynomial x("x0");
  auto       p1 = x * x + 2 * x + 3;
  auto       p2 = 2 * x + 3;

  Symbol sym1 = p1.to_sym_base();
  Symbol sym2 = p2.to_sym_base();
  auto   sym  = sym1 / sym2;

  constexpr auto var_index = 0;
  auto           diff_sym  = sym.get_diff_symbol(var_index);

  double         d      = 1.0;
  const auto     result = sym(&d);
  constexpr auto ref    = 1.2;
  EXPECT_DOUBLE_EQ(result, ref);
}

TEST(Symbol, absolute)
{
  const auto x = Polynomial("x0");
  const auto y = Polynomial("x1");
  const auto z = Polynomial("x2");

  auto   p   = x + y + z;
  Symbol sym = p.to_sym_base();
  sym.be_absolute();

  const std::vector<double> values = {4, -7, -8};
  const auto                result = sym(values.data());
  constexpr auto            ref    = 11;
  EXPECT_DOUBLE_EQ(result, ref);
}

Symbols cal_normal(const Polynomials& parametric_curve)
{
  const auto var_index = 0;

  auto tangent      = ms::sym::get_differentiate(parametric_curve, var_index);
  auto unit_tangent = ms::sym::get_normalize(tangent);

  return ms::sym::get_differentiate(unit_tangent, var_index);
}

Symbol cal_curvature(const Polynomials& parametric_curve)
{
  const auto var_index = 0;

  auto       tangent      = ms::sym::get_differentiate(parametric_curve, var_index);
  auto       unit_tangent = ms::sym::get_normalize(tangent);
  const auto diff_ut      = ms::sym::get_differentiate(unit_tangent, var_index);

  auto curvature = ms::sym::cal_L2_norm(diff_ut) / ms::sym::cal_L2_norm(tangent);
  return curvature;
}

} // namespace ms::sym
