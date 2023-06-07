#pragma once
#include <vector>

#include "Symbol_Base.h"

// forward declaration
namespace ms::sym
{
class Multiplied_Term;
class Symbol;
class Polynomial;
class Polynomials;
using Symbols     = std::vector<Symbol>;

} // namespace ms::sym

// class declration
namespace ms::sym
{

class Powered_Term
{
public:
  Powered_Term(const std::unique_ptr<Base>& term_ptr, const double exponenet = 1.0);
  Powered_Term(std::unique_ptr<Base>&& term, const double exponent = 1.0)
      : base_ptr_(std::move(term)),
        exponent_(exponent){};
  Powered_Term(const Powered_Term& other);

public:
  void operator=(const Powered_Term& other);

public:
  void multiply_assign_with_same_base(const Powered_Term& other);

public:
  double operator()(const double* input) const;
  double operator()(const std::pair<const double*, int>& input) const;

public:
  Multiplied_Term get_differentiate(const int variable_index) const;
  bool            is_constant(void) const;
  std::string     to_string(void) const;
  double          to_constant(void) const;

private:
  std::unique_ptr<Base> base_ptr_ = nullptr;
  double                exponent_ = 1.0;
};

// constant * powered term1 * powered term2 * ...
class Multiplied_Term
{
public:
  Multiplied_Term(void) = default;
  Multiplied_Term(const double constant, Powered_Term&& pterm);
  Multiplied_Term(Powered_Term&& pterm);

public:
  void operator*=(const double constant);
  void operator*=(const Multiplied_Term& other);

public:
  void multiply_assign_powred_term(const Powered_Term& pterm);
  void multiply_assign_powred_term(Powered_Term&& pterm);

public:
  Multiplied_Term operator*(const double constant) const;
  Multiplied_Term operator*(const Multiplied_Term& other) const;
  double          operator()(const double* input) const;
  double          operator()(const std::pair<const double*, int>& input) const;

public:
  Symbol      get_differentiate(const int variable_index) const;
  bool        is_constant(void) const;
  double      to_constant(void) const;
  std::string to_string(void) const;

private:
  int                       num_term_ = 0;
  double                    constant_ = 0.0;
  std::vector<Powered_Term> pterms_;
};

// constant + Multiplied Term1 + Multiplied Term2 + ...
class Symbol : public Base
{
public:
  Symbol(void) = default;
  Symbol(const double constant) : constant_(constant){};
  Symbol(const Sym_Base& base_ptr);
  Symbol(Powered_Term&& pterm);

public:
  double operator()(const double* input) const override;
  double operator()(const std::pair<const double*, int>& input) const override;

public:
  Sym_Base    copy(void) const override;
  Sym_Base    get_differentiate(const int var_index) const override;
  bool        is_constant(void) const override;
  bool        is_zero(void) const override;
  double      to_constant(void) const override;
  std::string to_string(void) const override;

public:
  void operator+=(const Symbol& other);
  void operator-=(const Symbol& other);
  void operator*=(const double constant);

public:
  void add_assign_mterm(const Multiplied_Term& mterm);
  void add_assign_mterm(Multiplied_Term&& mterm);
  void be_absolute(void);
  void minus_assign_mterm(const Multiplied_Term& mterm);
  void pow(const double exponent);

public:
  Symbol operator+(const Symbol& other) const;
  Symbol operator-(const Symbol& other) const;
  Symbol operator*(const double constant) const;
  Symbol operator*(const Symbol& other) const;
  Symbol operator/(const Symbol& other) const;

public:
  Symbol   get_diff_symbol(const int variable_index) const;
  Symbol   get_pow(const double exponent) const;
  Sym_Base to_base_ptr(void) const;

private:
  std::vector<Multiplied_Term> mterms_;
  double                       constant_    = 0.0;
  bool                         is_absolute_ = false;
};

Multiplied_Term operator*(const double constant, const Multiplied_Term& mterm);

} // namespace ms::sym

// free function declaration
namespace ms::sym
{
Symbol      cal_L2_norm(const Polynomials& polys);
Symbol      cal_L2_norm(const Symbols& symbols);
Polynomials get_differentiate(const Polynomials& polys, const int var_index);
Symbols     get_differentiate(const Symbols& symbols, const int var_index);
Symbols     get_normalize(Polynomials& polys);
Symbols     get_normalize(const Symbols& symbols);
} // namespace ms::sym
