#pragma once
#include <array>
#include <sstream>
#include <vector>

#include "Symbol_Base.h"

// forwrad declaration
namespace ms::sym
{
class Poly_Term;
class Polynomial;
using Polynomials = std::vector<Polynomial>;
} // namespace ms::sym

// class declration
namespace ms::sym
{

// Linear combination of monomial degree 1
class Simple_Poly_Term
{
public:
  Simple_Poly_Term(const double constant)
      : constant_(constant){};
  Simple_Poly_Term(const std::string& variable);
  Simple_Poly_Term(const std::vector<double>& coeffs, const double c = 0);
  Simple_Poly_Term(const Simple_Poly_Term& other);

public:
  Simple_Poly_Term& operator+=(const Simple_Poly_Term& other);
  Simple_Poly_Term& operator-=(const Simple_Poly_Term& other);
  Simple_Poly_Term& operator*=(const double constant);
  void              operator=(const Simple_Poly_Term& other);

public:
  Simple_Poly_Term operator+(const Simple_Poly_Term& other) const;
  Simple_Poly_Term operator*(const double constant) const;
  Poly_Term        operator*(const Simple_Poly_Term& other) const;
  double           operator()(const double* input) const;
  bool             operator==(const Simple_Poly_Term& other) const;
  bool             operator!=(const Simple_Poly_Term& other) const;
  bool             operator<(const Simple_Poly_Term& other) const;

public:
  double      to_constant(void) const;
  double      get_differentiate(const int variable_index) const;
  int         degree(void) const;
  int         domain_dimension(void) const;
  bool        is_constant(void) const;
  bool        is_zero(void) const;
  std::string to_string(void) const;

private:
  void change_domain_dimension(const int new_domain_dimension);
  bool is_small(void) const;

private:
  static constexpr int small_criterion_ = 4;

  int                                  domain_dim_         = 0;
  double                               constant_           = 0.0;
  double*                              coefficient_ptr_    = this->coefficient_buffer_.data();
  std::array<double, small_criterion_> coefficient_buffer_ = {0};
  std::vector<double>                  coefficients_;
};

class Powered_Poly_Term
{
public:
  Powered_Poly_Term(void) = default;
  Powered_Poly_Term(const double constant)
      : base_(constant){};
  Powered_Poly_Term(const Simple_Poly_Term& simple_poly_term)
      : base_(simple_poly_term){};
  Powered_Poly_Term(const Simple_Poly_Term& simple_poly_term, const int exponent)
      : base_(simple_poly_term), exponent_(exponent){};

public:
  void multiply_assign_with_same_base(const Powered_Poly_Term& other);

public:
  Poly_Term operator*(const double constant) const;
  double    operator()(const double* input) const;
  bool      operator==(const Powered_Poly_Term& other) const;
  bool      operator<(const Powered_Poly_Term& other) const; // for Poly_Term::has_same_form

public:
  Poly_Term        get_differentiate(const int variable_index) const;
  int              degree(void) const;
  int              domain_dimension(void) const;
  double           to_constant(void) const;
  Simple_Poly_Term get_simple(void) const;
  bool             has_same_base(const Powered_Poly_Term& other) const;
  bool             is_constant(void) const;
  bool             is_simple(void) const;
  std::string      to_string(void) const;

private:
  Simple_Poly_Term base_     = 0.0;
  int              exponent_ = 1;
};

// constant * powered poly term1 * powered poly term2 * ...
class Poly_Term
{
public:
  Poly_Term(const Simple_Poly_Term& simple_poly_term);
  Poly_Term(const Powered_Poly_Term& powered_poly_term);
  Poly_Term(const double constant, const Powered_Poly_Term& powered_poly_term = 1.0);
  Poly_Term(const Poly_Term& other);

public:
  Poly_Term& operator*=(const double constant);
  Poly_Term& operator*=(const Poly_Term& other);
  Poly_Term& operator=(const Poly_Term& other);

public:
  void add_assign_with_same_form(const Poly_Term& other);
  void minus_assign_with_same_form(const Poly_Term& other);

public:
  Poly_Term operator*(const double constant) const;
  Poly_Term operator*(const Poly_Term& other) const;
  double    operator()(const double* input) const;
  bool      operator==(const Poly_Term& other) const;
  bool      operator!=(const Poly_Term& other) const;

public:
  Simple_Poly_Term be_simple(void) const;
  int              degree(void) const;
  int              domain_dimension(void) const;
  Polynomial       get_differentiate(const int variable_index) const;
  double           to_constant(void) const;
  bool             has_same_form(const Poly_Term& other) const;
  bool             is_simple(void) const;
  bool             is_zero(void) const;
  std::string      to_string(void) const;

private:
  void add_term(const Powered_Poly_Term& powered_poly_term);
  void multiply_assign_powered_poly_term(const Powered_Poly_Term& power_poly_term);
  bool is_constant(void) const;
  bool is_small(void) const;

private:
  static constexpr int                            small_criterion_ = 3;
  int                                             num_term_        = 0;
  double                                          constant_        = 0.0;
  Powered_Poly_Term*                              term_ptr_        = small_buffer_.data();
  std::array<Powered_Poly_Term, small_criterion_> small_buffer_    = {0};
  std::vector<Powered_Poly_Term>                  terms_;
};

// simple poly term + poly term1 + poly term2 + ...
class Polynomial : public Base
{
public:
  Polynomial(void) = default;
  Polynomial(const double coeeficient)
      : simple_poly_term_(coeeficient){};
  Polynomial(const std::string& variable)
      : simple_poly_term_(variable){};
  Polynomial(const Simple_Poly_Term& simple_poly_term);
  Polynomial(const Poly_Term& poly_term);

public:
  double operator()(const double* input) const override;

public:
  Sym_Base    copy(void) const override;
  Sym_Base    get_differentiate(const int var_index) const override;
  bool        is_constant(void) const override;
  bool        is_zero(void) const override;
  double      to_constant(void) const override;
  std::string to_string(void) const override;

public:
  Polynomial& operator+=(const Polynomial& other);
  Polynomial& operator-=(const Polynomial& other);
  Polynomial& operator*=(const double constant);

public:
  Polynomial& differentiate(const int variable_index);

public:
  Polynomial operator+(const Polynomial& other) const;
  Polynomial operator-(const Polynomial& other) const;
  Polynomial operator*(const Polynomial& other) const;
  Polynomial operator*(const double constant) const;
  Polynomial operator^(const int power_index) const;
  bool       operator==(const Polynomial& other) const;

public:
  int                     degree(void) const;
  int                     domain_dimension(void) const;
  Polynomial              get_diff_polynomial(const int var_index) const;
  std::vector<Polynomial> gradient(void) const;
  std::vector<Polynomial> gradient(const int space_dimension) const;
  size_t                  num_term(void) const;
  Sym_Base                to_sym_base(void) const;

private:
  void add_assign_poly_term(const Poly_Term& term);
  void minus_assign_poly_term(const Poly_Term& term);

private:
  Simple_Poly_Term       simple_poly_term_ = 0.0;
  std::vector<Poly_Term> poly_terms_;
};

} // namespace ms::sym

// free function declaration
namespace ms::sym
{
Simple_Poly_Term                 operator*(const double constant, const Simple_Poly_Term& simple_poly_term);
Poly_Term                        operator*(const double constant, const Powered_Poly_Term& powered_poly_term);
Poly_Term                        operator*(const double constant, const Poly_Term& poly_term);
Polynomial                       operator*(const double constant, const Polynomial& polynomial);
Polynomial                       operator+(const double constant, const Polynomial& polynomial);
Polynomial                       operator-(const double constant, const Polynomial& polynomial);
std::ostream&                    operator<<(std::ostream& ostream, const Simple_Poly_Term& simple_poly_term);
std::ostream&                    operator<<(std::ostream& ostream, const Poly_Term& poly_term);
std::ostream&                    operator<<(std::ostream& ostream, const Polynomial& polynomial);
int                              combination(const int n, const int k);
int                              combination_with_repetition(const int n, const int k);
bool                             is_positive_odd_number(const double val);
bool                             is_natural_number(const double val);
bool                             compare_double(const double d1, const double d2, const size_t ULP_precision = 4);
std::vector<std::vector<double>> polynomial_compare_node_set(const int polynomial_order, const int domain_dimension);
} // namespace ms::sym