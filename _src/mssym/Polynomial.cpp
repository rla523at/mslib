#include "Polynomial.h"

#include "msexception/Exception.h"
#include <algorithm>
#include <iomanip>

namespace ms::sym
{

Simple_Poly_Term::Simple_Poly_Term(const std::string& variable)
{
  REQUIRE(variable.front() == 'x', "variable should be start with 'x'");

  constexpr int index_pos      = 1;
  const auto    variable_index = static_cast<int>(std::stoul(variable.substr(index_pos)));

  this->change_domain_dimension(variable_index + 1);
  this->coefficient_ptr_[variable_index] = 1.0;
}

Simple_Poly_Term::Simple_Poly_Term(const std::vector<double>& coefficients,
                                   const double               constant)
{
  this->change_domain_dimension(static_cast<int>(coefficients.size()));

  this->constant_ = constant;

  for (int i = 0; i < this->domain_dim_; ++i)
    this->coefficient_ptr_[i] += coefficients[i];
}

Simple_Poly_Term::Simple_Poly_Term(const Simple_Poly_Term& other)
{
  *this = other;
}

void Simple_Poly_Term::operator=(const Simple_Poly_Term& other)
{
  this->domain_dim_ = other.domain_dim_;
  this->constant_   = other.constant_;

  if (other.is_small())
  {
    this->coefficient_buffer_ = other.coefficient_buffer_;
    this->coefficient_ptr_    = this->coefficient_buffer_.data();
  }
  else
  {
    this->coefficients_    = other.coefficients_;
    this->coefficient_ptr_ = this->coefficients_.data();
  }
}

Simple_Poly_Term& Simple_Poly_Term::operator+=(const Simple_Poly_Term& other)
{
  if (this->domain_dim_ < other.domain_dim_)
    this->change_domain_dimension(other.domain_dim_);

  this->constant_ += other.constant_;

  for (int i = 0; i < other.domain_dim_; ++i)
    this->coefficient_ptr_[i] += other.coefficient_ptr_[i];

  if (this->is_constant())
    *this = this->constant_;

  return *this;
}

Simple_Poly_Term& Simple_Poly_Term::operator-=(const Simple_Poly_Term& other)
{
  if (this->domain_dim_ < other.domain_dim_)
    this->change_domain_dimension(other.domain_dim_);

  this->constant_ -= other.constant_;

  for (int i = 0; i < other.domain_dim_; ++i)
    this->coefficient_ptr_[i] -= other.coefficient_ptr_[i];

  return *this;
}

Simple_Poly_Term& Simple_Poly_Term::operator*=(const double constant)
{
  this->constant_ *= constant;

  for (int i = 0; i < this->domain_dim_; ++i)
    this->coefficient_ptr_[i] *= constant;

  return *this;
}

Simple_Poly_Term Simple_Poly_Term::operator+(
    const Simple_Poly_Term& other) const
{
  auto result = *this;
  return result += other;
}

Simple_Poly_Term Simple_Poly_Term::operator*(const double constant) const
{
  auto result = *this;
  return result *= constant;
}

Poly_Term Simple_Poly_Term::operator*(const Simple_Poly_Term& other) const
{
  Poly_Term result = *this;
  return result *= other;
}

bool Simple_Poly_Term::operator==(const Simple_Poly_Term& other) const
{
  if (this->domain_dim_ != other.domain_dim_)
    return false;

  if (this->constant_ != other.constant_)
    return false;

  for (int i = 0; i < this->domain_dim_; ++i)
  {
    if (this->coefficient_ptr_[i] != other.coefficient_ptr_[i])
      return false;
  }

  return true;
}

bool Simple_Poly_Term::operator!=(const Simple_Poly_Term& other) const
{
  return !(*this == other);
}

bool Simple_Poly_Term::operator<(const Simple_Poly_Term& other) const
{
  if (this->domain_dim_ != other.domain_dim_)
    return this->domain_dim_ < other.domain_dim_;

  for (int i = 0; i < this->domain_dim_; ++i)
  {
    if (this->coefficient_ptr_[i] != other.coefficient_ptr_[i])
      return this->coefficient_ptr_[i] < other.coefficient_ptr_[i];
  }

  return this->constant_ < other.constant_;
}

// bool SimplePolyTerm::operator>(const SimplePolyTerm& other) const
//{
//	if (this->coefficients_ == other.coefficients_)
//		return this->constant_ > other.constant_;
//	else
//		return this->coefficients_ > other.coefficients_;
// }

double Simple_Poly_Term::to_constant(void) const
{
  REQUIRE(this->is_constant(), "it should be constant");
  return this->constant_;
}

double Simple_Poly_Term::get_differentiate(const int variable_index) const
{
  if (this->domain_dim_ <= variable_index)
    return 0;
  else
    return this->coefficient_ptr_[variable_index];
}

int Simple_Poly_Term::degree(void) const
{
  if (this->is_constant())
    return 0;
  else
    return 1;
}

int Simple_Poly_Term::domain_dimension(void) const { return this->domain_dim_; }

bool Simple_Poly_Term::is_constant(void) const
{
  for (int i = 0; i < this->domain_dim_; ++i)
  {
    if (this->coefficient_ptr_[i] != 0)
      return false;
  }

  return true;
}

bool Simple_Poly_Term::is_zero(void) const
{
  return this->is_constant() && this->constant_ == 0.0;
}

std::string Simple_Poly_Term::to_string(void) const
{
  std::ostringstream os;
  os << std::setprecision(16) << std::showpos;

  os << "[";
  for (int i = 0; i < this->domain_dim_; ++i)
  {
    if (this->coefficient_ptr_[i] == 0.0)
      continue;
    else if (this->coefficient_ptr_[i] == 1.0)
      os << "+x" << std::to_string(i);
    else if (this->coefficient_ptr_[i] == -1.0)
      os << "-x" << std::to_string(i);
    else
      os << this->coefficient_ptr_[i] << "(x" << std::to_string(i) << ")";
  }

  if (this->constant_ == 0.0)
    os << "]";
  else
    os << this->constant_ << "]";

  auto str = os.str();

  constexpr int position = 1;
  constexpr int size     = 1;
  if (str.at(position) == '+')
    str.erase(position, size);

  return str;
}

void Simple_Poly_Term::change_domain_dimension(const int new_domain_dimension)
{
  if (this->is_small())
  {
    if (new_domain_dimension > this->small_criterion_)
    {
      this->coefficients_.resize(new_domain_dimension);
      this->coefficient_ptr_ = this->coefficients_.data(); // resize can cuase reallcoation
                                                           // -> ptr should be updated

      std::copy(this->coefficient_buffer_.begin(), coefficient_buffer_.end(),
                this->coefficients_.begin());
      this->coefficient_buffer_.fill(0);
    }
  }
  else
  {
    this->coefficients_.resize(new_domain_dimension);
    this->coefficient_ptr_ = this->coefficients_
                                 .data(); // resize can cuase reallcoation -> ptr should be updated
  }

  this->domain_dim_ = new_domain_dimension;
}

bool Simple_Poly_Term::is_small(void) const { return this->coefficients_.empty(); }

void Powered_Poly_Term::multiply_assign_with_same_base(
    const Powered_Poly_Term& other)
{
  this->exponent_ += other.exponent_;
}

Poly_Term Powered_Poly_Term::operator*(const double constant) const
{
  return {constant, *this};
}

bool Powered_Poly_Term::operator==(const Powered_Poly_Term& other) const
{
  return this->base_ == other.base_ && this->exponent_ == other.exponent_;
}

bool Powered_Poly_Term::operator<(const Powered_Poly_Term& other) const
{
  if (this->exponent_ == other.exponent_)
    return this->base_ < other.base_;
  else
    return this->exponent_ < other.exponent_;
}

// bool PoweredPolyTerm::operator>(const PoweredPolyTerm& other) const
//{
//	if (this->exponent_ == other.exponent_)
//		return this->base_ > other.base_;
//	else
//		return this->exponent_ > other.exponent_;
// }

double Powered_Poly_Term::to_constant(void) const
{
  return std::pow(this->base_.to_constant(), this->exponent_);
}

Simple_Poly_Term Powered_Poly_Term::get_simple(void) const
{
  return this->base_;
}

Poly_Term Powered_Poly_Term::get_differentiate(const int variable_index) const
{
  const auto base_derivative = this->base_.get_differentiate(variable_index);

  if (base_derivative == 0.0)
    return 0.0;

  if (this->exponent_ == 1)
    return base_derivative;
  else
    return this->exponent_ * Powered_Poly_Term(this->base_, this->exponent_ - 1) * base_derivative;
}

bool Powered_Poly_Term::has_same_base(const Powered_Poly_Term& other) const
{
  return this->base_ == other.base_;
}

int Powered_Poly_Term::domain_dimension(void) const
{
  return this->base_.domain_dimension();
}

bool Powered_Poly_Term::is_constant(void) const
{
  return this->base_.is_constant();
}

bool Powered_Poly_Term::is_simple(void) const { return this->exponent_ == 1; }

int Powered_Poly_Term::degree(void) const
{
  if (this->is_constant())
    return 0;
  else
    return this->exponent_;
}

std::string Powered_Poly_Term::to_string(void) const
{
  auto str = this->base_.to_string();
  if (this->exponent_ != 1)
    return str + "^" + std::to_string(this->exponent_);
  else
    return str;
}

// PolyTerm::PolyTerm(const double constant) {
//	this->constant_ = constant;
//	this->add_term(1.0);
// }

Poly_Term::Poly_Term(const Simple_Poly_Term& simple_poly_term)
{
  if (simple_poly_term.is_constant())
  {
    this->constant_ = simple_poly_term.to_constant();
    this->add_term(1.0);
  }
  else
  {
    this->constant_ = 1.0;
    this->add_term(simple_poly_term);
  }
}

Poly_Term::Poly_Term(const Powered_Poly_Term& powered_poly_term)
{
  if (powered_poly_term.is_constant())
  {
    this->constant_ = powered_poly_term.to_constant();
    this->add_term(1.0);
  }
  else
  {
    this->constant_ = 1.0;
    this->add_term(powered_poly_term);
  }
}

Poly_Term::Poly_Term(const double             constant,
                     const Powered_Poly_Term& powered_poly_term)
{
  if (powered_poly_term.is_constant())
  {
    this->constant_ = constant * powered_poly_term.to_constant();
    this->add_term(1.0);
  }
  else
  {
    this->constant_ = constant;
    this->add_term(powered_poly_term);
  }
}

Poly_Term::Poly_Term(const Poly_Term& other)
{
  this->constant_ = other.constant_;
  this->num_term_ = other.num_term_;
  if (other.is_small())
  {
    this->small_buffer_ = other.small_buffer_;
    this->term_ptr_     = this->small_buffer_.data();
  }
  else
  {
    this->terms_    = other.terms_;
    this->term_ptr_ = this->terms_.data();
  }
}

void Poly_Term::add_assign_with_same_form(const Poly_Term& other)
{
  this->constant_ += other.constant_;
  if (this->constant_ == 0.0)
    *this = 0.0;
}

void Poly_Term::minus_assign_with_same_form(const Poly_Term& other)
{
  this->constant_ -= other.constant_;
  if (this->constant_ == 0.0)
    *this = 0.0;
}

Poly_Term& Poly_Term::operator*=(const double constant)
{
  this->constant_ *= constant;

  if (this->constant_ == 0.0)
    *this = 0.0;

  return *this;
}

Poly_Term& Poly_Term::operator*=(const Poly_Term& other)
{
  *this *= other.constant_;

  for (int i = 0; i < other.num_term_; ++i)
    this->multiply_assign_powered_poly_term(other.term_ptr_[i]);

  // for Poly_Term::has_same_form
  if (this->is_small())
    std::sort(this->small_buffer_.begin(),
              this->small_buffer_.begin() + this->num_term_);
  else
    std::sort(this->terms_.begin(), this->terms_.end());

  return *this;
}

Poly_Term Poly_Term::operator*(const double constant) const
{
  auto result = *this;
  return result *= constant;
}

Poly_Term Poly_Term::operator*(const Poly_Term& other) const
{
  auto result = *this;
  return result *= other;
}

bool Poly_Term::operator==(const Poly_Term& other) const
{
  if (this->constant_ != other.constant_)
    return false;

  return this->has_same_form(other);
}

bool Poly_Term::operator!=(const Poly_Term& other) const
{
  return !(*this == other);
}

Poly_Term& Poly_Term::operator=(const Poly_Term& other)
{
  this->constant_ = other.constant_;
  this->num_term_ = other.num_term_;

  if (other.is_small())
  {
    this->small_buffer_ = other.small_buffer_;
    this->term_ptr_     = this->small_buffer_.data();
  }
  else
  {
    this->terms_    = other.terms_;
    this->term_ptr_ = this->terms_.data();
  }

  return *this;
}

Simple_Poly_Term Poly_Term::be_simple(void) const
{
  return this->term_ptr_[0].get_simple() * this->constant_;
}

Polynomial Poly_Term::get_differentiate(const int variable_index) const
{
  Polynomial result = 0.0;

  for (int i = 0; i < this->num_term_; ++i)
  {
    auto derivative = this->term_ptr_[i].get_differentiate(variable_index) * this->constant_;

    if (derivative.is_zero())
      continue;

    for (int j = 0; j < this->num_term_; ++j)
    {
      if (j == i)
        continue;
      else
        derivative.add_term(this->term_ptr_[j]);
    }

    if (derivative.is_simple())
      result += derivative.be_simple();
    else
      result += derivative;
  }

  return result;
}

double Poly_Term::to_constant(void) const { return this->constant_; }

int Poly_Term::degree(void) const
{
  int result = 0;

  for (int i = 0; i < this->num_term_; ++i)
    result += this->term_ptr_[i].degree();

  return result;
}

int Poly_Term::domain_dimension(void) const
{
  int domain_dimension = 0;

  for (int i = 0; i < this->num_term_; ++i)
    domain_dimension = std::max(domain_dimension, this->term_ptr_[i].domain_dimension());

  return domain_dimension;
}

bool Poly_Term::has_same_form(const Poly_Term& other) const
{
  if (this->num_term_ != other.num_term_)
    return false;

  if (this->is_small())
    return this->small_buffer_ == other.small_buffer_;
  else
    return this->terms_ == other.terms_;
}

bool Poly_Term::is_simple(void) const
{
  if (this->num_term_ != 1)
    return false;

  return this->term_ptr_[0].is_simple();
}

bool Poly_Term::is_zero(void) const { return this->constant_ == 0.0; }

std::string Poly_Term::to_string(void) const
{
  std::ostringstream oss;
  oss << std::setprecision(16) << std::showpos;

  if (std::abs(this->constant_) != 1.0)
    oss << this->constant_;
  else if (this->constant_ == 1.0)
    oss << "+";
  else
    oss << "-";

  for (int i = 0; i < this->num_term_; ++i)
    oss << this->term_ptr_[i].to_string();

  return oss.str();
}

void Poly_Term::add_term(const Powered_Poly_Term& powered_poly_term)
{
  if (this->is_constant())
  {
    *this->term_ptr_ = powered_poly_term;
    return;
  }

  const auto new_term_pos = this->num_term_;
  this->num_term_++;

  if (this->is_small())
  {
    if (this->small_criterion_ < this->num_term_)
    {
      this->terms_.resize(this->num_term_); // resize can cuase reallcoation ->
                                            // ptr should be updated
      this->term_ptr_ = this->terms_.data();

      std::copy(this->small_buffer_.begin(), this->small_buffer_.end(),
                this->terms_.begin());
      this->small_buffer_.fill(0.0);
    }
  }
  else
  {
    this->terms_.resize(this->num_term_); // resize can cuase reallcoation ->
                                          // ptr should be updated
    this->term_ptr_ = this->terms_.data();
  }

  this->term_ptr_[new_term_pos] = powered_poly_term;
}

void Poly_Term::multiply_assign_powered_poly_term(
    const Powered_Poly_Term& power_poly_term)
{
  if (power_poly_term == 0.0)
    *this = 0.0;

  if (power_poly_term.is_constant())
    return;

  for (int i = 0; i < this->num_term_; ++i)
  {
    if (this->term_ptr_[i].has_same_base(power_poly_term))
    {
      this->term_ptr_[i].multiply_assign_with_same_base(power_poly_term);
      return;
    }
  }

  this->add_term(power_poly_term);
}

bool Poly_Term::is_constant(void) const
{
  return this->num_term_ == 1 && this->term_ptr_->is_constant();
}

bool Poly_Term::is_small(void) const { return this->terms_.empty(); }

Polynomial::Polynomial(const Simple_Poly_Term& simple_poly_term)
{
  this->simple_poly_term_ += simple_poly_term;
}

Polynomial::Polynomial(const Poly_Term& poly_term)
{
  this->add_assign_poly_term(poly_term);
}

Sym_Base Polynomial::copy(void) const
{
  return std::make_unique<Polynomial>(*this);
}

Sym_Base Polynomial::get_differentiate(const int var_index) const
{
  return std::make_unique<Polynomial>(this->get_diff_polynomial(var_index));
}

Polynomial& Polynomial::operator+=(const Polynomial& other)
{
  this->simple_poly_term_ += other.simple_poly_term_;
  for (const auto& poly_term : other.poly_terms_)
    this->add_assign_poly_term(poly_term);

  return *this;
}

Polynomial& Polynomial::operator-=(const Polynomial& other)
{
  this->simple_poly_term_ -= other.simple_poly_term_;
  for (const auto& poly_term : other.poly_terms_)
    this->minus_assign_poly_term(poly_term);

  return *this;
}

Polynomial& Polynomial::operator*=(const double constant)
{
  if (constant == 0.0)
    return *this = 0.0;

  this->simple_poly_term_ *= constant;
  for (auto& poly_term : this->poly_terms_)
    poly_term *= constant;

  return *this;
}

Polynomial Polynomial::operator+(const Polynomial& other) const
{
  Polynomial result(*this);
  return result += other;
}

Polynomial Polynomial::operator-(const Polynomial& other) const
{
  Polynomial result(*this);
  return result -= other;
}

Polynomial Polynomial::operator*(const Polynomial& other) const
{
  Polynomial result = 0.0;

  const auto num_this_term  = this->poly_terms_.size();
  const auto num_other_term = other.poly_terms_.size();
  for (int i = 0; i < num_this_term; ++i)
    for (int j = 0; j < num_other_term; ++j)
      result.add_assign_poly_term(this->poly_terms_[i] * other.poly_terms_[j]);

  if (this->simple_poly_term_ != 0.0)
  {
    for (int j = 0; j < num_other_term; ++j)
      result.add_assign_poly_term(other.poly_terms_[j] * this->simple_poly_term_);
  }

  if (other.simple_poly_term_ != 0.0)
  {
    for (int i = 0; i < num_this_term; ++i)
      result.add_assign_poly_term(this->poly_terms_[i] * other.simple_poly_term_);
  }

  if (this->simple_poly_term_ != 0.0 && other.simple_poly_term_ != 0.0)
  {
    if (this->simple_poly_term_.is_constant())
    {
      result.simple_poly_term_ = other.simple_poly_term_ * this->simple_poly_term_.to_constant();
    }
    else if (other.simple_poly_term_.is_constant())
    {
      result.simple_poly_term_ = this->simple_poly_term_ * other.simple_poly_term_.to_constant();
    }
    else
    {
      result.add_assign_poly_term(Poly_Term(this->simple_poly_term_) * other.simple_poly_term_);
    }
  }

  return result;
}

Polynomial Polynomial::operator*(const double constant) const
{
  Polynomial result = *this;
  return result *= constant;
}

Polynomial Polynomial::operator^(const int power_index) const
{
  if (power_index == 0)
    return 1.0;

  auto result = *this;
  for (int i = 1; i < power_index; ++i)
    result = std::move(result * *this);
  // result *= *this;

  return result;
}

bool Polynomial::operator==(const Polynomial& other) const
{
  const auto max_degree           = std::max(this->degree(), other.degree());
  const auto max_domain_dimension = std::max(this->domain_dimension(), other.domain_dimension());
  const auto compare_node_set     = ms::sym::polynomial_compare_node_set(max_degree, max_domain_dimension);

  for (const auto& compare_node : compare_node_set)
  {
    if (!ms::sym::compare_double((*this)(compare_node.data()), other(compare_node.data())))
      return false;
  }

  return true;
}

Polynomial& Polynomial::differentiate(const int variable_index)
{
  auto result  = this->get_diff_polynomial(variable_index);
  return *this = std::move(result);
};

int Polynomial::degree(void) const
{
  int result = this->simple_poly_term_.degree();
  for (const auto& term : this->poly_terms_)
    result = std::max(result, term.degree());
  return result;
}

int Polynomial::domain_dimension(void) const
{
  int domain_dimension = this->simple_poly_term_.domain_dimension();

  for (const auto& term : poly_terms_)
    domain_dimension = std::max(domain_dimension, term.domain_dimension());

  return domain_dimension;
}

Polynomial Polynomial::get_diff_polynomial(const int var_index) const
{
  Polynomial result = 0.0;

  if (this->domain_dimension() <= var_index)
    return result;

  result.simple_poly_term_ = this->simple_poly_term_.get_differentiate(var_index);

  for (const auto& poly_term : this->poly_terms_)
  {
    result += poly_term.get_differentiate(var_index);
  }

  return result;
}

std::vector<Polynomial> Polynomial::gradient(void) const
{
  const auto              domain_dimension = this->domain_dimension();
  std::vector<Polynomial> gradient(domain_dimension);

  for (int i = 0; i < domain_dimension; ++i)
  {
    gradient[i] = this->get_diff_polynomial(i);
  }

  return gradient;
}

std::vector<Polynomial> Polynomial::gradient(const int space_dimension) const
{
  std::vector<Polynomial> gradient(space_dimension);

  for (int i = 0; i < space_dimension; ++i)
  {
    gradient[i] = this->get_diff_polynomial(i);
  }

  return gradient;
}

size_t Polynomial::num_term(void) const
{
  return this->poly_terms_.size();
};

Sym_Base Polynomial::to_sym_base(void) const
{
  return std::make_unique<Polynomial>(*this);
}

double Polynomial::to_constant(void) const
{
  REQUIRE(this->is_constant(), "it should be constant");

  return this->simple_poly_term_.to_constant();
};

std::string Polynomial::to_string(void) const
{
  if (this->poly_terms_.empty())
  {
    return this->simple_poly_term_.to_string();
  }

  std::string str;
  for (const auto& poly_term : this->poly_terms_)
  {
    str += poly_term.to_string();
  }

  if (!(this->simple_poly_term_ == 0))
  {
    str += "+" + this->simple_poly_term_.to_string();
  }

  if (str.front() == '+')
  {
    str.erase(str.begin());
  }

  return str;
}

bool Polynomial::is_constant(void) const
{
  return (this->poly_terms_.empty() && this->simple_poly_term_.is_constant());
}

bool Polynomial::is_zero(void) const
{
  for (const auto& term : this->poly_terms_)
  {
    if (not term.is_zero())
      return false;
  }

  return simple_poly_term_.is_zero();
}

void Polynomial::add_assign_poly_term(const Poly_Term& term)
{
  auto iter = this->poly_terms_.begin();
  for (; iter != this->poly_terms_.end(); ++iter)
  {
    if (iter->has_same_form(term))
    {
      iter->add_assign_with_same_form(term);

      if (iter->is_zero())
        this->poly_terms_.erase(iter);
      return;
    }
  }

  this->poly_terms_.push_back(term);
}

void Polynomial::minus_assign_poly_term(const Poly_Term& term)
{
  for (auto iter = this->poly_terms_.begin(); iter != this->poly_terms_.end();
       ++iter)
  {
    if (iter->has_same_form(term))
    {
      iter->minus_assign_with_same_form(term);
      if (iter->is_zero())
        this->poly_terms_.erase(iter);
      return;
    }
  }
  this->poly_terms_.push_back(-1 * term);
}

std::ostream& operator<<(std::ostream&           ostream,
                         const Simple_Poly_Term& simple_poly_term)
{
  return ostream << simple_poly_term.to_string();
}

std::ostream& operator<<(std::ostream& ostream, const Poly_Term& poly_term)
{
  return ostream << poly_term.to_string();
}

std::ostream& operator<<(std::ostream& ostream, const Polynomial& polynomial)
{
  return ostream << polynomial.to_string();
}

Simple_Poly_Term operator*(const double            constant,
                           const Simple_Poly_Term& simple_poly_term)
{
  return simple_poly_term * constant;
}

Poly_Term operator*(const double             constant,
                    const Powered_Poly_Term& powered_poly_term)
{
  return powered_poly_term * constant;
}

Poly_Term operator*(const double constant, const Poly_Term& poly_term)
{
  return poly_term * constant;
}

Polynomial operator*(const double constant, const Polynomial& polynomial)
{
  return polynomial * constant;
}

Polynomial operator+(const double constant, const Polynomial& polynomial)
{
  return polynomial + constant;
}

Polynomial operator-(const double constant, const Polynomial& polynomial)
{
  return (-1.0 * polynomial) + constant;
}

double Simple_Poly_Term::operator()(const double* input) const
{
  auto result = this->constant_;
  for (int i = 0; i < this->domain_dim_; ++i)
    result += this->coefficient_ptr_[i] * input[i];

  return result;
}

double Powered_Poly_Term::operator()(const double* input) const
{
  return std::pow(this->base_(input), this->exponent_);
}

double Poly_Term::operator()(const double* input) const
{
  auto result = this->constant_;
  for (int i = 0; i < this->num_term_; ++i)
    result *= this->term_ptr_[i](input);
  return result;
}

double Polynomial::operator()(const double* input) const
{
  auto result = this->simple_poly_term_(input);
  for (const auto& poly_term : this->poly_terms_)
  {
    result += poly_term(input);
  }

  return result;
}

int combination(const int n, const int k)
{
  // calculate nCk
  // the combination of n things taken k at a time without repetition.
  if (n == k || k == 0)
    return 1;
  else
    return combination(n - 1, k - 1) + combination(n - 1, k);
}

int combination_with_repetition(const int n, const int k)
{
  // calculate nHk
  // the combination of n things taken k at a time with repetition.
  return combination(n + k - 1, k);
}

bool is_positive_odd_number(const double val)
{
  if (val < 0)
    return false;

  if (val - std::floor(val) == 0)
    return static_cast<size_t>(val) % 2 == 0;
  else
    return false;
}

bool is_natural_number(const double val)
{
  if (val < 0)
    return false;

  if (val - std::floor(val) == 0)
    return true;
  else
    return false;
}

bool compare_double(const double d1, const double d2,
                    const size_t ULP_precision)
{
  const auto lower_ULP = d1 - std::nextafter(d1, std::numeric_limits<double>::lowest());
  const auto upper_ULP = std::nextafter(d1, std::numeric_limits<double>::max()) - d1;

  return d1 - ULP_precision * lower_ULP <= d2 && d2 <= d1 + ULP_precision * upper_ULP;
}

std::vector<std::vector<double>> polynomial_compare_node_set(
    const int polynomial_order, const int domain_dimension)
{
  const auto num_node = ms::sym::combination_with_repetition(polynomial_order + 1, domain_dimension);

  std::vector<std::vector<double>> compare_node_set;
  compare_node_set.reserve(num_node);

  std::vector<double> compare_node(domain_dimension);
  if (domain_dimension == 0)
  {
    compare_node_set.push_back(compare_node);
    return compare_node_set;
  }

  while (true)
  {
    auto iter = std::find(compare_node.begin(), compare_node.end(), polynomial_order);
    if (iter != compare_node.end())
    {
      compare_node_set.push_back(compare_node);

      if (iter == compare_node.begin())
        break;

      std::fill(compare_node.begin(), compare_node.end(), 0);
      (*(--iter))++;

      if (compare_node.front() == polynomial_order)
      {
        compare_node_set.push_back(compare_node);
        break;
      }
    }

    double component_sum = 0;
    for (const auto& val : compare_node)
      component_sum += val;

    if (component_sum == polynomial_order)
    {
      compare_node_set.push_back(compare_node);
      const auto is_zero = [](const double i)
      { return i == 0; };
      auto iter2 = std::find_if_not(compare_node.rbegin(), compare_node.rend(), is_zero);
      *iter2     = 0;
      (*(++iter2))++;
      continue;
    }

    compare_node_set.push_back(compare_node);
    compare_node.back()++;
  }

  return compare_node_set;
}

} // namespace ms::sym
