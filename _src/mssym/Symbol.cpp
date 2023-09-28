#include "mssym/Symbol.h"

#include "msexception/Exception.h"
#include "mssym/Polynomial.h"
#include <iomanip>

namespace ms::sym
{

Powered_Term::Powered_Term(const std::unique_ptr<Base>& term_ptr, const double exponenet)
{
  this->base_ptr_ = term_ptr->copy();
  this->exponent_ = exponenet;
}

Powered_Term::Powered_Term(const Powered_Term& other)
{
  this->base_ptr_ = other.base_ptr_->copy();
  this->exponent_ = other.exponent_;
}

void Powered_Term::operator=(const Powered_Term& other)
{
  this->base_ptr_ = other.base_ptr_->copy();
  this->exponent_ = other.exponent_;
}

void Powered_Term::multiply_assign_with_same_base(const Powered_Term& other)
{
  this->exponent_ += other.exponent_;
}

double Powered_Term::operator()(const double* input) const
{
  return std::pow((*this->base_ptr_)(input), this->exponent_);
}

double Powered_Term::operator()(const std::pair<const double*, int>& input) const
{
  return std::pow((*this->base_ptr_)(input), this->exponent_);
}

double Powered_Term::operator()(const ms::math::Vector_View input) const
{
  return std::pow((*this->base_ptr_)(input), this->exponent_);
}

Multiplied_Term Powered_Term::get_differentiate(const int var_index) const
{
  if (this->is_constant())
  {
    return Multiplied_Term();
  }

  auto base_derivative = this->base_ptr_->get_differentiate(var_index);

  if (this->exponent_ == 1)
  {
    Powered_Term result(std::move(base_derivative));
    return result;
  }
  else
  {
    Powered_Term temp(*this);
    temp.exponent_ = this->exponent_ - 1;

    Multiplied_Term result(this->exponent_, std::move(temp));
    result.multiply_assign_powred_term(std::move(base_derivative));
    return result;
  }
}

std::string Powered_Term::to_string(void) const
{
  auto str = this->base_ptr_->to_string();
  if (this->exponent_ != 1)
  {
    return str + "^" + std::to_string(this->exponent_);
  }
  else
  {
    return str;
  }
}

double Powered_Term::to_constant(void) const
{
  if (this->exponent_ == 0.0)
    return 1.0;

  return std::pow(this->base_ptr_->to_constant(), this->exponent_);
}

bool Powered_Term::is_constant(void) const
{
  return this->base_ptr_->is_constant() || this->exponent_ == 0.0;
}

Multiplied_Term::Multiplied_Term(const double constant, Powered_Term&& pterm)
{
  this->_constant = constant;
  pterms_.push_back(std::move(pterm));
  this->num_term_++;
}

Multiplied_Term::Multiplied_Term(Powered_Term&& powered_term)
{
  this->_constant = 1.0;
  pterms_.push_back(std::move(powered_term));
  this->num_term_++;
}

void Multiplied_Term::operator*=(const double constant)
{
  if (constant == 0.0)
  {
    *this = Multiplied_Term();
  }
  else
  {
    this->_constant *= constant;
  }
}

void Multiplied_Term::operator*=(const Multiplied_Term& other)
{
  if (other.is_constant())
  {
    this->_constant *= other.to_constant();
  }
  else
  {
    this->_constant *= other._constant;

    this->num_term_ += other.num_term_;
    this->pterms_.reserve(this->num_term_);
    this->pterms_.insert(this->pterms_.end(), other.pterms_.begin(),
                         other.pterms_.end());
  }
}

void Multiplied_Term::multiply_assign_powred_term(const Powered_Term& pterm)
{
  if (pterm.is_constant())
  {
    this->_constant *= pterm.to_constant();
  }
  else
  {
    this->pterms_.push_back(pterm);
    this->num_term_++;
  }
}

void Multiplied_Term::multiply_assign_powred_term(Powered_Term&& pterm)
{
  if (pterm.is_constant())
  {
    this->_constant *= pterm.to_constant();
  }
  else
  {
    this->pterms_.push_back(std::move(pterm));
    this->num_term_++;
  }
}

Multiplied_Term Multiplied_Term::operator*(const double constant) const
{
  if (constant == 0.0)
  {
    return Multiplied_Term();
  }

  auto result = *this;
  result *= constant;
  return result;
}

Multiplied_Term Multiplied_Term::operator*(const Multiplied_Term& other) const
{
  auto result = *this;
  result *= other;
  return result;
}

double Multiplied_Term::operator()(const double* input) const
{
  double result = this->_constant;
  for (auto& term : this->pterms_)
  {
    result *= term(input);
  }

  return result;
}

double Multiplied_Term::operator()(const std::pair<const double*, int>& input) const
{
  double result = this->_constant;
  for (auto& term : this->pterms_)
  {
    result *= term(input);
  }

  return result;
}

double Multiplied_Term::operator()(const ms::math::Vector_View input) const
{
  double result = this->_constant;
  for (auto& term : this->pterms_)
  {
    result *= term(input);
  }

  return result;
}

Symbol Multiplied_Term::get_differentiate(const int variable_index) const
{
  if (this->is_constant())
    return Symbol();

  Symbol result;

  for (int i = 0; i < this->num_term_; ++i)
  {
    const auto& pterm = this->pterms_[i];

    if (pterm.is_constant())
      continue;

    auto derivative = pterm.get_differentiate(variable_index);

    derivative *= this->_constant;

    for (int j = 0; j < this->num_term_; ++j)
    {
      if (j == i)
        continue;

      derivative.multiply_assign_powred_term(this->pterms_[j]);
    }

    result.add_assign_mterm(std::move(derivative));
  }

  return result;
}

double Multiplied_Term::to_constant(void) const
{
  auto result = this->_constant;
  for (const auto& term : this->pterms_)
  {
    result *= term.to_constant();
  }

  return result;
}

std::string Multiplied_Term::to_string(void) const
{
  std::string result;

  if (this->_constant != 1.0)
  {
    result += std::to_string(_constant) + " * ";
  }

  for (const auto& pterm : this->pterms_)
  {
    result += "(" + pterm.to_string() + ") * ";
  }

  result.pop_back();
  result.pop_back();
  result.pop_back();

  return result;
}

bool Multiplied_Term::is_constant(void) const
{
  for (const auto& term : this->pterms_)
  {
    if (not term.is_constant())
      return false;
  }

  return true;
}

Symbol::Symbol(const Sym_Base& base_ptr)
{
  if (base_ptr->is_constant())
  {
    this->_constant = base_ptr->to_constant();
    return;
  }

  Powered_Term    pterm = base_ptr;
  Multiplied_Term mterm = std::move(pterm);
  this->mterms_.push_back(mterm);
}

Symbol::Symbol(Powered_Term&& pterm)
{
  if (pterm.is_constant())
  {
    this->_constant = pterm.to_constant();
  }
  else
  {
    this->mterms_.push_back(std::move(pterm));
  }
};

Sym_Base Symbol::copy(void) const
{
  return this->to_base_ptr();
}

double Symbol::operator()(const double* input) const
{
  auto result = this->_constant;

  for (const auto& term : this->mterms_)
  {
    result += term(input);
  }

  if (this->is_absolute_)
  {
    return std::abs(result);
  }
  else
  {
    return result;
  }
}

double Symbol::operator()(const std::pair<const double*, int>& input) const
{
  auto result = this->_constant;

  for (const auto& term : this->mterms_)
  {
    result += term(input);
  }

  if (this->is_absolute_)
  {
    return std::abs(result);
  }
  else
  {
    return result;
  }
}

double Symbol::operator()(const ms::math::Vector_View input) const
{
  auto result = this->_constant;

  for (const auto& term : this->mterms_)
  {
    result += term(input);
  }

  if (this->is_absolute_)
  {
    return std::abs(result);
  }
  else
  {
    return result;
  }
}

Sym_Base Symbol::get_differentiate(const int var_index) const
{
  return std::make_unique<Symbol>(this->get_diff_symbol(var_index));
}

double Symbol::to_constant(void) const { return this->_constant; }

std::string Symbol::to_string(void) const
{
  if (this->is_constant())
    return std::to_string(this->_constant);

  std::string result;

  if (this->_constant != 0.0)
  {
    result += std::to_string(this->_constant) + " + ";
  }

  for (const auto& mterm : this->mterms_)
  {
    result += mterm.to_string() + " + ";
  }

  result.pop_back();
  result.pop_back();
  result.pop_back();

  return result;
}

bool Symbol::is_constant(void) const { return this->mterms_.empty(); }

bool Symbol::is_zero(void) const
{
  for (const auto& term : this->mterms_)
  {
    if (not term.is_constant())
      return false;
  }

  return this->_constant == 0.0;
}

void Symbol::operator+=(const Symbol& other)
{
  if (other.is_constant())
  {
    this->_constant += other.to_constant();
  }
  else
  {
    this->_constant += other._constant;

    for (const auto& mterm : other.mterms_)
    {
      this->add_assign_mterm(mterm);
    }
  }
}

void Symbol::operator-=(const Symbol& other)
{
  if (other.is_constant())
  {
    this->_constant -= other.to_constant();
  }
  else
  {
    this->_constant -= other._constant;

    for (const auto& mterm : other.mterms_)
    {
      this->minus_assign_mterm(mterm);
    }
  }
}

void Symbol::operator*=(const double constant)
{
  if (constant == 0.0)
  {
    *this = Symbol();
  }
  else
  {
    this->_constant *= constant;

    for (auto& mterm : this->mterms_)
    {
      mterm *= constant;
    }
  }
}

void Symbol::pow(const double exponent)
{
  if (exponent == 1.0)
    return;

  if (this->is_constant())
  {
    this->_constant = std::pow(this->_constant, exponent);
    return;
  }

  Powered_Term temp(this->to_base_ptr(), exponent);
  (*this) = std::move(temp);
}

void Symbol::add_assign_mterm(const Multiplied_Term& mterm)
{
  if (mterm.is_constant())
  {
    this->_constant += mterm.to_constant();
  }
  else
  {
    this->mterms_.push_back(mterm);
  }
}

void Symbol::add_assign_mterm(Multiplied_Term&& mterm)
{
  if (mterm.is_constant())
  {
    this->_constant += mterm.to_constant();
  }
  else
  {
    this->mterms_.push_back(std::move(mterm));
  }
}

void Symbol::minus_assign_mterm(const Multiplied_Term& mterm)
{
  if (mterm.is_constant())
  {
    this->_constant -= mterm.to_constant();
  }
  else
  {
    auto minus_mterm = mterm * -1.0;
    this->add_assign_mterm(minus_mterm);
  }
}
Symbol Symbol::operator+(const Symbol& other) const
{
  Symbol result = *this;
  result += other;
  return result;
}

Symbol Symbol::operator-(const Symbol& other) const
{
  Symbol result = *this;
  result -= other;
  return result;
}

Symbol Symbol::operator*(const double constant) const
{
  Symbol result = *this;
  result *= constant;
  return result;
}

Symbol Symbol::operator*(const Symbol& other) const
{
  if (other.is_constant())
  {
    const auto constant = other.to_constant();

    Symbol result = *this;
    result *= constant;
    return result;
  }

  Symbol result = this->_constant * other._constant;

  const auto num_this_term = this->mterms_.size();
  const auto num_oter_term = other.mterms_.size();

  for (int i = 0; i < num_this_term; ++i)
  {
    for (int j = 0; j < num_oter_term; ++j)
    {
      result.add_assign_mterm(this->mterms_[i] * other.mterms_[j]);
    }
  }

  if (this->_constant != 0)
  {
    for (int i = 0; i < num_oter_term; ++i)
    {
      result.add_assign_mterm(this->_constant * other.mterms_[i]);
    }
  }

  if (other._constant != 0)
  {
    for (int i = 0; i < num_this_term; ++i)
    {
      result.add_assign_mterm(other._constant * this->mterms_[i]);
    }
  }

  return result;
}

Symbol Symbol::operator/(const Symbol& other) const
{
  if (other.is_constant())
  {
    const auto constant = other.to_constant();
    REQUIRE(constant != 0.0, "divider should not be zero");

    const auto divider = 1.0 / constant;
    return (*this) * divider;
  }
  else
  {
    return (*this) * other.get_pow(-1.0);
  }
}

Symbol Symbol::get_diff_symbol(const int variable_index) const
{
  REQUIRE(not this->is_absolute_, "absolute function is not a differentiable");

  Symbol result;

  for (const auto& mterm : this->mterms_)
  {
    result += mterm.get_differentiate(variable_index);
  }

  return result;
}

Symbol Symbol::get_pow(const double exponent) const
{
  auto result = *this;
  result.pow(exponent);
  return result;
}

Sym_Base Symbol::to_base_ptr(void) const
{
  return std::make_unique<Symbol>(*this);
}

Multiplied_Term operator*(const double constant, const Multiplied_Term& mterm)
{
  return mterm * constant;
}

void Symbol::be_absolute(void)
{
  this->is_absolute_ = true;
}

Polynomials get_differentiate(const Polynomials& polys, const int var_index)
{
  const auto num_polys = polys.size();

  Polynomials result(num_polys);
  for (int i = 0; i < num_polys; ++i)
  {
    result[i] = polys[i].get_diff_polynomial(var_index);
  }

  return result;
}

Symbol cal_L2_norm(const Polynomials& polys)
{
  const auto num_polys = polys.size();

  Polynomial result_poly;
  for (int i = 0; i < num_polys; ++i)
  {
    result_poly += polys[i] * polys[i];
  }

  Symbol result = result_poly.to_sym_base();
  result.pow(0.5);
  return result;
}

Symbols get_normalize(Polynomials& polys)
{
  const auto L2 = cal_L2_norm(polys);

  const auto num_polys = polys.size();

  Symbols result(num_polys);
  for (int i = 0; i < num_polys; ++i)
  {
    result[i] = Symbol(polys[i].to_sym_base()) / L2;
  }

  return result;
}

Symbols get_differentiate(const Symbols& symbols, const int var_index)
{
  const auto num_symbol = symbols.size();

  Symbols result(num_symbol);
  for (int i = 0; i < num_symbol; ++i)
  {
    result[i] = symbols[i].get_diff_symbol(var_index);
  }

  return result;
}

Symbol cal_L2_norm(const Symbols& symbols)
{
  const auto num_symbol = symbols.size();

  Symbol L2;
  for (int i = 0; i < num_symbol; ++i)
  {
    L2 += symbols[i].get_pow(2.0);
  }

  L2.pow(0.5);
  return L2;
}
Symbols get_normalize(const Symbols& symbols)
{
  const auto L2 = cal_L2_norm(symbols);

  const auto num_symbol = symbols.size();

  Symbols result(num_symbol);
  for (int i = 0; i < num_symbol; ++i)
  {
    result[i] = symbols[i] / L2;
  }

  return result;
}

} // namespace ms::sym
