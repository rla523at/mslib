#pragma once
#include <memory>
#include <string>

namespace ms::sym
{

class Base
{
public:
  virtual ~Base(void) = default;

public:
  virtual double operator()(const double* input) const = 0;

public:
  virtual std::unique_ptr<Base> copy(void) const                                  = 0;
  virtual std::unique_ptr<Base> get_differentiate(const int variable_index) const = 0;
  virtual bool                  is_constant(void) const                           = 0;
  virtual bool                  is_zero(void) const                               = 0;
  virtual double                to_constant(void) const                           = 0;
  virtual std::string           to_string(void) const                             = 0;
};

using Sym_Base = std::unique_ptr<Base>;

} // namespace ms::sym
