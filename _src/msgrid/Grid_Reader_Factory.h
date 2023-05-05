#pragma once
#include "Reader.h"

#include <memory>

namespace ms::grid
{

class Reader_Factory
{
public:
  static std::unique_ptr<Reader> make(const std::string_view grid_type, const int dimension);
};

} // namespace ms
