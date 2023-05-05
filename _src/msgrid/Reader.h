#pragma once
#include <string_view>

// forward declaration
namespace ms::grid
{
struct Data;
}

namespace ms::grid
{
class Reader
{
public:
  virtual Data read(const std::string_view file_path) const = 0;
};

} // namespace ms::grid