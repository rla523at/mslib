#pragma once
#include <string_view>

// forward declaration
namespace ms::grid
{
struct Grid_Data;
}

namespace ms::grid
{
class Reader
{
public:
  virtual ~Reader() = default;

public:
  virtual Grid_Data read(const std::string_view file_path) const = 0;
};

} // namespace ms::grid