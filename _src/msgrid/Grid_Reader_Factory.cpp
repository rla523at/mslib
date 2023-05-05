#include "Grid_Reader_Factory.h"

#include "Gmsh_Reader.h"

#include "msexception/Exception.h"
#include "msstring/string.h"

namespace ms::grid
{

std::unique_ptr<Reader> Reader_Factory::make(const std::string_view grid_type, const int dimension)
{
  if (ms::string::contain_icase(grid_type, "gmsh"))
  {
    return std::make_unique<Gmsh_Reader>(dimension);
  }
  else
  {
    EXCEPTION(std::string(grid_type) + " is not supported grid type");
    return nullptr;
  }
}

} // namespace ms