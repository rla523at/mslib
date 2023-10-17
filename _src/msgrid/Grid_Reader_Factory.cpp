#include "msgrid/Grid_Reader_Factory.h"

#include "msexception/Exception.h"
#include "msgrid/Gmsh_Reader.h"
#include "msstring/string.h"

namespace ms::grid
{

std::unique_ptr<Reader> Reader_Factory::make(const std::string_view grid_type, const int dimension)
{ 
  const auto KEY = ms::string::upper_case(grid_type);

  if (ms::string::contain_icase(KEY, "GMSH"))
  {
    return std::make_unique<Gmsh_Reader>(dimension);
  }
  else
  {
    EXCEPTION(std::string(grid_type) + " is not supported grid type");
    return nullptr;
  }
}

} // namespace ms::grid