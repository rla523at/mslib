#include "Reference_Geometry_Container.h"

#include "Figure.h"

#include "msexception/Exception.h"

namespace ms::geo
{

const Reference_Geometry& Reference_Geometry_Container::get(const Figure fig)
{
  switch (fig)
  {
  case Figure::POINT:
    return _ref_point;
  case Figure::LINE:
    return _ref_line;
  default:
    EXCEPTION("unsupported figure");
    return _ref_point;
  }
}

} // namespace ms::geo
