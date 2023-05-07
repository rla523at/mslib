#include "Reference_Geometry_Container.h"

#include "Figure.h"

#include "msexception/Exception.h"

namespace ms::geo
{

const Reference_Geometry& RGeo_Container::get(const Figure fig)
{
  switch (fig)
  {
  case Figure::LINE:
    return ref_line_;

  default:
    EXCEPTION("unsupported figure");
    return ref_line_;
  }
}

} // namespace ms::geo
