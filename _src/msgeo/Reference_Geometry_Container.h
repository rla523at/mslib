#pragma once
#include "Reference_Geometry_Impl.h"

namespace ms::geo
{
enum class Figure;
}

namespace ms::geo
{

// static class
class Reference_Geometry_Container
{
public:
  static const Reference_Geometry& get(const Figure fig);

private:
  inline static Reference_Point _ref_point;
  inline static Reference_Line  _ref_line;

private:
  //Static class doesn't need to create object.
  //Thus, access to the constructor is restricted to prevent unnecessary creation.
  Reference_Geometry_Container(void) = delete;
};

} // namespace ms::geo
