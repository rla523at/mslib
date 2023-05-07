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
private:
  Reference_Geometry_Container(void) = delete;

public:
  static const Reference_Geometry& get(const Figure fig);

private:
  inline static Reference_Line ref_line_;
};

using RGeo_Container = Reference_Geometry_Container;

}

