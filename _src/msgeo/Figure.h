#pragma once

namespace ms::geo
{
enum class Figure
{
  NOT_FIGURE   = -1,
  POINT         = 0,
  LINE          = 1,
  TRIANGLE      = 2,
  QUADRILATERAL = 3,
  TETRAHEDRAL   = 4,
  HEXAHEDRAL    = 5,
  PRISM         = 6,
  PYRAMID       = 7,
  NUM_FIGURES,
};
}