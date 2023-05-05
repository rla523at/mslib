#include "Gmsh_Reader.h"

#include "Data.h"

#include "msexception/Exception.h"
#include "mspath/path.h"
#include "msstring/string.h"
#include <fstream>

namespace ms::grid
{

Element_Type str_to_element_type(std::string_view str)
{
  if (ms::string::contain_icase(str, "unspecified"))
  {
    return Element_Type::CELL;
  }
  else if (ms::string::contain_icase(str, "slip", "wall"))
  {
    return Element_Type::SLIP_WALL;
  }
  else if (ms::string::contain_icase(str, "superSonic", "inlet"))
  {
    return Element_Type::SUPERSONIC_INLET;
  }
  else if (ms::string::contain_icase(str, "SuperSonic", "outlet"))
  {
    return Element_Type::SUPERSONIC_OUTLET;
  }
  else if (ms::string::contain_icase(str, "periodic"))
  {
    return Element_Type::PERIODIC;
  }
  else if (ms::string::contain_icase(str, "initial", "constant"))
  {
    return Element_Type::INTIAL_CONSTANT;
  }
  else
  {
    EXCEPTION(std::string(str) + " is not supproted element type");
    return Element_Type::NOT_IN_LIST;
  }
}

enum class Gmsh_Figure_Type
{
  POINT   = 0,
  LINE_P1 = 1,
  LINE_P2 = 8,
  LINE_P3 = 26,
  LINE_P4 = 27,
  LINE_P5 = 28,
  LINE_P6 = 62,
  TRIS_P1 = 2,
  TRIS_P2 = 9,
  TRIS_P3 = 21,
  TRIS_P4 = 23,
  TRIS_P5 = 25,
  QUAD_P1 = 3,
  QUAD_P2 = 10,
  QUAD_P3 = 36,
  QUAD_P4 = 37,
  QUAD_P5 = 38,
  QUAD_P6 = 47,
  TETS_P1 = 4,
  TETS_P2 = 11,
  TETS_P3 = 29,
  TETS_P4 = 30,
  TETS_P5 = 31,
  HEXA_P1 = 5,
  HEXA_P2 = 12,
  HEXA_P3 = 92,
  HEXA_P4 = 93,
  HEXA_P5 = 94,
  PRIS_P1 = 6,
  PRIS_P2 = 13,
  PRIS_P3 = 90,
  PRIS_P4 = 91,
  PRIS_P5 = 106,
  PYRA_P1 = 7,
  PYRA_P2 = 14,
  PYRA_P3 = 118,
  PYRA_P4 = 119
};

Figure index_to_figure_type(const int figure_type_index)
{
  switch (static_cast<Gmsh_Figure_Type>(figure_type_index))
  {
  case Gmsh_Figure_Type::POINT: return Figure::POINT;
  case Gmsh_Figure_Type::LINE_P1:
  case Gmsh_Figure_Type::LINE_P2:
  case Gmsh_Figure_Type::LINE_P3:
  case Gmsh_Figure_Type::LINE_P4:
  case Gmsh_Figure_Type::LINE_P5:
  case Gmsh_Figure_Type::LINE_P6: return Figure::LINE;
  case Gmsh_Figure_Type::TRIS_P1:
  case Gmsh_Figure_Type::TRIS_P2:
  case Gmsh_Figure_Type::TRIS_P3:
  case Gmsh_Figure_Type::TRIS_P4:
  case Gmsh_Figure_Type::TRIS_P5: return Figure::TRIANGLE;
  case Gmsh_Figure_Type::QUAD_P1:
  case Gmsh_Figure_Type::QUAD_P2:
  case Gmsh_Figure_Type::QUAD_P3:
  case Gmsh_Figure_Type::QUAD_P4:
  case Gmsh_Figure_Type::QUAD_P5:
  case Gmsh_Figure_Type::QUAD_P6: return Figure::QUADRILATERAL;
  case Gmsh_Figure_Type::TETS_P1:
  case Gmsh_Figure_Type::TETS_P2:
  case Gmsh_Figure_Type::TETS_P3:
  case Gmsh_Figure_Type::TETS_P4:
  case Gmsh_Figure_Type::TETS_P5: return Figure::TETRAHEDRAL;
  case Gmsh_Figure_Type::HEXA_P1:
  case Gmsh_Figure_Type::HEXA_P2:
  case Gmsh_Figure_Type::HEXA_P3:
  case Gmsh_Figure_Type::HEXA_P4:
  case Gmsh_Figure_Type::HEXA_P5: return Figure::HEXAHEDRAL;
  case Gmsh_Figure_Type::PRIS_P1:
  case Gmsh_Figure_Type::PRIS_P2:
  case Gmsh_Figure_Type::PRIS_P3:
  case Gmsh_Figure_Type::PRIS_P4:
  case Gmsh_Figure_Type::PRIS_P5: return Figure::PRISM;
  case Gmsh_Figure_Type::PYRA_P1:
  case Gmsh_Figure_Type::PYRA_P2:
  case Gmsh_Figure_Type::PYRA_P3:
  case Gmsh_Figure_Type::PYRA_P4: return Figure::PYRAMID;
  default:
    EXCEPTION("invalid element type index");
    return Figure::NOT_IN_LIST;
  }
}

short index_to_figure_order(const int figure_type_index)
{
  switch (static_cast<Gmsh_Figure_Type>(figure_type_index))
  {
  case Gmsh_Figure_Type::POINT: return 0;
  case Gmsh_Figure_Type::LINE_P1:
  case Gmsh_Figure_Type::TRIS_P1:
  case Gmsh_Figure_Type::QUAD_P1:
  case Gmsh_Figure_Type::TETS_P1:
  case Gmsh_Figure_Type::HEXA_P1:
  case Gmsh_Figure_Type::PRIS_P1:
  case Gmsh_Figure_Type::PYRA_P1: return 1;
  case Gmsh_Figure_Type::LINE_P2:
  case Gmsh_Figure_Type::TRIS_P2:
  case Gmsh_Figure_Type::QUAD_P2:
  case Gmsh_Figure_Type::TETS_P2:
  case Gmsh_Figure_Type::HEXA_P2:
  case Gmsh_Figure_Type::PRIS_P2:
  case Gmsh_Figure_Type::PYRA_P2: return 2;
  case Gmsh_Figure_Type::LINE_P3:
  case Gmsh_Figure_Type::TRIS_P3:
  case Gmsh_Figure_Type::QUAD_P3:
  case Gmsh_Figure_Type::TETS_P3:
  case Gmsh_Figure_Type::HEXA_P3:
  case Gmsh_Figure_Type::PRIS_P3:
  case Gmsh_Figure_Type::PYRA_P3: return 3;
  case Gmsh_Figure_Type::LINE_P4:
  case Gmsh_Figure_Type::TRIS_P4:
  case Gmsh_Figure_Type::QUAD_P4:
  case Gmsh_Figure_Type::TETS_P4:
  case Gmsh_Figure_Type::HEXA_P4:
  case Gmsh_Figure_Type::PRIS_P4:
  case Gmsh_Figure_Type::PYRA_P4: return 4;
  case Gmsh_Figure_Type::LINE_P5:
  case Gmsh_Figure_Type::TRIS_P5:
  case Gmsh_Figure_Type::QUAD_P5:
  case Gmsh_Figure_Type::TETS_P5:
  case Gmsh_Figure_Type::HEXA_P5:
  case Gmsh_Figure_Type::PRIS_P5: return 5;
  case Gmsh_Figure_Type::LINE_P6:
  case Gmsh_Figure_Type::QUAD_P6: return 6;
  default:
    EXCEPTION("invalid element type index");
    return -1;
  }
}

Data Gmsh_Reader::read(const std::string_view file_path) const
{
  REQUIRE(path::is_exist_file(file_path), "grid file should be exist");

  std::ifstream file(file_path.data());
  REQUIRE(file.is_open(), "grid file should be open");

  auto node_data    = this->extract_block_type_node_data(file);
  auto element_data = this->extract_element_datas(file);

  return {node_data, element_data};
}

Nodes_Data Gmsh_Reader::extract_block_type_node_data(std::ifstream& file) const
{
  constexpr auto delimiter = ' ';

  std::string str;
  while (std::getline(file, str))
  {
    if (str == "$Nodes") break;
  }
  REQUIRE(str == "$Nodes", "fail to file node data in grid file");

  std::getline(file, str);
  const auto num_nodes = ms::string::str_to_value<int>(str);

  std::vector<double> coordinates(num_nodes * this->_dimension);

  for (int i = 0; i < num_nodes; ++i)
  {
    std::getline(file, str);

    const auto parsed_strs = ms::string::parse_by(str, delimiter);

    // save coordinates as block style
    for (int j = 0; j < this->_dimension; ++j)
    {
      auto index = j * num_nodes + i;

      // parsed_strs[0] refers to node index
      coordinates[index] = ms::string::str_to_value<double>(parsed_strs[j + 1]);
    }
  }

  Nodes_Data node_data;
  node_data.type        = Coordinate_Type::BLOCK;
  node_data.dimension   = this->_dimension;
  node_data.num_nodes   = num_nodes;
  node_data.coordinates = std::move(coordinates);

  return node_data;
}

std::vector<Element_Data> Gmsh_Reader::extract_element_datas(std::ifstream& file) const
{
  constexpr auto delimiter        = ' ';
  constexpr auto num_type_indexes = 5;

  std::string str;
  while (std::getline(file, str))
  {
    if (str == "$Elements") break;
  }
  REQUIRE(str == "$Elements", "fail to file element data in grid file");

  std::getline(file, str);
  const auto num_elements = ms::string::str_to_value<int>(str);

  std::vector<int>              figure_type_indexes(num_elements);
  std::vector<int>              physical_group_indexes(num_elements);
  std::vector<std::vector<int>> consisting_node_indexess(num_elements);

  for (int i = 0; i < num_elements; ++i)
  {
    std::getline(file, str);

    const auto parsed_strs = ms::string::parse_by(str, delimiter);
    // parsed_strs[0]  : element index
    // parsed_strs[1]  : figure type index
    // parsed_strs[2]  : tag index
    // parsed_strs[3]  : physical group index
    // parsed_strs[4]  : element group index
    // parsed_strs[5-] : consisting node indexes

    figure_type_indexes[i]    = ms::string::str_to_value<int>(parsed_strs[1]);
    physical_group_indexes[i] = ms::string::str_to_value<int>(parsed_strs[3]);

    const auto num_consisting_nodes = parsed_strs.size() - num_type_indexes;

    auto& indexes = consisting_node_indexess[i];
    indexes.resize(num_consisting_nodes);

    for (int j = 0; j < num_consisting_nodes; ++j)
    {
      indexes[j] = ms::string::str_to_value<int>(parsed_strs[num_type_indexes + j]);
    }
  }

  const auto physical_group_index_to_element_type = this->extract_physical_group_index_to_element_type(file);

  std::vector<Element_Data> element_datas(num_elements);

  for (int i = 0; i < num_elements; ++i)
  {
    auto& data = element_datas[i];

    data.element_type = physical_group_index_to_element_type.at(physical_group_indexes[i]);
    data.figure       = ms::grid::index_to_figure_type(figure_type_indexes[i]);
    data.figure_order = ms::grid::index_to_figure_order(figure_type_indexes[i]);
    data.node_indexes = std::move(consisting_node_indexess[i]);
  }

  return element_datas;
}

std::map<int, Element_Type> Gmsh_Reader::extract_physical_group_index_to_element_type(std::ifstream& file) const
{
  constexpr auto delimiter = ' ';

  std::string str;
  while (std::getline(file, str))
  {
    if (str == "$PhysicalNames") break;
  }
  REQUIRE(str == "$PhysicalNames", "fail to file physical name data in grid file");

  std::getline(file, str);
  const auto num_physical_groups = ms::string::str_to_value<int>(str);

  std::map<int, Element_Type> physical_group_index_to_element_type;

  for (int i = 0; i < num_physical_groups; ++i)
  {
    std::getline(file, str);

    const auto parsed_strs = ms::string::parse_by(str, delimiter);
    // parsed_strs[0] : space_dimension
    // parsed_strs[1] : physical group index
    // parsed_strs[2] : "name"

    const auto physical_group_index = ms::string::str_to_value<int>(parsed_strs[1]);
    const auto name                 = ms::string::remove(parsed_strs[2], '"');

    const auto element_type = ms::grid::str_to_element_type(name);
    physical_group_index_to_element_type.emplace(physical_group_index, element_type);
  }

  return physical_group_index_to_element_type;
}

} // namespace ms::grid::reader