#include "msgrid/Gmsh_Reader.h"

#include "msexception/Exception.h"
#include "msfilesystem/filesystem.h"
#include "msgrid/Data.h"
#include "msstring/string.h"
#include <fstream>
#include <unordered_map>

namespace ms::grid
{

struct Gmsh_Nodes_Data
{
  std::vector<int>    numbers;
  std::vector<double> coordinates;
};

struct Gmsh_Elements_Data
{
  std::vector<int>              elem_numbers;
  std::vector<int>              figure_type_numbers;
  std::vector<int>              physical_group_numbers;
  std::vector<std::vector<int>> node_numbers_s;
};

struct Gmsh_Physical_Data
{
  std::map<int, std::string> physical_group_number_to_name;
};

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
  POINT   = 15,
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

Figure index_to_figure_type(const int figure_type_numbers)
{
  switch (static_cast<Gmsh_Figure_Type>(figure_type_numbers))
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

} // namespace ms::grid

/*




*/

namespace ms::grid
{

Grid_Data Gmsh_Reader::read(const std::string_view file_path) const
{
  REQUIRE(filesystem::is_exist_file(file_path), "grid file should be exist");

  std::ifstream file(file_path.data());
  REQUIRE(file.is_open(), "grid file should be open");

  Gmsh_Nodes_Data    node_data;
  Gmsh_Elements_Data element_data;
  Gmsh_Physical_Data physical_data;

  bool is_not_read_node_data = true;
  bool is_not_read_elem_data = true;
  bool is_not_read_phys_data = true;

  std::string str;
  while (std::getline(file, str))
  {
    if (is_not_read_node_data)
    {
      if (str == "$Nodes")
      {
        this->read_node_data(file, node_data);
        is_not_read_node_data = false;
      }
    }

    if (is_not_read_elem_data)
    {
      if (str == "$Elements")
      {
        this->read_elem_data(file, element_data);
        is_not_read_elem_data = false;
      }
    }

    if (is_not_read_phys_data)
    {
      if (str == "$PhysicalNames")
      {
        this->read_phys_data(file, physical_data);
        is_not_read_phys_data = false;
      }
    }
  }

  return this->convert(std::move(node_data), std::move(element_data), physical_data);
}

Grid_Data Gmsh_Reader::convert(Gmsh_Nodes_Data&& node_data, Gmsh_Elements_Data&& elem_data, const Gmsh_Physical_Data& phys_data) const
{
  constexpr auto comma = ',';

  Grid_Data grid_data;

  // convert Gmsh_Node_Data to Grid_Nodes_Data
  const auto num_nodes = node_data.numbers.size();

  auto& grid_nodes_data   = grid_data.nodes_data;
  grid_nodes_data.nodes   = ms::geo::Nodes(num_nodes, this->_dimension, std::move(node_data.coordinates));
  grid_nodes_data.numbers = std::move(node_data.numbers);

  // convert Gmsh_Element_Data & Gmsh_Physical_Data to Grid_Element_Data & Grid_Peridoic_Data
  auto& grid_element_datas  = grid_data.element_datas;
  auto& grid_periodic_datas = grid_data.periodic_datas;

  const auto& elem_numbers                  = elem_data.elem_numbers;
  const auto& figure_type_numbers           = elem_data.figure_type_numbers;
  auto&       node_numbers_s                = elem_data.node_numbers_s;
  const auto& physical_group_numbers        = elem_data.physical_group_numbers;
  const auto& physical_group_number_to_name = phys_data.physical_group_number_to_name;

  const auto num_elements = elem_numbers.size();
  grid_element_datas.reserve(num_elements);
  grid_periodic_datas.reserve(num_elements);

  for (int i = 0; i < num_elements; ++i)
  {
    const auto  elem_number           = elem_numbers[i];
    const auto  physical_group_number = physical_group_numbers[i];
    const auto& name                  = physical_group_number_to_name.at(physical_group_number);
    const auto  elem_type             = str_to_element_type(name);

    Grid_Element_Data elem_data;
    elem_data.number       = elem_number;
    elem_data.type         = elem_type;
    elem_data.figure       = ms::grid::index_to_figure_type(figure_type_numbers[i]);
    elem_data.node_numbers = std::move(node_numbers_s[i]);

    // 이 부분 테스트하기!
    //// The node numbers in the Grid_Element_Data start from 0.
    //// Since Gmsh's node numbering starts from 1, subtract 1 from it.
    // for (auto& node_number : elem_data.node_numbers)
    //{
    //   node_number--;
    // }

    if (elem_type == Element_Type::PERIODIC)
    {
      const auto parsed_strs = ms::string::parse_by(name, comma);

      const auto          num_periodic_direction_components = parsed_strs.size() - 1;
      std::vector<double> periodic_direction(num_periodic_direction_components);

      for (auto j = 0; j < num_periodic_direction_components; ++j)
      {
        periodic_direction[j] = ms::string::str_to_value<double>(parsed_strs[j]);
      }

      Grid_Peridoic_Data periodic_data = {std::move(elem_data), std::move(periodic_direction)};
      grid_periodic_datas.push_back(std::move(periodic_data));
    }
    else
    {
      grid_element_datas.push_back(std::move(elem_data));
    }
  }

  grid_element_datas.shrink_to_fit();
  grid_periodic_datas.shrink_to_fit();

  return grid_data;
}

void Gmsh_Reader::read_node_data(std::ifstream& file, Gmsh_Nodes_Data& data) const
{
  constexpr auto delimiter = ' ';

  auto& numbers     = data.numbers;
  auto& coordinates = data.coordinates;

  std::string str;
  std::getline(file, str);

  const auto num_nodes = ms::string::str_to_value<int>(str);
  numbers.resize(num_nodes);
  coordinates.resize(num_nodes * this->_dimension);

  for (int i = 0; i < num_nodes; ++i)
  {
    std::getline(file, str);

    const auto parsed_strs = ms::string::parse_by(str, delimiter);
    // parsed_strs[0]  : node number
    // parsed_strs[1]  : x coordinate
    // parsed_strs[2]  : y coordinate
    // parsed_strs[3]  : z coordinate

    numbers[i] = ms::string::str_to_value<int>(parsed_strs[0]);

    auto coord_start_index = i * this->_dimension;
    for (int j = 0; j < this->_dimension; ++j)
    {
      coordinates[coord_start_index + j] = ms::string::str_to_value<double>(parsed_strs[j + 1]);
    }
  }
}

void Gmsh_Reader::read_elem_data(std::ifstream& file, Gmsh_Elements_Data& data) const
{
  constexpr auto space            = ' ';
  constexpr auto num_type_indexes = 5;

  auto& elem_numbers           = data.elem_numbers;
  auto& figure_type_numbers    = data.figure_type_numbers;
  auto& node_numbers_s         = data.node_numbers_s;
  auto& physical_group_numbers = data.physical_group_numbers;

  std::string str;
  std::getline(file, str);
  const auto num_elements = ms::string::str_to_value<int>(str);

  elem_numbers.resize(num_elements);
  figure_type_numbers.resize(num_elements);
  physical_group_numbers.resize(num_elements);
  node_numbers_s.resize(num_elements);

  for (int i = 0; i < num_elements; ++i)
  {
    std::getline(file, str);

    const auto parsed_strs = ms::string::parse_by(str, space);
    // parsed_strs[0]  : element number
    // parsed_strs[1]  : figure type number
    // parsed_strs[2]  : tag number
    // parsed_strs[3]  : physical group number
    // parsed_strs[4]  : element group number
    // parsed_strs[5-] : node numbers composing element

    elem_numbers[i]           = ms::string::str_to_value<int>(parsed_strs[0]);
    figure_type_numbers[i]    = ms::string::str_to_value<int>(parsed_strs[1]);
    physical_group_numbers[i] = ms::string::str_to_value<int>(parsed_strs[3]);

    const auto num_nodes = parsed_strs.size() - num_type_indexes;

    auto& node_numbers = node_numbers_s[i];
    node_numbers.resize(num_nodes);

    for (int j = 0; j < num_nodes; ++j)
    {
      node_numbers[j] = ms::string::str_to_value<int>(parsed_strs[num_type_indexes + j]);
    }
  }
}

void Gmsh_Reader::read_phys_data(std::ifstream& file, Gmsh_Physical_Data& data) const
{
  constexpr auto delimiter = ' ';

  std::string str;
  std::getline(file, str);
  const auto num_physical_groups = ms::string::str_to_value<int>(str);

  for (int i = 0; i < num_physical_groups; ++i)
  {
    std::getline(file, str);

    const auto parsed_strs = ms::string::parse_by(str, delimiter);
    // parsed_strs[0] : space_dimension
    // parsed_strs[1] : physical group index
    // parsed_strs[2] : "name"

    const auto physical_group_numbers = ms::string::str_to_value<int>(parsed_strs[1]);
    const auto name                   = ms::string::remove(parsed_strs[2], '"');

    data.physical_group_number_to_name.emplace(physical_group_numbers, name);
  }
}

// void Gmsh_Reader::extract_block_type_node_data(std::ifstream& file, Grid_Data& data) const
//{
//   constexpr auto delimiter = ' ';
//
//   std::string str;
//   std::getline(file, str);
//   const auto num_nodes = ms::string::str_to_value<int>(str);
//
//   std::vector<double> coordinates(num_nodes * this->_dimension);
//
//   for (int i = 0; i < num_nodes; ++i)
//   {
//     std::getline(file, str);
//
//     const auto parsed_strs = ms::string::parse_by(str, delimiter);
//     // parsed_strs[0]  : node number
//     // parsed_strs[1]  : x coordinate
//     // parsed_strs[2]  : y coordinate
//     // parsed_strs[3]  : z coordinate
//
//     // save coordinates as block style
//     for (int j = 0; j < this->_dimension; ++j)
//     {
//       auto index = j * num_nodes + i;
//
//       coordinates[index] = ms::string::str_to_value<double>(parsed_strs[j + 1]);
//     }
//   }
//
//   auto& nodes_data       = data.nodes_data;
//   nodes_data.type        = Coordinate_Type::BLOCK;
//   nodes_data.dimension   = this->_dimension;
//   nodes_data.num_nodes   = num_nodes;
//   nodes_data.coordinates = std::move(coordinates);
// }
//
// void Gmsh_Reader::extract_element_datas(std::ifstream& file, Grid_Data& data) const
//{
//   constexpr auto space            = ' ';
//   constexpr auto comma            = ',';
//   constexpr auto num_type_indexes = 5;
//
//   std::string str;
//   while (std::getline(file, str))
//   {
//     if (str == "$Elements") break;
//   }
//   REQUIRE(str == "$Elements", "fail to file element grid_data in grid file");
//
//   std::getline(file, str);
//   const auto num_elements = ms::string::str_to_value<int>(str);
//
//   std::vector<int>              elem_index_to_figure_type_index(num_elements);
//   std::vector<int>              elem_index_to_physical_group_index(num_elements);
//   std::vector<std::vector<int>> elem_index_to_node_numbers(num_elements);
//
//   for (int i = 0; i < num_elements; ++i)
//   {
//     std::getline(file, str);
//
//     const auto parsed_strs = ms::string::parse_by(str, space);
//     // parsed_strs[0]  : element index
//     // parsed_strs[1]  : figure type index
//     // parsed_strs[2]  : tag index
//     // parsed_strs[3]  : physical group index
//     // parsed_strs[4]  : element group index
//     // parsed_strs[5-] : node numbers composing element
//
//     elem_index_to_figure_type_index[i]    = ms::string::str_to_value<int>(parsed_strs[1]);
//     elem_index_to_physical_group_index[i] = ms::string::str_to_value<int>(parsed_strs[3]);
//
//     const auto num_nodes = parsed_strs.size() - num_type_indexes;
//
//     auto& numbers = elem_index_to_node_numbers[i];
//     numbers.resize(num_nodes);
//
//     for (int j = 0; j < num_nodes; ++j)
//     {
//       // The node numbers in the Element_Data must follow the rule starting from 0.
//       // Since Gmsh's node numbering starts from 1, subtract 1 from it.
//       numbers[j] = ms::string::str_to_value<int>(parsed_strs[num_type_indexes + j]) - 1;
//     }
//   }
//
//   const auto physical_group_index_to_name = this->extract_physical_group_index_to_name(file);
//
//   auto& element_datas  = data.element_datas;
//   auto& periodic_datas = data.periodic_datas;
//   element_datas.reserve(num_elements);
//   periodic_datas.reserve(num_elements);
//
//   for (int i = 0; i < num_elements; ++i)
//   {
//     const auto  physical_group_index = elem_index_to_physical_group_index[i];
//     const auto& name                 = physical_group_index_to_name.at(physical_group_index);
//     const auto  element_type         = str_to_element_type(name);
//
//     Element_Data elem_data;
//     elem_data.element_type = element_type;
//     elem_data.figure       = ms::grid::index_to_figure_type(elem_index_to_figure_type_index[i]);
//     elem_data.node_numbers = std::move(elem_index_to_node_numbers[i]);
//
//     if (element_type == Element_Type::PERIODIC)
//     {
//       const auto parsed_strs = ms::string::parse_by(name, comma);
//
//       const auto          num_periodic_direction_components = parsed_strs.size() - 1;
//       std::vector<double> periodic_direction(num_periodic_direction_components);
//
//       for (auto j = 0; j < num_periodic_direction_components; ++j)
//       {
//         periodic_direction[j] = ms::string::str_to_value<double>(parsed_strs[j]);
//       }
//
//       periodic_datas.push_back({std::move(elem_data), std::move(periodic_direction)});
//     }
//     else
//     {
//       element_datas.push_back(std::move(elem_data));
//     }
//   }
//
//   element_datas.shrink_to_fit();
//   periodic_datas.shrink_to_fit();
// }
//
// std::map<int, std::string> Gmsh_Reader::extract_physical_group_index_to_name(std::ifstream& file) const
//{
//   constexpr auto delimiter = ' ';
//
//   std::string str;
//   while (std::getline(file, str))
//   {
//     if (str == "$PhysicalNames") break;
//   }
//   REQUIRE(str == "$PhysicalNames", "fail to file physical name grid_data in grid file");
//
//   std::getline(file, str);
//   const auto num_physical_groups = ms::string::str_to_value<int>(str);
//
//   std::map<int, std::string> physical_group_index_to_name;
//
//   for (int i = 0; i < num_physical_groups; ++i)
//   {
//     std::getline(file, str);
//
//     const auto parsed_strs = ms::string::parse_by(str, delimiter);
//     // parsed_strs[0] : space_dimension
//     // parsed_strs[1] : physical group index
//     // parsed_strs[2] : "name"
//
//     const auto physical_group_index = ms::string::str_to_value<int>(parsed_strs[1]);
//     const auto name                 = ms::string::remove(parsed_strs[2], '"');
//
//     physical_group_index_to_name.emplace(physical_group_index, name);
//   }
//
//   return physical_group_index_to_name;
// }

} // namespace ms::grid