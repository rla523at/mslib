#include "mstecplot/Writer.h"

#include "msexception/Exception.h"
#include "msfilesystem/filesystem.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>

namespace ms::tecplot
{

void Writer::set_configuration(Configuration&& configuration)
{
  REQUIRE(ms::filesystem::is_folder_path(configuration.folder_path), "folder path in configuration file is not valid");

  THIS::_configuration                   = std::move(configuration);
  THIS::_is_configuration_setup_complete = true;

  THIS::_dimension          = 0;
  THIS::_num_nodes          = 0;
  THIS::_num_elements       = 0;
  THIS::_is_grid_data_saved = false;
}

void Writer::write_grid_file(const Grid_Data& data)
{
  REQUIRE(THIS::_is_configuration_setup_complete, "Configuration has not been setup.");

  THIS::save_part_of_grid_data(data);

  const auto grid_folder_path = THIS::_configuration.folder_path + THIS::_configuration.title + "/";
  if (!ms::filesystem::is_exist_folder(grid_folder_path))
  {
    std::cout << "msTecplot :: Given folder path is not exist! \n";
    ms::filesystem::make_folder(grid_folder_path);
    std::cout << "msTecplot :: " << grid_folder_path << " is created! \n";
  }

  const auto grid_file_path = grid_folder_path + "/_grid.plt";

  if (THIS::_configuration.file_format == File_Format::ASCII)
  {
    std::ofstream out_file(grid_file_path);
    REQUIRE(out_file.is_open(), "fail to open grid file.");

    THIS::write_ASCII_grid_file_header(out_file);
    THIS::write_ASCII_grid_file_data(out_file, data.coordinates_ptr);
    THIS::write_ASCII_connectivity(out_file, data.connectivities);
    out_file.close();
  }
  else if (THIS::_configuration.file_format == File_Format::Binary)
  {
    Binary_File out_file(grid_file_path);
    THIS::write_binary_grid_file_header(out_file);
    THIS::write_binary_grid_file_data(out_file, data.coordinates_ptr);
    THIS::write_binary_connectivity(out_file, data.connectivities);
    out_file.close();
  }
  else
  {
    EXCEPTION("not supported file format");
  }

  THIS::_strand_id = 0;
}

void Writer::write_solution_file(const Solution_Data& info)
{
  REQUIRE(THIS::_is_configuration_setup_complete, "Configuration has not been setup.");
  REQUIRE(THIS::_is_grid_data_saved, "Grid file data is not saved. Please call the write grid file function first.");

  const auto solution_folder_path = THIS::_configuration.folder_path + THIS::_configuration.title + "/";
  if (!ms::filesystem::is_exist_folder(solution_folder_path))
  {
    std::cout << "msTecplot :: Given folder path is not exist! \n";
    ms::filesystem::make_folder(solution_folder_path);
    std::cout << "msTecplot :: " << solution_folder_path << " is created! \n";
  }

  const auto solution_file_path = solution_folder_path + "solution" + std::to_string(THIS::_strand_id) + ".plt";

  if (THIS::_configuration.file_format == File_Format::ASCII)
  {
    std::ofstream out_file(solution_file_path);
    REQUIRE(out_file.is_open(), "fail to open solution file.");

    THIS::write_ASCII_solution_file_header(out_file, info);
    THIS::write_ASCII_solution_file_data(out_file, info);
    out_file.close();
  }
  else if (THIS::_configuration.file_format == File_Format::Binary)
  {
    Binary_File out_file(solution_file_path);
    THIS::write_binary_solution_file_header(out_file, info);
    THIS::write_binary_solution_file_data(out_file, info);
    out_file.close();
  }
  else
  {
    EXCEPTION("not supported file format");
  }

  ++THIS::_strand_id;
}

void Writer::write_ASCII_connectivity(std::ofstream& out_file, const std::vector<std::vector<int>>& connectivities)
{
  out_file << "# CONNECTIVITY \n";
  for (const auto& connectivity : connectivities)
  {
    for (const auto node_index : connectivity)
    {
      out_file << node_index << " ";
    }
    out_file << "\n";
  }
}

void Writer::write_ASCII_grid_file_header(std::ofstream& out_file)
{
  REQUIRE(THIS::_is_grid_data_saved, "Grid file data should be saved first");

  const auto zone_title = THIS::_configuration.title + "_GRID";

  out_file << "# FILE HEADER \n";
  out_file << "TITLE = " << THIS::_configuration.title << "\n";
  out_file << "FILETYPE = GRID\n";
  out_file << "VARIABLES = " << THIS::make_ASCII_grid_variable_string() << "\n\n";

  out_file << "# ZONE HEADER \n";
  out_file << "ZONE \n";
  out_file << "T = " + zone_title << "\n";
  out_file << "ZONETYPE =  " + THIS::zone_type_str() << "\n";
  out_file << "NODES = " + std::to_string(THIS::_num_nodes) << "\n";
  out_file << "ELEMENTS = " + std::to_string(THIS::_num_elements) << "\n";
  out_file << "DATAPACKING = BLOCK \n\n";
}

void Writer::write_ASCII_grid_file_data(std::ofstream& out_file, const double* coordinates_ptr)
{
  constexpr auto num_decimical_points = 6;
  constexpr auto data_space           = 14;
  constexpr auto data_in_line         = 10;

  REQUIRE(THIS::_is_grid_data_saved, "Grid file data should be saved first");

  out_file << "# VARIABLE DATA \n";
  out_file << std::scientific << std::setprecision(num_decimical_points) << std::left;

  for (int i = 0; i < THIS::_dimension; ++i)
  {
    int index = i * THIS::_num_nodes;
    int count = 1;

    for (int j = 0; j < THIS::_num_nodes; ++j, ++count)
    {
      out_file << std::setw(data_space) << coordinates_ptr[index + j];

      if (data_in_line == count)
      {
        count = 0;
        out_file << "\n";
      }
    }

    out_file << "\n\n";
  }
}

void Writer::write_ASCII_solution_file_header(std::ofstream& out_file, const Solution_Data& info)
{
  REQUIRE(THIS::_is_grid_data_saved, "Grid file data should be saved first");

  const auto straind_id = std::to_string(THIS::_strand_id);
  const auto zone_title = THIS::_configuration.title + "_SOL" + straind_id;

  out_file << "# FILE HEADER \n";
  out_file << "TITLE = " << THIS::_configuration.title << "\n";
  out_file << "FILETYPE = SOLUTION\n";
  out_file << "VARIABLES = " << THIS::make_ASCII_variable_string(info.variable_strs) << "\n\n";

  out_file << "# ZONE HEADER \n";
  out_file << "ZONE \n";
  out_file << "T = " + zone_title << "\n";
  out_file << "ZONETYPE =  " + THIS::zone_type_str() << "\n";
  out_file << "NODES = " + std::to_string(THIS::_num_nodes) << "\n";
  out_file << "ELEMENTS = " + std::to_string(THIS::_num_elements) << "\n";
  out_file << "DATAPACKING = BLOCK \n";
  out_file << "VARLOCATION = " + info.var_location_str() << "\n";
  out_file << "STRANDID = " + straind_id << "\n";
  out_file << "SOLUTIONTIME = " + std::to_string(info.solution_time) << "\n\n";
}

void Writer::write_ASCII_solution_file_data(std::ofstream& out_file, const Solution_Data& info)
{
  constexpr auto decimical_point = 10;
  constexpr auto data_width      = 18;
  constexpr auto data_in_line    = 10;

  out_file << "# VARIABLE DATA \n";
  out_file << std::scientific << std::setprecision(decimical_point) << std::left;

  auto index = 0;
  for (int i = 0; i < info.num_variables(); ++i)
  {
    const auto num_datas = THIS::num_solution_datas(info.var_locations_ptr[i]);

    int count = 1;

    for (int j = 0; j < num_datas; ++j, ++count, ++index)
    {
      out_file << std::setw(data_width) << info.values_ptr[index];

      if (data_in_line == count)
      {
        count = 0;
        out_file << "\n";
      }
    }

    out_file << "\n\n";
  }
}

std::string Writer::make_ASCII_variable_string(const std::vector<std::string>& variable_strs)
{
  std::string result;

  for (const auto& variable_str : variable_strs)
  {
    result += variable_str + ",";
  }

  result.pop_back();

  return result;
}

std::string Writer::make_ASCII_grid_variable_string(void)
{
  switch (THIS::_dimension)
  {
  case 1: return "X";
  case 2: return "X,Y";
  case 3: return "X,Y,Z";
  default: EXCEPTION("unsupported dimension"); return "";
  }
}

std::vector<int> Writer::make_binary_grid_variable_string(void)
{
  std::vector<int> binary_grid_name;

  const auto x = THIS::to_tecplot_binary_format("X");
  binary_grid_name.insert(binary_grid_name.end(), x.begin(), x.end());

  if (2 <= THIS::_dimension)
  {
    const auto y = THIS::to_tecplot_binary_format("Y");
    binary_grid_name.insert(binary_grid_name.end(), y.begin(), y.end());
  }

  if (3 == THIS::_dimension)
  {
    const auto z = THIS::to_tecplot_binary_format("Z");
    binary_grid_name.insert(binary_grid_name.end(), z.begin(), z.end());
  }

  return binary_grid_name;
}

std::vector<int> Writer::make_binary_variable_string(const std::vector<std::string>& variable_strs)
{
  std::vector<int> binary_variable_str;

  for (const auto& str : variable_strs)
  {
    const auto bin_str = THIS::to_tecplot_binary_format(str);
    binary_variable_str.insert(binary_variable_str.end(), bin_str.begin(), bin_str.end());
  }

  return binary_variable_str;
}

int Writer::num_solution_datas(const Variable_Location variable_location)
{
  if (variable_location == Variable_Location::NODE)
  {
    return THIS::_num_nodes;
  }
  else
  {
    return THIS::_num_elements;
  }
}

void Writer::save_part_of_grid_data(const Grid_Data& info)
{
  THIS::_dimension          = info.dimension;
  THIS::_num_nodes          = info.num_nodes;
  THIS::_num_elements       = info.num_elements;
  THIS::_is_grid_data_saved = true;
}

std::string Writer::zone_type_str(void)
{
  switch (THIS::_dimension)
  {
  case 1: return "FELINESEG";
  case 2: return "FETRIANGLE";
  default: EXCEPTION("unsupported diemsnion"); return "";
  }
}

void Writer::write_binary_grid_file_header(Binary_File& out_file)
{
  REQUIRE(THIS::_is_grid_data_saved, "Grid file data should be saved first");

  const auto zone_title = THIS::_configuration.title + "_GRID";

  // Setting
  constexpr auto version = "#!TDV112";

  constexpr auto byte_order = 1; // default

  constexpr auto file_type             = 1; // 0=FULL, 1=GRID, 2=SOLUTION
  const auto     binary_title          = THIS::to_tecplot_binary_format(THIS::_configuration.title);
  const auto     num_variables         = THIS::_dimension;
  const auto     binary_variable_names = THIS::make_binary_grid_variable_string();

  constexpr auto zone_marker                          = 299.0f;
  const auto     binary_zone_title                    = THIS::to_tecplot_binary_format(zone_title);
  constexpr auto parent_zone                          = -1;  // default
  constexpr auto strand_id                            = 0;   // default for grid
  constexpr auto solution_time                        = 0.0; // default for grid
  constexpr auto zone_color                           = 1;   // default
  const auto     zone_type                            = THIS::binary_zone_type();
  constexpr auto specify_var_location                 = 0; // 0 = All NODAL, 1 = Specify, default for grid
  constexpr auto one_to_one_face_neighbor             = 0; // default
  constexpr auto user_define_face_neighbor_connection = 0; // default
  constexpr auto idim                                 = 0; // default
  constexpr auto jdim                                 = 0; // default
  constexpr auto kdim                                 = 0; // default
  constexpr auto auxilarily_name_index_pair           = 0; // default
  constexpr auto EOH_makrer                           = 357.0f;

  // I. HEADER SECTION

  // i. version number
  out_file << version;

  // ii. Integer value of 1.
  out_file << byte_order;

  // iii. Title and variable names
  out_file << file_type;
  out_file << binary_title;
  out_file << num_variables;
  out_file << binary_variable_names;

  // iiii. zones
  out_file << zone_marker;
  out_file << binary_zone_title;
  out_file << parent_zone;
  out_file << strand_id;
  out_file << solution_time;
  out_file << zone_color;
  out_file << zone_type;
  out_file << specify_var_location;
  out_file << one_to_one_face_neighbor;
  out_file << user_define_face_neighbor_connection;
  out_file << THIS::_num_nodes;
  out_file << THIS::_num_elements;
  out_file << idim;
  out_file << jdim;
  out_file << kdim;
  out_file << auxilarily_name_index_pair;
  out_file << EOH_makrer;

  // if (this->specify_variable_location_ == 0) // specify var location
  //   writer << 0;                             // Var location : Node
  // else
  //   writer << 1 << 1; // Var location : Cell center
}

void Writer::write_binary_grid_file_data(Binary_File& binary_writer, const double* coordinates_ptr)
{
  REQUIRE(THIS::_is_grid_data_saved, "Grid file data should be saved first");

  // Setting
  constexpr auto zone_marker                         = 299.0f;
  constexpr auto variable_data_format                = 2;  // double : 2
  constexpr auto has_passive_variable                = 0;  // default
  constexpr auto has_variable_sharing                = 0;  // default
  constexpr auto zone_number_to_shared_connecitivity = -1; // default

  const auto num_variables = THIS::_dimension;

  // II. DATA SECTION

  // zone
  binary_writer << zone_marker;

  for (int i = 0; i < num_variables; ++i)
  {
    binary_writer << variable_data_format;
  }

  binary_writer << has_passive_variable;
  binary_writer << has_variable_sharing;
  binary_writer << zone_number_to_shared_connecitivity;

  for (int i = 0; i < num_variables; ++i)
  {
    const auto* begin = coordinates_ptr + i * THIS::_num_nodes;
    const auto* end   = coordinates_ptr + (i + 1) * THIS::_num_nodes;

    const auto min_value = *std::min_element(begin, end);
    const auto max_value = *std::max_element(begin, end);
    binary_writer << min_value << max_value; // min,max value of each variable
  }

  const auto num_total_data = num_variables * THIS::_num_nodes;
  for (int i = 0; i < num_total_data; ++i)
  {
    binary_writer << coordinates_ptr[i];
  }
}

void Writer::write_binary_connectivity(Binary_File& out_file, const std::vector<std::vector<int>>& connectivities)
{
  for (const auto& connectivity : connectivities)
  {
    for (const auto node_index : connectivity)
    {
      // In binary files, node indices start at 0
      out_file << node_index - 1;
    }
  }
}

void Writer::write_binary_solution_file_header(Binary_File& out_file, const Solution_Data& info)
{
  REQUIRE(THIS::_is_grid_data_saved, "Grid file data should be saved first");

  const auto straind_id = std::to_string(THIS::_strand_id);
  const auto zone_title = THIS::_configuration.title + "_SOL" + straind_id;

  // Setting
  constexpr auto version = "#!TDV112";

  constexpr auto byte_order = 1; // default

  constexpr auto file_type             = 2; // 0=FULL, 1=GRID, 2=SOLUTION
  const auto     binary_title          = THIS::to_tecplot_binary_format(THIS::_configuration.title);
  const auto     num_variables         = info.num_variables();
  const auto     binary_variable_names = THIS::make_binary_variable_string(info.variable_strs);

  constexpr auto zone_marker                          = 299.0f;
  const auto     binary_zone_title                    = THIS::to_tecplot_binary_format(zone_title);
  constexpr auto parent_zone                          = -1; // default
  const auto     strand_id                            = THIS::_strand_id;
  const auto     solution_time                        = info.solution_time;
  constexpr auto zone_color                           = 1; // default
  const auto     zone_type                            = THIS::binary_zone_type();
  constexpr auto specify_var_location                 = 1; // 1 = Specify, default for solution
  constexpr auto one_to_one_face_neighbor             = 0; // default
  constexpr auto user_define_face_neighbor_connection = 0; // default
  constexpr auto idim                                 = 0; // default
  constexpr auto jdim                                 = 0; // default
  constexpr auto kdim                                 = 0; // default
  constexpr auto auxilarily_name_index_pair           = 0; // default
  constexpr auto EOH_makrer                           = 357.0f;

  // I. HEADER SECTION

  // i. version number
  out_file << version;

  // ii. Integer value of 1.
  out_file << byte_order;

  // iii. Title and variable names
  out_file << file_type;
  out_file << binary_title;
  out_file << num_variables;
  out_file << binary_variable_names;

  // iiii. zones
  out_file << zone_marker;
  out_file << binary_zone_title;
  out_file << parent_zone;
  out_file << strand_id;
  out_file << solution_time;
  out_file << zone_color;
  out_file << zone_type;
  out_file << specify_var_location;

  for (int i = 0; i < num_variables; ++i)
  {
    const auto var_location = info.var_locations_ptr[i];

    if (var_location == Variable_Location::NODE)
    {
      out_file << 0;
    }
    else if (var_location == Variable_Location::CELL_CENTER)
    {
      out_file << 1;
    }
    else
    {
      EXCEPTION("not supproted var var_index");
    }
  }

  out_file << one_to_one_face_neighbor;
  out_file << user_define_face_neighbor_connection;
  out_file << THIS::_num_nodes;
  out_file << THIS::_num_elements;
  out_file << idim;
  out_file << jdim;
  out_file << kdim;
  out_file << auxilarily_name_index_pair;
  out_file << EOH_makrer;
}

void Writer::write_binary_solution_file_data(Binary_File& out_file, const Solution_Data& info)
{
  REQUIRE(THIS::_is_grid_data_saved, "Grid file data should be saved first");

  // Setting
  constexpr auto zone_marker                         = 299.0f;
  constexpr auto variable_data_format                = 2;  // double : 2
  constexpr auto has_passive_variable                = 0;  // default
  constexpr auto has_variable_sharing                = 0;  // default
  constexpr auto zone_number_to_shared_connecitivity = -1; // default

  const auto num_variables = info.num_variables();

  // II. DATA SECTION

  out_file << zone_marker;

  for (int i = 0; i < num_variables; ++i)
  {
    out_file << variable_data_format;
  }

  out_file << has_passive_variable;
  out_file << has_variable_sharing;
  out_file << zone_number_to_shared_connecitivity;

  const auto* iter1 = info.values_ptr;
  const auto* iter2 = info.values_ptr;

  for (int i = 0; i < num_variables; ++i)
  {
    const auto var_location = info.var_locations_ptr[i];
    const auto num_datas    = THIS::num_solution_datas(var_location);

    iter2 += num_datas;

    const auto min_value = *std::min_element(iter1, iter2);
    const auto max_value = *std::max_element(iter1, iter2);
    out_file << min_value << max_value; // min,max value of each variable

    iter1 = iter2;
  }

  for (auto iter3 = info.values_ptr; iter3 != iter2; ++iter3)
  {
    out_file << *iter3;
  }
}

int Writer::binary_zone_type(void)
{
  // ZoneType
  // 0 = ORDERED,
  // 1 = FELINESEG,
  // 2 = FETRIANGLE,
  // 3 = FEQUADRILATERAL,
  // 4 = FETETRAHEDRON,
  // 5 = FEBRICK,

  switch (THIS::_dimension)
  {
  case 1: return 1; // "FELINESEG"
  case 2: return 2; // "FETRIANGLE"
  default: EXCEPTION("unsupported diemsnion"); return -1;
  }
}

std::vector<int> Writer::to_tecplot_binary_format(const std::string_view str)
{
  std::vector<int> tecplot_binary_format;
  tecplot_binary_format.insert(tecplot_binary_format.end(), str.begin(), str.end());
  tecplot_binary_format.push_back(0); // null terminated
  return tecplot_binary_format;
}

} // namespace ms::tecplot



