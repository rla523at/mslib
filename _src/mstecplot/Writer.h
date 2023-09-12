#pragma once
#include <string>
#include <vector>
#include "Data.h"
#include "Binary_File.h"


// class declaration
namespace ms::tecplot
{

enum class File_Format
{
  ASCII,
  Binary,
};

struct Configuration
{
  File_Format file_format;
  std::string folder_path;
  std::string title;
};

class Writer
{
public:
  static void set_configuration(Configuration&& configuration);
  static void write_grid_file(const Grid_Data& info);
  static void write_solution_file(const Solution_Data& info);
  // void        write_full_file(void) const;

private:
  // ASCII
  static void write_ASCII_grid_file_header(std::ofstream& out_file);
  static void write_ASCII_grid_file_data(std::ofstream& out_file, const double* coordinates_ptr);
  static void write_ASCII_connectivity(std::ofstream& out_file, const std::vector<std::vector<int>>& connectivities);
  static void write_ASCII_solution_file_header(std::ofstream& out_file, const Solution_Data& info);
  static void write_ASCII_solution_file_data(std::ofstream& out_file, const Solution_Data& info);

  static std::string make_ASCII_variable_string(const std::vector<std::string>& variable_strs);
  static std::string make_ASCII_grid_variable_string(void);

  // Binary
  static void write_binary_grid_file_header(Binary_File& out_file);
  static void write_binary_grid_file_data(Binary_File& out_file, const double* coordinates_ptr);
  static void write_binary_connectivity(Binary_File& out_file, const std::vector<std::vector<int>>& connectivities);
  static void write_binary_solution_file_header(Binary_File& out_file, const Solution_Data& info);
  static void write_binary_solution_file_data(Binary_File& out_file, const Solution_Data& info);

  static int              binary_zone_type(void);
  static std::vector<int> to_tecplot_binary_format(const std::string_view str);
  static std::vector<int> make_binary_grid_variable_string(void);
  static std::vector<int> make_binary_variable_string(const std::vector<std::string>& variable_strs);

  static int         num_solution_datas(const Variable_Location variable_location);
  static void        save_part_of_grid_data(const Grid_Data& info);
  static std::string zone_type_str(void);

private:
  inline static Configuration _configuration;
  inline static bool          _is_configuration_setup_complete = false;
  inline static int           _dimension                       = 0;
  inline static int           _num_nodes                       = 0;
  inline static int           _num_elements                    = 0;
  inline static bool          _is_grid_data_saved              = false;
  inline static int           _strand_id                       = 0;

private:
  using THIS = Writer;
};

} // namespace ms::tecplot
