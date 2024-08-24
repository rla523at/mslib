#pragma once
#include "Config_Data.h"
#include <fstream>
#include <string_view>

namespace ms::config
{

// static class
class Reader
{
public:
  static Data read(std::string_view file_path);

private:
  static bool is_data_set_start_format(std::string_view str);
  static bool is_data_set_end_format(std::string_view str);
  static bool is_key_value_foramt(std::string_view str);
  static Data make_sub_data(std::ifstream& fstream);
  static void set_key_value(Data& data_set, std::string_view key_value_str);

private:
  Reader(void) = delete;

private:
  using THIS = Reader;
};

} // namespace ms::config