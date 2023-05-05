#pragma once
#include <fstream>
#include <string_view>

// forward declaration
namespace ms::config
{
class Data;
} // namespace ms::config

// class declaration
namespace ms::config
{

class Reader
{
public:
  Data read(std::string_view file_path) const;

private:
  bool       is_data_set_start_format(std::string_view str) const;
  bool       is_data_set_end_format(std::string_view str) const;
  bool       is_key_value_foramt(std::string_view str) const;
  Data   make_sub_data_set(std::ifstream& fstream) const;
  void       set_key_value(Data& data_set, std::string_view key_value_str) const;
};

} // namespace ms::config