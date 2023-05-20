#include "Reader.h"

#include <iostream>

namespace
{

enum class Value_Type
{
  BOOL,
  NUMBER,
  STRING,
  DATA,
  VOID
};

Value_Type evaluate_value_type(std::string_view value_str)
{
  if (value_str.front() == '"' && value_str.back() == '"')
  {
    return Value_Type::STRING;
  }
  else if (ms::string::compare_icase(value_str, "true") || ms::string::compare_icase(value_str, "false"))
  {
    return Value_Type::BOOL;
  }
  else if (ms::string::is_real_number(value_str))
  {
    return Value_Type::NUMBER;
  }
  else
  {
    EXCEPTION(std::string(value_str) + " has unsupproted value type");
    return Value_Type::VOID;
  }
}

} // namespace

namespace ms::config
{

Data Reader::read(std::string_view file_path)
{
  constexpr auto comment_tag = "//";
  constexpr auto space       = ' ';
  constexpr auto tab         = '\t';

  Data data;

  std::ifstream file(file_path.data());
  REQUIRE(file.is_open(), "configuration file should be succesfully open");

  std::string str;
  while (std::getline(file, str))
  {
    ms::string::remove_after_inplace(str, comment_tag);
    ms::string::remove_inplace(str, space);
    ms::string::remove_inplace(str, tab);

    if (str.empty()) continue;

    // set data
    if (THIS::is_key_value_foramt(str))
    {
      THIS::set_key_value(data, str);
    }
    else
    {
      try
      {
        auto sub_data_set = THIS::make_sub_data(file);
        data.set_sub_data(str, std::move(sub_data_set));
      }
      catch (const std::exception& except)
      {
        std::cout << "fail to set " << str << "\n";
        std::cout << except.what();
        std::exit(523);
      }
    }
  }

  return data;
}

bool Reader::is_key_value_foramt(std::string_view str)
{
  return ms::string::contain(str, ":");
}

bool Reader::is_data_set_start_format(std::string_view str)
{
  return str == "{";
}

bool Reader::is_data_set_end_format(std::string_view str)
{
  return str == "}";
}

Data Reader::make_sub_data(std::ifstream& fstream)
{
  constexpr auto comment_tag = "//";

  Data sub_data;

  bool not_found_start_format = true;

  std::string str;
  while (std::getline(fstream, str))
  {
    ms::string::remove_after_inplace(str, comment_tag);
    ms::string::remove_inplace(str, ' ');
    ms::string::remove_inplace(str, "\t");

    if (str.empty()) continue;

    // check start format
    if (not_found_start_format)
    {
      if (THIS::is_data_set_start_format(str))
      {
        not_found_start_format = false;
        continue;
      }
      else
      {
        EXCEPTION("data set should be start with { \n" + str + " is appeared before { ");
      }
    }

    // check end format
    if (THIS::is_data_set_end_format(str))
    {
      return sub_data;
    }

    // set data
    if (THIS::is_key_value_foramt(str))
    {
      THIS::set_key_value(sub_data, str);
    }
    else
    {
      const auto& key = str;

      try
      {
        auto ssub_data = THIS::make_sub_data(fstream);
        sub_data.set_sub_data(str, std::move(ssub_data));
      }
      catch (const std::exception& except)
      {
        std::cout << "fail to set " << str << "\n";
        std::cout << except.what();
        std::exit(523);
      }
    }
  }

  EXCEPTION("can't not find }");
  return sub_data;
}

void Reader::set_key_value(Data& data, std::string_view key_value_str)
{
  const auto parsed_str = ms::string::parse_by(key_value_str, ':');
  REQUIRE(parsed_str.size() == 2, "key-value str should have 2 part");

  const auto& key_str   = parsed_str[0];
  const auto& value_str = parsed_str[1];

  const auto value_type = evaluate_value_type(value_str);

  switch (value_type)
  {
  case Value_Type::BOOL:
    data.set_bool(key_str, ms::string::str_to_value<bool>(value_str));
    break;
  case Value_Type::NUMBER:
    data.set_number(key_str, std::string(value_str));
    break;
  case Value_Type::STRING:
    data.set_string(key_str, std::string(value_str));
    break;
  default:
    EXCEPTION("fail to set data from " + std::string(key_value_str));
    break;
  }
}

} // namespace ms::config