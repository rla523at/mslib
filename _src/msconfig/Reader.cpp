#include "Reader.h"

#include "Data.h"

#include <iostream>

namespace
{

enum class Value_Type
{
  BOOL,
  INT,
  DOUBLE,
  STRING,
  DATA_SET,
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
  else if (ms::string::is_integer(value_str))
  {
    return Value_Type::INT;
  }
  else if (ms::string::is_real_number(value_str))
  {
    return Value_Type::DOUBLE;
  }
  else
  {
    EXCEPTION(std::string(value_str) + " has unsupproted format");
    return Value_Type::VOID;
  }
}

} // namespace

namespace ms::config
{

Data Reader::read(std::string_view file_path) const
{
  Data data;

  std::ifstream fstream(file_path.data());
  REQUIRE(fstream.is_open(), "configuration fstream should be succesfully open");

  std::string str;
  while (std::getline(fstream, str))
  {
    // remove comment
    constexpr auto comment_tag = "//";
    str                        = ms::string::remove_after(str, comment_tag);

    // remove space and tab
    ms::string::remove_inplace(str, ' ');
    ms::string::remove_inplace(str, "\t");

    // ignore empty line
    if (str.empty())
    {
      continue;
    }

    // set data
    if (this->is_key_value_foramt(str))
    {
      this->set_key_value(data, str);
    }
    else
    {
      try
      {
        auto sub_data_set = this->make_sub_data_set(fstream);
        data.set<Data>(str, sub_data_set);
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

bool Reader::is_key_value_foramt(std::string_view str) const
{
  return ms::string::contain(str, ":");
}

bool Reader::is_data_set_start_format(std::string_view str) const
{
  return str == "{";
}

bool Reader::is_data_set_end_format(std::string_view str) const
{
  return str == "}";
}

Data Reader::make_sub_data_set(std::ifstream& fstream) const
{
  Data sub_data;

  bool not_found_start_format = true;

  std::string str;
  while (std::getline(fstream, str))
  {
    // remove comment
    constexpr auto comment_tag = "//";
    str                        = ms::string::remove_after(str, comment_tag);

    // remove space and tab
    ms::string::remove_inplace(str, ' ');
    ms::string::remove_inplace(str, "\t");

    // ignore empty str
    if (str.empty())
    {
      continue;
    }

    // check start format
    if (not_found_start_format)
    {
      if (this->is_data_set_start_format(str))
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
    if (this->is_data_set_end_format(str))
    {
      return sub_data;
    }

    // set data
    if (this->is_key_value_foramt(str))
    {
      this->set_key_value(sub_data, str);
    }
    else
    {
      const auto& key = str;

      try
      {
        auto ssub_data_set = this->make_sub_data_set(fstream);
        sub_data.set<Data>(str, ssub_data_set);
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

void Reader::set_key_value(Data& data_set, std::string_view key_value_str) const
{
  const auto parsed_str = ms::string::parse_by(key_value_str, ':');
  REQUIRE(parsed_str.size() == 2, "str should have 2 part, key-value");

  const auto& key_str   = parsed_str[0];
  const auto& value_str = parsed_str[1];

  const auto value_type = evaluate_value_type(value_str);

  switch (value_type)
  {
  case Value_Type::BOOL:
    data_set.set<bool>(key_str, ms::string::str_to_value<bool>(value_str));
    break;
  case Value_Type::INT:
    data_set.set<int>(key_str, ms::string::str_to_value<int>(value_str));
    break;
  case Value_Type::DOUBLE:
    data_set.set<double>(key_str, ms::string::str_to_value<double>(value_str));
    break;
  case Value_Type::STRING:
    data_set.set<std::string>(key_str, std::string(value_str));
    break;
  default:
    EXCEPTION("fail to set data from " + std::string(key_value_str));
    break;
  }
}

} // namespace ms::config