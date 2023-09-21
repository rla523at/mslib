#pragma once
#include "msstring/string.h"

#include <map>
#include <string>
#include <string_view>
#include <type_traits>

namespace ms::config
{

class Data
{
public:
  void set_string(std::string_view key, std::string_view str);
  void set_number(std::string_view key, std::string_view number);
  void set_bool(std::string_view key, const bool boolean);
  void set_sub_data(std::string_view key, Data&& sub_data);

public:
  template <typename T>
  const T& get_data(std::string_view key) const
  {
    const auto KEY = ms::string::upper_case(key);

    if constexpr (std::is_same_v<T, bool>)
    {
      REQUIRE(this->_key_to_bool.contains(KEY), std::string(KEY) + "doesn't exist");
      return this->_key_to_bool.at(KEY);
    }
    else if constexpr (std::is_same_v<T, int> || std::is_same_v<T, double>)
    {
      REQUIRE(this->_key_to_number.contains(KEY), std::string(KEY) + "doesn't exist");
      return ms::string::str_to_value<T>(this->_key_to_number.at(KEY));
    }
    else if constexpr (std::is_same_v<T, std::string>)
    {
      REQUIRE(this->_key_to_str.contains(KEY), std::string(KEY) + "doesn't exist");
      return this->_key_to_str.at(KEY);
    }
    else if constexpr (std::is_same_v<T, Data>)
    {
      REQUIRE(this->_key_to_data.contains(KEY), std::string(KEY) + "doesn't exist");
      return this->_key_to_data.at(KEY);
    }
    else
    {
      EXCEPTION("Not supported type");
      return NULL;
    }
  }

private:
  std::map<std::string, bool>        _key_to_bool;
  std::map<std::string, std::string> _key_to_number;
  std::map<std::string, std::string> _key_to_str;
  std::map<std::string, Data>        _key_to_data;
};
} // namespace ms::config