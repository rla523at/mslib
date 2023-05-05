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
  template <typename T>
  void set(std::string_view key, const T& value)
  {
    // key shold be upper case!
    const auto KEY = ms::string::upper_case(key);

    if constexpr (std::is_same_v<T, bool>)
    {
      if (this->key_to_bool_value_.contains(KEY.data()))
      {
        EXCEPTION(std::string(key) + "is duplicated");
      }

      this->key_to_bool_value_.insert({KEY, value});
    }
    else if constexpr (std::is_same_v<T, int>)
    {
      if (this->key_to_integer_value_.contains(KEY.data()))
      {
        EXCEPTION(std::string(key) + "is duplicated");
      }

      this->key_to_integer_value_.insert({KEY, value});
    }
    else if constexpr (std::is_same_v<T, double>)
    {
      if (this->key_to_double_value_.contains(KEY.data()))
      {
        EXCEPTION(std::string(key) + "is duplicated");
      }

      this->key_to_double_value_.insert({KEY, value});
    }
    else if constexpr (std::is_same_v<T, std::string>)
    {
      if (this->key_to_str_value_.contains(KEY.data()))
      {
        EXCEPTION(std::string(key) + "is duplicated");
      }

      this->key_to_str_value_.insert({KEY, value});
    }
    else if constexpr (std::is_same_v<T, Data>)
    {
      if (this->key_to_data_set_.contains(KEY.data()))
      {
        EXCEPTION(std::string(key) + "is duplicated");
      }

      this->key_to_data_set_.insert({KEY, value});
    }
    else
    {
      EXCEPTION("Not supported type");
    }
  }

public:
  template <typename T>
  T extract_data(std::string_view key) const
  {
    const auto KEY = ms::string::upper_case(key);

    if constexpr (std::is_same_v<T, bool>)
    {
      REQUIRE(this->key_to_bool_value_.contains(KEY), std::string(KEY) + "doesn't exist");
      return this->key_to_bool_value_.at(KEY);
    }
    else if constexpr (std::is_same_v<T, int>)
    {
      REQUIRE(this->key_to_integer_value_.contains(KEY), std::string(KEY) + "doesn't exist");
      return this->key_to_integer_value_.at(KEY);
    }
    else if constexpr (std::is_same_v<T, double>)
    {
      REQUIRE(this->key_to_double_value_.contains(KEY), std::string(KEY) + "doesn't exist");
      return this->key_to_double_value_.at(KEY);
    }
    else if constexpr (std::is_same_v<T, std::string>)
    {
      REQUIRE(this->key_to_str_value_.contains(KEY), std::string(KEY) + "doesn't exist");
      return this->key_to_str_value_.at(KEY);
    }
    else if constexpr (std::is_same_v<T, Data>)
    {
      REQUIRE(this->key_to_data_set_.contains(KEY), std::string(KEY) + "doesn't exist");
      return this->key_to_data_set_.at(KEY);
    }
    else
    {
      EXCEPTION("Not supported type");
      return NULL;
    }
  }

private:
  std::map<std::string, bool>        key_to_bool_value_;
  std::map<std::string, int>         key_to_integer_value_;
  std::map<std::string, double>      key_to_double_value_;
  std::map<std::string, std::string> key_to_str_value_;
  std::map<std::string, Data>    key_to_data_set_;
};
} // namespace ms::config