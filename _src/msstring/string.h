#pragma once
#include "msexception/Exception.h"

#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace ms::string
{

inline constexpr int fail_to_find = -19940523;

bool compare_icase(const int c1, const int c2);
bool compare_icase(std::string_view str1, std::string_view str2);
bool contain(std::string_view str, const char target); // 형변환이 필요하면 template이 우선순위가 있기때문에 필요
bool contain(std::string_view str, const int target);
bool contain(std::string_view str, const char* target);        // 형변환이 필요하면 template이 우선순위가 있기때문에 필요
bool contain(std::string_view str, const std::string& target); // 형변환이 필요하면 template이 우선순위가 있기때문에 필요
bool contain(std::string_view str, std::string_view target);
bool contain_icase(std::string_view str, const char target); // 형변환이 필요하면 template이 우선순위가 있기때문에 필요
bool contain_icase(std::string_view str, const int target);
bool contain_icase(std::string_view str, const char* target); // 형변환이 필요하면 template이 우선순위가 있기때문에 필요
bool contain_icase(std::string_view str, std::string_view target);
bool contain_icase(std::string_view str, const std::string& target); // 형변환이 필요하면 template이 우선순위가 있기때문에 필요

// When fail to find target in str, return ms::sel::failt_to_find
int find_nth_position(std::string_view str, const char target, const int n);
int find_nth_position(std::string_view str, const int target, const int n) = delete;
int find_nth_position(std::string_view str, std::string_view target, const int n);
int find_nth_position_icase(std::string_view str, const char target, const int n);
int find_nth_position_icase(std::string_view str, const int target, const int n) = delete;
int find_nth_position_icase(std::string_view str, std::string_view target, const int n);
int find_r_nth_position(std::string_view str, std::string_view target, const int n);

bool is_natural_number(std::string_view str);
bool is_integer(std::string_view str);
bool is_real_number(std::string_view str);

// It doesn't make empty str
std::vector<std::string_view> parse_by(std::string_view str, const char delimiter);

std::string      remove(std::string_view str, const char target);
std::string      remove(std::string_view str, std::string_view target);
std::string_view remove_after(std::string_view str, std::string_view target);
void             remove_inplace(std::string& str, const char target);
void             remove_inplace(std::string& str, std::string_view target);

std::string replace(std::string_view str, const char target, const char replacement);
std::string replace(std::string_view str, std::string_view target, std::string_view replacement);
void        replace_inplace(std::string& str, const char target, const char replacement);
void        replace_inplace(std::string& str, std::string_view target, std::string_view replacement);

int         upper_case(const int c);
std::string upper_case(std::string_view str);
void        upper_case_inplace(int& c);
void        upper_case_inplace(std::string& str);

template <typename... Args>
bool contain(std::string_view str, const Args... args)
{
  return (ms::string::contain(str, args) && ...);
};

template <typename... Args>
bool contain_icase(std::string_view str, const Args... args)
{
  return (ms::string::contain_icase(str, args) && ...);
};

template <typename... Options>
std::string double_to_str(const double value, const Options... options)
{
  std::stringstream sstream;
  (sstream << ... << options);
  sstream << value;
  return sstream.str();
}

template <typename T>
T str_to_value(std::string_view str)
{
  if constexpr (std::is_same_v<T, std::string>)
  {
    return str.data();
  }
  else if constexpr (std::is_same_v<T, bool>)
  {
    if (ms::string::compare_icase(str, "true"))
    {
      return true;
    }
    else if (ms::string::compare_icase(str, "false"))
    {
      return false;
    }
    else
    {
      EXCEPTION(std::string(str) + " doesn't have boolean format");
      return NULL;
    }
  }
  else
  {
    std::istringstream iss(str.data());
    T                  value;
    iss >> value;
    return value;
  }

};

} // namespace ms::string