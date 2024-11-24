#include "msstring.h"

#include <algorithm>
#include <functional>
#include <iomanip>
#include <sstream>

namespace ms::string
{

class Case_Insensitive_Comparator
{
public:
  bool operator()(const char c1, const char c2) const
  {
    return ms::string::compare_icase(c1, c2);
  }
};

bool compare_icase(const int c1, const int c2)
{
  return ms::string::upper_case(c1) == ms::string::upper_case(c2);
}

bool compare_icase(std::string_view str1, std::string_view str2)
{
  const auto str_size = str1.size();

  if (str_size != str2.size())
  {
    return false;
  }

  for (size_t i = 0; i < str_size; i++)
  {
    if (ms::string::upper_case(str1[i]) != ms::string::upper_case(str2[i]))
    {
      return false;
    }
  }

  return true;
}

bool contain(std::string_view str, const char target)
{
  return str.find(target) != std::string_view::npos;
}

bool contain(std::string_view str, const int target)
{
  return str.find(static_cast<char>(target)) != std::string_view::npos;
};

bool contain(std::string_view str, const char* target)
{
  const auto target_begin = target;
  const auto target_end   = target + std::strlen(target);

  std::boyer_moore_searcher searcher(target_begin, target_end);
  return std::search(str.begin(), str.end(), searcher) != str.end();
}

bool contain(std::string_view str, const std::string& target)
{
  std::boyer_moore_searcher searcher(target.begin(), target.end());
  return std::search(str.begin(), str.end(), searcher) != str.end();
}

bool contain(std::string_view str, std::string_view target)
{
  std::boyer_moore_searcher searcher(target.begin(), target.end());
  return std::search(str.begin(), str.end(), searcher) != str.end();
};

bool contain_icase(std::string_view str, const char target)
{
  return contain_icase(str, static_cast<int>(target));
}

bool contain_icase(std::string_view str, const int target)
{
  for (size_t i = 0; i < str.size(); ++i)
  {
    if (ms::string::compare_icase(str[i], target))
    {
      return true;
    }
  }

  return false;
}

bool contain_icase(std::string_view str, const char* target)
{
  const auto            target_begin = target;
  const auto            target_end   = target + std::strlen(target);
  std::default_searcher searcher(target_begin, target_end, Case_Insensitive_Comparator());

  if (std::search(str.begin(), str.end(), searcher) == str.end())
  {
    return false;
  }
  else
  {
    return true;
  }
}

bool contain_icase(std::string_view str, std::string_view target)
{
  return ms::string::contain_icase(str, target.data());
}

bool contain_icase(std::string_view str, const std::string& target)
{
  return ms::string::contain_icase(str, target.data());
}

int find_nth_position(std::string_view str, const char target, const int n)
{
  if (str.empty())
  {
    return ms::string::fail_to_find;
  }

  if (n < 1)
  {
    return ms::string::fail_to_find;
  }

  auto count = 0;
  for (int i = 0; i < str.size(); i++)
  {
    if (str[i] == target)
    {
      if (n == ++count)
      {
        return i;
      }
    }
  }

  return ms::string::fail_to_find;
}

int find_nth_position(std::string_view str, std::string_view target, const int n)
{
  if (str.empty() || target.empty())
  {
    return ms::string::fail_to_find;
  }

  if (n < 1)
  {
    return ms::string::fail_to_find;
  }

  std::boyer_moore_searcher searcher(target.begin(), target.end());

  const auto str_begin   = str.begin();
  const auto target_size = target.size();

  auto       search_start_iter = str_begin;
  const auto search_end_iter   = str.end();

  // iter := str iter pointing ith target begin position
  auto iter = str_begin;

  for (size_t i = 0; i < n; ++i)
  {
    iter = std::search(search_start_iter, search_end_iter, searcher);

    if (iter == search_end_iter)
    {
      return ms::string::fail_to_find;
    }

    // str iter pointing ith target end position
    search_start_iter = iter + target_size;
  }

  const auto pos = static_cast<int>(iter - str_begin);
  return pos;
}

int find_nth_position_icase(std::string_view str, const char target, const int n)
{
  if (str.empty())
  {
    return ms::string::fail_to_find;
  }

  if (n < 1)
  {
    return ms::string::fail_to_find;
  }

  auto count = 0;
  for (int i = 0; i < str.size(); ++i)
  {
    if (ms::string::compare_icase(str[i], target))
    {
      if (n == ++count)
      {
        return i;
      }
    }
  }

  return ms::string::fail_to_find;
}

int find_nth_position_icase(std::string_view str, std::string_view target, const int n)
{
  if (str.empty() || target.empty())
  {
    return ms::string::fail_to_find;
  }

  if (n < 1)
  {
    return ms::string::fail_to_find;
  }

  std::default_searcher searcher(target.begin(), target.end(), Case_Insensitive_Comparator());

  const auto str_begin   = str.begin();
  const auto target_size = target.size();

  auto       search_start_iter = str_begin;
  const auto search_end_iter   = str.end();

  // iter := str iter pointing ith target begin position
  auto iter = str_begin;

  for (size_t i = 0; i < n; ++i)
  {
    iter = std::search(search_start_iter, search_end_iter, searcher);

    if (iter == search_end_iter)
    {
      return ms::string::fail_to_find;
    }

    // search_start_iter := str iter pointing ith target end position
    search_start_iter = iter + target_size;
  }

  const auto pos = static_cast<int>(iter - str_begin);
  return pos;
}

int find_r_nth_position(std::string_view str, std::string_view target, const int n)
{
  if (str.empty() || target.empty())
  {
    return ms::string::fail_to_find;
  }

  if (n < 1)
  {
    return ms::string::fail_to_find;
  }

  const auto target_size = target.size();
  const auto str_rbegin  = str.rbegin();
  const auto str_rend    = str.rend();

  std::boyer_moore_searcher searcher(target.rbegin(), target.rend());

  auto rpos = 0;
  for (size_t i = 0; i < n; ++i)
  {
    // str riter pointing ith target rbegin position
    const auto iter1 = std::search(str_rbegin + rpos, str_rend, searcher);

    if (iter1 == str_rend)
    {
      return ms::string::fail_to_find;
    }

    // str riter pointing ith target rend position
    const auto iter2 = iter1 + target_size;

    // ith target rpos
    rpos = static_cast<int>(iter2 - str_rbegin);
  }

  const auto pos = static_cast<int>(str.size()) - rpos;
  return pos;
}

std::vector<std::wstring_view> parse_by( std::wstring_view wstr, const wchar_t delimiter )
{
  std::vector<std::wstring_view> result;

  if ( wstr.empty() )
  {
    result.push_back( wstr );
    return result;
  }

  const auto str_begin = wstr.begin();
  const auto str_end   = wstr.end();

  auto front_iter = str_begin;
  auto back_iter  = str_begin;

  while ( true )
  {
    if ( *back_iter == delimiter )
    {
      const auto pos   = front_iter - str_begin;
      const auto count = back_iter - front_iter;

      if ( count != 0 )
      {
        result.push_back( wstr.substr( pos, count ) );
      }

      front_iter = ++back_iter;
    }
    else
    {
      ++back_iter;
    }

    if ( back_iter == wstr.end() )
    {
      const auto pos   = front_iter - str_begin;
      const auto count = back_iter - front_iter;

      if ( count != 0 )
      {
        result.push_back( wstr.substr( pos, count ) );
      }

      break;
    }
  }

  return result;
}

std::vector<std::string_view> parse_by( std::string_view str, const char delimiter )
{
  std::vector<std::string_view> result;

  if (str.empty())
  {
    return result;
  }

  const auto str_begin = str.begin();
  const auto str_end   = str.end();

  auto front_iter = str_begin;
  auto back_iter  = str_begin;

  while (true)
  {
    if (*back_iter == delimiter)
    {
      const auto pos   = front_iter - str_begin;
      const auto count = back_iter - front_iter;

      if (count != 0)
      {
        result.push_back(str.substr(pos, count));
      }

      front_iter = ++back_iter;
    }
    else
    {
      ++back_iter;
    }

    if (back_iter == str.end())
    {
      const auto pos   = front_iter - str_begin;
      const auto count = back_iter - front_iter;

      if (count != 0)
      {
        result.push_back(str.substr(pos, count));
      }

      break;
    }
  }

  return result;
};

std::string remove(std::string_view str, const char target)
{
  std::string result;
  result.reserve(str.size());

  for (const auto c : str)
  {
    if (c != target)
    {
      result.push_back(c);
    }
  }

  return result;
}

std::string remove(std::string_view str, std::string_view target)
{
  std::string result;
  result.reserve(str.size());

  while (true)
  {
    const auto target_start_pos = ms::string::find_nth_position(str, target, 1);

    if (target_start_pos == ms::string::fail_to_find)
    {
      result.append(str);
      break;
    }
    else
    {
      result.append(str.substr(0, target_start_pos));
      str.remove_prefix(target_start_pos + target.size());
    }
  }

  result.shrink_to_fit();
  return result;
}

std::string_view remove_after(std::string_view str, std::string_view target)
{
  constexpr auto n   = 1;
  const auto     pos = ms::string::find_nth_position(str, target, n);

  if (pos == ms::string::fail_to_find)
  {
    return str;
  }

  const auto num_remove = str.size() - pos;
  str.remove_suffix(num_remove);
  return str;
}

void remove_after_inplace(std::string& str, std::string_view target)
{
  constexpr auto n   = 1;
  const auto     pos = ms::string::find_nth_position(str, target, n);

  if (pos == ms::string::fail_to_find) return;

  str.erase(pos);
}

std::string_view remove_before(std::string_view str, std::string_view target)
{
  constexpr auto n   = 1;
  const auto     pos = ms::string::find_nth_position(str, target, n);

  if (pos == ms::string::fail_to_find)
  {
    return str;
  }

  const auto num_remove = pos;
  str.remove_prefix(num_remove);
  return str;
}

std::string_view remove_up_to(std::string_view str, std::string_view target)
{
  constexpr auto n   = 1;
  const auto     pos = ms::string::find_nth_position(str, target, n);

  if (pos == ms::string::fail_to_find)
  {
    return str;
  }

  const auto num_remove = pos + target.size();
  str.remove_prefix(num_remove);
  return str;
}

void remove_inplace(std::string& str, const char target)
{
  str = ms::string::remove(str, target);
};

void remove_inplace(std::string& str, std::string_view target)
{
  str = ms::string::remove(str, target);
}

std::string replace(std::string_view str, const char target, const char replacement)
{
  std::string result;
  result.reserve(str.size());

  for (size_t i = 0; i < str.size(); i++)
  {
    if (str[i] == target)
    {
      result.push_back(replacement);
    }
    else
    {
      result.push_back(str[i]);
    }
  }

  return result;
}

std::string replace(std::string_view str, std::string_view target, std::string_view replacement)
{
  const auto expected_maximum_size = str.size() * 2;

  std::string result;
  result.reserve(expected_maximum_size);

  while (true)
  {
    const auto target_start_pos = ms::string::find_nth_position(str, target, 1);

    if (target_start_pos == ms::string::fail_to_find)
    {
      result.append(str);
      break;
    }
    else
    {
      result.append(str.substr(0, target_start_pos));
      result.append(replacement);
      str.remove_prefix(target_start_pos + target.size());
    }
  }

  result.shrink_to_fit();
  return result;
}

void replace_inplace(std::string& str, const char target, const char replacement)
{
  for (auto& c : str)
  {
    if (c == target)
    {
      c = replacement;
    }
  }
}

void replace_inplace(std::string& str, std::string_view target, std::string_view replacement)
{
  str = ms::string::replace(str, target, replacement);
}

int upper_case(const int c)
{
  return std::toupper(c);
}

std::string upper_case(std::string_view str)
{
  const auto size = str.size();

  std::string result;
  result.resize(size);

  for (size_t i = 0; i < size; i++)
  {
    result[i] = static_cast<char>(ms::string::upper_case(str[i]));
  }

  return result;
}

void upper_case_inplace(int& c)
{
  c = std::toupper(c);
}

void upper_case_inplace(std::string& str)
{
  for (auto& c : str)
  {
    c = static_cast<char>(ms::string::upper_case(c));
  }
}

bool is_natural_number(std::string_view str)
{
  if (str.empty())
  {
    return false;
  }

  if (str[0] == '+')
  {
    str.remove_prefix(1);
  }

  if (str.empty())
  {
    return false;
  }

  return std::all_of(str.begin(), str.end(), ::isdigit); // 왜 ::isdigit해야되지??
}

bool is_integer(std::string_view str)
{
  if (str.empty())
  {
    return false;
  }

  if (str[0] == '+' || str[0] == '-')
  {
    str.remove_prefix(1);
  }

  if (str.empty())
  {
    return false;
  }

  return std::all_of(str.begin(), str.end(), ::isdigit);
}

bool is_real_number(std::string_view str)
{
  // check empty
  if (str.empty()) return false;

  // check and remove sign
  if (str[0] == '+' || str[0] == '-')
  {
    str.remove_prefix(1);
  }

  // check empty
  if (str.empty()) return false;

  bool is_first_dot = true;

  for (const auto c : str)
  {
    if (std::isdigit(c)) continue;

    if (c == '.' && is_first_dot)
    {
      is_first_dot = false;
    }
    else
    {
      return false;
    }

  }

  return true;
}

} // namespace ms::string