#pragma once

#include <algorithm>
#include <cassert>
#include <chrono>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#if !defined( BLOCK )
#define BLOCK( Description )
#endif

namespace mslib
{
  using uint64 = unsigned long long;
  using utf8   = char;
} // namespace mslib
