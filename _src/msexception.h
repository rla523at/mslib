#pragma once
#include <cassert>
#include <format>
#include <iostream>
#include <string_view>

#ifdef _DEBUG

#define REQUIRE( requirement, format_string, ... )                                                      \
  if ( !( requirement ) )                                                                               \
  {                                                                                                     \
    setlocale( LC_ALL, "" );                                                                            \
                                                                                                        \
    std::string_view function_call_file_name = __FILE__;                                                \
    const size_t     num_remove              = function_call_file_name.rfind( "\\" ) + 1;               \
    function_call_file_name.remove_prefix( num_remove );                                                \
                                                                                                        \
    std::cout << "\n==============================EXCEPTION========================================\n"; \
    std::cout << std::format( "{:10}: {}", "File", function_call_file_name ) << "\n";                   \
    std::cout << std::format( "{:10}: {}", "Function", __FUNCTION__ ) << "\n";                          \
    std::cout << std::format( "{:10}: {}", "Line", __LINE__ ) << "\n";                                  \
    ms::exception::print_message( format_string __VA_OPT__(, __VA_ARGS__ ) );                           \
    std::cout << "==============================EXCEPTION========================================\n\n"; \
                                                                                                        \
    assert( false );                                                                                    \
  }                                        


#define EXCEPTION( format_string, ... ) REQUIRE( false, format_string, __VA_ARGS__ )

#else
#define REQUIRE( requirement, format_string, ... )
#define EXCEPTION( format_string, ... )
#endif

namespace ms::exception
{
  template <typename... Ts>
  void print_message( const std::string_view format_string, const Ts&... args )
  {
    std::cout << std::format( "{:10}: {}", "Message", std::vformat( format_string, std::make_format_args( args... ) ) ) << "\n";
  }

  template <typename... Ts>
  void print_message( const std::wstring_view format_string, const Ts&... args )
  {
    std::wcout << std::format( L"{:10}: {}", L"Message", std::vformat( format_string, std::make_wformat_args( args... ) ) ) << "\n";
  }

} // namespace ms::exception



