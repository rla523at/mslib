#include "msfilesystem.h"

#include "msexception.h"
#include "msstring.h"

#include <fstream>

namespace ms::filesystem
{
  std::wstring_view base_name_view( std::wstring_view file_path )
  {
    const size_t folder_path_end_pos = file_path.rfind( L'/' );

    if ( folder_path_end_pos != std::wstring_view::npos )
      file_path.remove_prefix( folder_path_end_pos + 1 );

    const size_t base_name_end_pos = file_path.find( L'.' );

    const std::wstring_view result = file_path.substr( 0, base_name_end_pos );
    return result;
  }

  void copy_file( const std::string_view from_file_path, const std::string_view to_file_path, const std::filesystem::copy_options option )
  {
    REQUIRE( !is_folder_path( from_file_path ), std::string( from_file_path ) + " should be a file" );
    REQUIRE( !is_folder_path( to_file_path ), std::string( to_file_path ) + " should be a file" );
    REQUIRE( is_exist_file( from_file_path ), std::string( from_file_path ) + " is not exist file" );
    REQUIRE( is_exist_folder( extract_folder_path( to_file_path ) ), extract_folder_path( to_file_path ) + " is not exist folder" );
    std::filesystem::copy_file( from_file_path, to_file_path, option );
  }

  void replace_file( const std::string_view from_file_path, const std::string_view to_file_path )
  {
    REQUIRE( !is_folder_path( from_file_path ), std::string( from_file_path ) + " should be a file" );
    REQUIRE( !is_folder_path( to_file_path ), std::string( to_file_path ) + " should be a file" );
    REQUIRE( is_exist_file( from_file_path ), std::string( from_file_path ) + " is not exist file" );
    REQUIRE( is_exist_file( to_file_path ), std::string( to_file_path ) + " is not exist file" );
    std::filesystem::copy_file( from_file_path, to_file_path, std::filesystem::copy_options::overwrite_existing );
  }

  std::wstring extract_folder_path( const std::wstring_view file_path )
  {
    REQUIRE( !is_folder_path( file_path ), L"input({}) should be file path", file_path );
    return std::filesystem::path( file_path ).parent_path().wstring() + L"/";
  }

  std::string extract_folder_path( const std::string_view file_path )
  {
    REQUIRE( !is_folder_path( file_path ), "input should be file path" );
    return std::filesystem::path( file_path ).parent_path().string() + "/";
  }

  std::wstring extract_file_name( const std::wstring_view file_path )
  {
    REQUIRE( !is_folder_path( file_path ), L"input({}) should be file path", file_path );

    return std::filesystem::path( file_path ).filename().wstring();
  }

  std::string extract_file_name( const std::string_view file_path )
  {
    REQUIRE( !is_folder_path( file_path ), "input should be file path" );

    return std::filesystem::path( file_path ).filename().string();
  }

  std::wstring_view file_name_view( std::wstring_view file_path )
  {
    REQUIRE( !is_folder_path( file_path ), L"input({}) should be file path", file_path );

    const size_t folder_path_end_pos = file_path.rfind( L'/' );

    if ( folder_path_end_pos != std::wstring_view::npos )
      file_path.remove_prefix( folder_path_end_pos + 1 );

    return file_path;
  }

  std::vector<std::string> file_names_in_folder( const std::string_view folder_path )
  {
    REQUIRE( is_folder_path( folder_path ), "folder_path should be end with /" );

    std::vector<std::string> file_names;

    std::filesystem::directory_iterator iter( folder_path );
    while ( iter != std::filesystem::end( iter ) )
    {
      const auto& entry = *iter;

      if ( !entry.is_directory() )
      {
        file_names.push_back( entry.path().filename().string() );
      }

      iter++;
    }

    return file_names;
  }

  std::wstring_view folder_path_view( std::wstring_view file_path )
  {
    REQUIRE( !is_folder_path( file_path ), L"input({}) should be file path", file_path );

    const size_t folder_path_end_pos = file_path.rfind( L'/' );

    const bool there_are_no_folder_path = folder_path_end_pos == std::wstring_view::npos;
    if ( there_are_no_folder_path )
      return {};

    const std::wstring_view result = file_path.substr( 0, folder_path_end_pos + 1 );
    return result;
  }

  std::vector<std::string> folder_names_in_folder( const std::string_view folder_path )
  {
    REQUIRE( is_folder_path( folder_path ), "folder_path should be end with /" );

    std::vector<std::string> folder_names;

    std::filesystem::directory_iterator iter( folder_path );
    while ( iter != std::filesystem::end( iter ) )
    {
      const auto& entry = *iter;

      if ( entry.is_directory() )
      {
        folder_names.push_back( entry.path().filename().string() );
      }

      iter++;
    }

    return folder_names;
  }

  std::vector<std::string> file_paths_in_folder( const std::string& folder_path )
  {
    REQUIRE( is_folder_path( folder_path ), "folder_path should be end with /" );

    std::vector<std::string> file_name_text;

    std::filesystem::directory_iterator iter( folder_path );
    while ( iter != std::filesystem::end( iter ) )
    {
      const auto& entry = *iter;

      if ( entry.is_directory() )
      {
        iter++;
        continue;
      }

      file_name_text.push_back( entry.path().string() );
      iter++;
    }

    return file_name_text;
  }

  bool has_this_extension( const std::string_view file_path, const std::string_view extension )
  {
    REQUIRE( extension.front() == '.', std::string( extension ) + "is not extension format. extension should be start with ." );
    return file_path.ends_with( extension );
  }

  bool is_exist_folder( const std::wstring_view folder_path )
  {
    REQUIRE( is_folder_path( folder_path ), "folder_path should be end with /" );

    std::filesystem::path p( folder_path );
    return std::filesystem::exists( p );
  }

  bool is_exist_folder( const std::string_view folder_path )
  {
    REQUIRE( is_folder_path( folder_path ), "folder_path should be end with /" );

    std::filesystem::path p( folder_path );
    return std::filesystem::exists( p );
  }

  bool is_exist_file( const std::string_view file_path )
  {
    REQUIRE( !is_folder_path( file_path ), "file path should not be end with /" );

    std::filesystem::path p( file_path );
    return std::filesystem::exists( p );
  }

  bool is_folder_path( const std::string_view folder_path )
  {
    const auto last_char = folder_path.back();
    return last_char == '/';
  }

  bool is_folder_path( const std::wstring_view folder_path )
  {
    const auto last_char = folder_path.back();
    return last_char == '/';
  }

  void make_folder( const std::string_view folder_path )
  {
    REQUIRE( is_folder_path( folder_path ), "folder_path should be end with /" );

    if ( folder_path.empty() ) return;

    if ( ms::filesystem::is_exist_folder( folder_path ) ) return;

    std::filesystem::path p( folder_path );
    std::filesystem::create_directories( p );
  }

  void move_file( const std::string_view file_path, const std::string_view new_folder_path )
  {
    REQUIRE( ms::filesystem::is_exist_file( file_path ), " file should be exist." );
    REQUIRE( ms::filesystem::is_exist_folder( new_folder_path ), " new folder path should be exist." ); // �ݵ�� �����ϴ� �������� ��!

    const auto file_name     = ms::filesystem::extract_file_name( file_path );
    const auto new_file_path = new_folder_path.data() + file_name;

    std::filesystem::path old_p( file_path );
    std::filesystem::path new_p( new_file_path );
    std::filesystem::rename( old_p, new_p );
  }

  void remove_empty_folder( const std::string_view folder_path )
  {
    REQUIRE( is_folder_path( folder_path ), "folder_path should be end with /" );
    std::filesystem::path p( folder_path );
    std::filesystem::remove( p );
  }

  void remove_folder( const std::string_view folder_path )
  {
    REQUIRE( is_folder_path( folder_path ), "folder_path should be end with /" );
    std::filesystem::path p( folder_path );
    std::filesystem::remove_all( p );
  }

  void rename( const std::string& folder_path, const std::string& old_name, const std::string& new_name )
  {
    REQUIRE( is_folder_path( folder_path ), "folder_path should be end with /" );
    std::filesystem::rename( folder_path + old_name, folder_path + new_name );
  }

  void replace_file_content( const std::string_view file_path, const std::string_view old_content, const std::string_view new_contnet )
  {
    std::ifstream ifs( file_path.data() );
    REQUIRE( ifs.is_open(), "file should be exist" );

    ifs.seekg( 0, std::ios::end );
    const auto file_size = ifs.tellg();
    ifs.seekg( 0 );

    std::string new_contents;
    new_contents.reserve( file_size );

    std::string temp;
    while ( std::getline( ifs, temp ) )
    {
      ms::string::replace_inplace( temp, old_content, new_contnet );
      new_contents += temp + "\n";
      temp.clear();
    }

    std::ofstream of( file_path.data() );
    of << new_contents;
  }

  void remove_file( const std::string_view file_path )
  {
    REQUIRE( ms::filesystem::is_exist_file( file_path ), std::string( file_path ) + " is not exist file" );
    std::filesystem::path p( file_path );
    std::filesystem::remove( p );
  }

} // namespace ms::filesystem