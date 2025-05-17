#pragma once
/*
convention

* denominator in path should be '/'
* folder path always ended with '/'
* file/folder paths are not case-sensitive
* base_name of "folder1/abc.txt" --> "abc"
* file name of "folder1/abc.txt" --> "abc.txt"
* file path of "folder1/abc.txt" --> "folder1/"
* path : A/B/C/D/E/
  * E : depth 0
  * D : depth 1
  * C : depth 2
  ...

*/

namespace mslib::filesystem
{
  std::wstring_view        base_name_view( std::wstring_view file_path );
  void                     copy_file( const std::string_view from_file_path, const std::string_view to_file_path, const std::filesystem::copy_options option = std::filesystem::copy_options::none );
  std::wstring             excutable_file_path_utf16( void );
  std::string              excutable_file_path( void );
  std::wstring             extract_folder_path( const std::wstring_view file_path );
  std::string              extract_folder_path( const std::string_view file_path );
  std::wstring             extract_file_name( const std::wstring_view file_path );
  std::string              extract_file_name( const std::string_view file_path );
  std::wstring_view        file_name_view( std::wstring_view file_path );
  std::vector<std::string> file_names_in_folder( const std::string_view folder_path );
  std::wstring_view        folder_path_view( std::wstring_view file_path );
  std::vector<std::string> folder_names_in_folder( const std::string_view folder_path );
  std::vector<std::string> file_paths_in_folder( const std::string& path );
  bool                     has_this_extension( const std::string_view file_path, const std::string_view extension );
  bool                     is_folder_path( const std::wstring_view folder_path );
  bool                     is_folder_path( const std::string_view folder_path );
  bool                     is_exist_folder( const std::wstring_view folder_path );
  bool                     is_exist_folder( const std::string_view folder_path );
  bool                     is_exist_file( const std::wstring_view file_path );
  bool                     is_exist_file( const std::string_view file_path );
  bool                     is_modified_later( const std::wstring_view reference_file_path, const std::wstring_view target_file_path );
  bool                     is_modified_later( const std::string_view reference_file_path, const std::string_view target_file_path );
  void                     make_folder( const std::string_view folder_path );
  void                     move_file( const std::string_view file_path, const std::string_view new_folder_path );
  std::wstring_view        parent_path_view( const std::wstring_view path, const int depth );
  std::string_view         parent_path_view( const std::string_view path, const int depth );
  void                     replace_file( const std::string_view from_file_path, const std::string_view to_file_path );
  void                     rename( const std::string& path, const std::string& old_name, const std::string& new_name );
  void                     remove_empty_folder( const std::string_view folder_path );
  void                     remove_folder( const std::string_view folder_path );
  void                     replace_file_content( const std::string_view file_path, const std::string_view old_content, const std::string_view new_contnet );
  void                     remove_file( const std::string_view file_path );

} // namespace mslib::filesystem
