#pragma once
#include <string>
#include <string_view>
#include <vector>

// convention
// denominator in path should be '/'
// folder path always ended with '/'
// file/folder paths are not case-sensitive
namespace ms::path
{

std::string              extract_folder_path(const std::string_view file_path);
std::vector<std::string> file_names_in_folder(const std::string_view folder_path);
std::vector<std::string> folder_names_in_folder(const std::string_view folder_path);
std::vector<std::string> file_paths_in_folder(const std::string& path);
bool                     is_folder_path(const std::string_view folder_path);
bool                     is_exist_folder(const std::string_view folder_path);
bool                     is_exist_file(const std::string_view file_path);
void                     make_folder(const std::string_view folder_path);
void                     rename(const std::string& path, const std::string& old_name, const std::string& new_name);
void                     remove_empty_folder(const std::string_view folder_path);
void                     remove_folder(const std::string_view folder_path);

} // namespace ms::path