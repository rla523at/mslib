#pragma once
#include <string>
#include <string_view>
#include <vector>
#include <filesystem>

// convention
// denominator in path should be '/'
// folder path always ended with '/'
// file/folder paths are not case-sensitive
namespace ms::filesystem
{
void										 copy_file(const std::string_view from_file_path, const std::string_view to_file_path, const std::filesystem::copy_options option = std::filesystem::copy_options::none);
std::string              extract_folder_path(const std::string_view file_path);
std::string              extract_file_name(const std::string_view file_path);
std::vector<std::string> file_names_in_folder(const std::string_view folder_path);
std::vector<std::string> folder_names_in_folder(const std::string_view folder_path);
std::vector<std::string> file_paths_in_folder(const std::string& path);
bool                     has_this_extension(const std::string_view file_path, const std::string_view extension);
bool                     is_folder_path(const std::string_view folder_path);
bool                     is_exist_folder(const std::string_view folder_path);
bool                     is_exist_file(const std::string_view file_path);
void                     make_folder(const std::string_view folder_path);
void                     rename(const std::string& path, const std::string& old_name, const std::string& new_name);
void                     remove_empty_folder(const std::string_view folder_path);
void                     remove_folder(const std::string_view folder_path);
void                     replace_file_content(const std::string_view file_path, const std::string_view old_content, const std::string_view new_contnet);
void                     remove_file(const std::string_view file_path);

} // namespace ms::filesystem