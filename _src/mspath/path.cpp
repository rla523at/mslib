#include "path.h"

#include "msexception/Exception.h"

#include <filesystem>

namespace ms::path
{

std::string extract_folder_path(const std::string_view file_path)
{
  return std::filesystem::path(file_path).parent_path().string() + "/";
}

std::vector<std::string> file_names_in_folder(const std::string_view folder_path)
{
  REQUIRE(is_folder_path(folder_path), "folder_path should be end with /");

  std::vector<std::string> file_names;

  std::filesystem::directory_iterator iter(folder_path);
  while (iter != std::filesystem::end(iter))
  {
    const auto& entry = *iter;

    if (!entry.is_directory())
    {
      file_names.push_back(entry.path().filename().string());
    }

    iter++;
  }

  return file_names;
}

std::vector<std::string> folder_names_in_folder(const std::string_view folder_path)
{
  REQUIRE(is_folder_path(folder_path), "folder_path should be end with /");

  std::vector<std::string> folder_names;

  std::filesystem::directory_iterator iter(folder_path);
  while (iter != std::filesystem::end(iter))
  {
    const auto& entry = *iter;

    if (entry.is_directory())
    {
      folder_names.push_back(entry.path().filename().string());
    }

    iter++;
  }

  return folder_names;
}

std::vector<std::string> file_paths_in_folder(const std::string& folder_path)
{
  REQUIRE(is_folder_path(folder_path), "folder_path should be end with /");

  std::vector<std::string> file_name_text;

  std::filesystem::directory_iterator iter(folder_path);
  while (iter != std::filesystem::end(iter))
  {
    const auto& entry = *iter;

    if (entry.is_directory())
    {
      iter++;
      continue;
    }

    file_name_text.push_back(entry.path().string());
    iter++;
  }

  return file_name_text;
}

bool is_exist_folder(const std::string_view folder_path)
{
  REQUIRE(is_folder_path(folder_path), "folder_path should be end with /");

  std::filesystem::path p(folder_path);
  return std::filesystem::exists(p);
}

bool is_exist_file(const std::string_view file_path)
{
  REQUIRE(!is_folder_path(file_path), "file path should not be end with /");

  std::filesystem::path p(file_path);
  return std::filesystem::exists(p);
}

bool is_folder_path(const std::string_view folder_path)
{
  const auto last_char = folder_path.back();
  return last_char == '/';
}

void make_folder(const std::string_view folder_path)
{
  REQUIRE(is_folder_path(folder_path), "folder_path should be end with /");

  if (folder_path.empty()) return;

  if (ms::path::is_exist_folder(folder_path)) return;

  std::filesystem::path p(folder_path);
  std::filesystem::create_directories(p);
}

void remove_empty_folder(const std::string_view folder_path)
{
  REQUIRE(is_folder_path(folder_path), "folder_path should be end with /");
  std::filesystem::path p(folder_path);
  std::filesystem::remove(p);
}

void remove_folder(const std::string_view folder_path)
{
  REQUIRE(is_folder_path(folder_path), "folder_path should be end with /");
  std::filesystem::path p(folder_path);
  std::filesystem::remove_all(p);
}

void rename(const std::string& folder_path, const std::string& old_name, const std::string& new_name)
{
  REQUIRE(is_folder_path(folder_path), "folder_path should be end with /");
  std::filesystem::rename(folder_path + old_name, folder_path + new_name);
}

} // namespace ms::path