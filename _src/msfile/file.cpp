#include "file.h"

#include "msexception/Exception.h"
#include "msstring/string.h"
#include "mspath/path.h"
#include <fstream>
#include <filesystem>


namespace ms::file
{
	void replace_file_content(const std::string_view file_path, const std::string_view old_content, const std::string_view new_contnet)
	{
    std::ifstream ifs(file_path.data());
    REQUIRE(ifs.is_open(), "file should be exist");

    ifs.seekg(0, std::ios::end);
    const auto file_size = ifs.tellg();
    ifs.seekg(0);

    std::string new_contents;
    new_contents.reserve(file_size);

    std::string temp;
    while (std::getline(ifs, temp))
    {
      ms::string::replace_inplace(temp, old_content, new_contnet);
      new_contents += temp + "\n";
      temp.clear();
    }

    std::ofstream of(file_path.data());
    of << new_contents;
	}

  void remove_file(const std::string_view file_path)
  {
    REQUIRE(ms::path::is_exist_file(file_path), std::string(file_path) + " is not exist file");
    std::filesystem::path p(file_path);
    std::filesystem::remove(p);
  }
}