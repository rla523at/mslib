#pragma once

#include <fstream>
#include <string_view>
#include <vector>

namespace ms::tecplot
{

class Binary_File
{
public:
  Binary_File(const std::string_view file_path);
  Binary_File(const std::string_view file_path, std::ios_base::openmode mode);

public: // Command
  template <typename T>
  Binary_File& operator<<(const T value)
  {
    this->_file.write(reinterpret_cast<const char*>(&value), sizeof(T));
    return *this;
  };
  template <typename T>
  Binary_File& operator<<(const std::vector<T>& values)
  {
    for (const auto value : values)
    {
      this->_file.write(reinterpret_cast<const char*>(&value), sizeof(T));
    }
    return *this;
  };
  template <>
  Binary_File& operator<<(const char* value);
  template <>
  Binary_File& operator<<(const std::string& str);

  void close(void);

private:
  std::ofstream _file;
};

} // namespace ms::tecplot
