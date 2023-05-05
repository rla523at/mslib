#include "Binary_Writer.h"

#include "msexception/Exception.h"
#include <iostream>

namespace ms::tecplot
{

Binary_File::Binary_File(const std::string_view file_path)
{
  this->_file.open(file_path.data(), std::ios::binary);
  REQUIRE(this->_file.is_open(), "file should be opened");
}

Binary_File::Binary_File(const std::string_view file_path, std::ios_base::openmode mode)
{
  this->_file.open(file_path.data(), std::ios::binary | mode);
  REQUIRE(this->_file.is_open(), "file should be opened");
}

template <>
Binary_File& Binary_File::operator<<(const char* value)
{
  this->_file.write(reinterpret_cast<const char*>(value), sizeof(value));
  return *this;
}

template <>
Binary_File& Binary_File::operator<<(const std::string& str)
{
  this->_file << str;
  return *this;
}

void Binary_File::close(void)
{
  this->_file.close();
}

} // namespace ms::tecplot