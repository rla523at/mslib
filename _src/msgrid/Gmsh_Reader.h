#pragma once
#include "Reader.h"

#include <fstream>
#include <map>
#include <vector>

// class declaration
namespace ms::grid
{

class Gmsh_Reader : public Reader
{
public:
  Gmsh_Reader(const int dimension) : _dimension(dimension){};

public:
  Data read(const std::string_view file_path) const override;

private:
  void                       extract_block_type_node_data(std::ifstream& file, Data& data) const;
  void                       extract_element_datas(std::ifstream& file, Data& data) const;
  std::map<int, std::string> extract_physical_group_index_to_name(std::ifstream& file) const;

private:
  int _dimension;
};

} // namespace ms::grid