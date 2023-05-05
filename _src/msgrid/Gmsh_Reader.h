#pragma once
#include "Reader.h"

#include <fstream>
#include <map>
#include <vector>

// forward declaration
namespace ms::grid
{
struct Nodes_Data;
struct Element_Data;
enum class Element_Type;
} // namespace ms::grid::reader

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
  Nodes_Data                  extract_block_type_node_data(std::ifstream& file) const;
  std::vector<Element_Data>   extract_element_datas(std::ifstream& file) const;
  std::map<int, Element_Type> extract_physical_group_index_to_element_type(std::ifstream& file) const;

private:
  int _dimension;
};

} // namespace ms::grid::reader