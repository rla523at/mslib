#pragma once
#include "Reader.h"

#include <fstream>
#include <map>
#include <vector>

// forward declaration
namespace ms::grid
{
struct Gmsh_Nodes_Data;
struct Gmsh_Elements_Data;
struct Gmsh_Physical_Data;
} // namespace ms::grid

/*




*/

// class declaration
namespace ms::grid
{

class Gmsh_Reader : public Reader
{
public:
  Gmsh_Reader(const int dimension) : _dimension(dimension){};

public:
  Grid_Data read(const std::string_view file_path) const override;

private:
  Grid_Data convert(Gmsh_Nodes_Data&& node_data, Gmsh_Elements_Data&& elem_data, const Gmsh_Physical_Data& phys_data) const;
  void read_node_data(std::ifstream& file, Gmsh_Nodes_Data& data) const;
  void read_elem_data(std::ifstream& file, Gmsh_Elements_Data& data) const;
  void read_phys_data(std::ifstream& file, Gmsh_Physical_Data& data) const;

//private:
//  void                       extract_block_type_node_data(std::ifstream& file, Grid_Data& data) const;
//  void                       extract_element_datas(std::ifstream& file, Grid_Data& data) const;
//  std::map<int, std::string> extract_physical_group_index_to_name(std::ifstream& file) const;

private:
  int _dimension;
};

} // namespace ms::grid