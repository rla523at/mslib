// #include "Data.h"
//
// #include "msexception/Exception.h"
//
// namespace ms::grid
//{
//
// const double* Grid_Nodes_Data::coordinates_ptr(const int node_index) const
//{
//   //Asuume that node index in Grid File start with 1
//   REQUIRE(1 <= node_index, "node index should be greater than 1");
//
//   auto index = node_index - 1;
//
//   if (this->type == Coordinate_Type::BLOCK)
//   {
//     return this->coordinates.data() + index;
//   }
//   else if (this->type == Coordinate_Type::NODAL)
//   {
//     return this->coordinates.data() + this->dimension * index;
//   }
//   else
//   {
//     EXCEPTION("unsupported coordinate type");
//     return nullptr;
//   }
// }
//
// int Grid_Nodes_Data::stride(void) const
//{
//   if (this->type == Coordinate_Type::BLOCK)
//   {
//     return this->num_nodes;
//   }
//   else if (this->type == Coordinate_Type::NODAL)
//   {
//     return 1;
//   }
//   else
//   {
//     EXCEPTION("unsupported coordinate type");
//     return -1;
//   }
// }
//
//
// } // namespace ms
