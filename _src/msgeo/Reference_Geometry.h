#pragma once
#include <vector>

// forward declaration
namespace ms::sym
{
class Polynomial;
class Symbol;
} // namespace ms::sym

namespace ms::geo
{
class Node_Const_Wrapper;
class Nodes_Const_Wrapper;
using Polynomial  = ms::sym::Polynomial;
using Symbol      = ms::sym::Symbol;
using Polynomials = std::vector<ms::sym::Polynomial>;
using Symbols     = std::vector<ms::sym::Symbol>;

} // namespace ms::geo

/*





*/

// class declaration
namespace ms::geo
{

class Reference_Geometry
{
protected:
  Reference_Geometry(void)          = default;
  virtual ~Reference_Geometry(void) = default;

public:
  virtual Polynomials         cal_normal_functions(const Polynomials& parametric_functions) const                    = 0;
  virtual Polynomials         cal_parametric_functions(const std::vector<Node_Const_Wrapper>& consisting_nodes) const = 0;
  virtual Node_Const_Wrapper  center_point(void) const                                                               = 0;
  virtual Nodes_Const_Wrapper quadrature_points(const int integrand_degree) const                                     = 0;
};

} // namespace ms::geo

// free functions
namespace ms::geo
{

Symbol      cal_curvature(const Polynomials& curve);
Symbols     cal_principal_normal(const Polynomials& curve);
Polynomials cal_paramteric_curve_normal_functions(const Polynomials& curve);

} // namespace ms::geo

// const std::vector<Euclidean_Vector>& get_post_points(
//     const int post_order) const;
// const std::vector<std::vector<uint>>& get_connectivities(
//     const int post_order) const;

// std::vector<int> vertex_node_index_sequneces(void) const;

// virtual std::vector<std::shared_ptr<const Reference_Geometry>>
// face_reference_geometries(void) const = 0;
// virtual Figure figure(void) const = 0;
// virtual bool is_simplex(void) const = 0;
// virtual bool is_line(void) const = 0;
// virtual int num_vertices(void) const = 0;
// virtual int num_post_nodes(const int post_order) const = 0;
// virtual int num_post_elements(const int post_order) const = 0;
// virtual Vector_Function<Polynomial> make_normal_vector_function(
//     const Vector_Function<Polynomial>& mapping_function) const = 0;
// virtual Euclidean_Vector random_point(void) const = 0;
// virtual std::vector<std::vector<int>> set_of_face_vertex_index_sequences(
//     void) const = 0;
// virtual std::vector<std::vector<int>> set_of_face_node_index_sequences(
//     void) const = 0;
// virtual std::vector<std::shared_ptr<const Reference_Geometry>>
// sub_simplex_reference_geometries(void) const = 0;
// virtual std::vector<std::vector<int>>
// set_of_sub_simplex_vertex_index_sequences(void) const = 0;
// virtual Irrational_Function scale_function(
//     const Vector_Function<Polynomial>& mapping_function) const = 0;
// virtual int scale_function_order(void) const = 0;

// Matrix make_inverse_mapping_monomial_matrix(void) const;
// std::vector<std::vector<uint>> quadrilateral_connectivities(
//     const std::array<uint, 4>& node_indexes) const;
//  std::vector<std::vector<uint>> sliced_hexahedral_connectivities(const
//  std::array<uint, 7>& node_indexes) const; std::vector<std::vector<uint>>
//  hexahedral_connectivities(const std::array<uint, 8>& node_indexes) const;

// virtual Quadrature_Rule make_quadrature_rule(
//     const int integrand_order) const = 0;
// virtual std::vector<Euclidean_Vector> make_mapping_points(
//     void) const = 0;
// virtual std::vector<Euclidean_Vector> make_post_points(
//     const int post_order) const = 0;
// virtual std::vector<std::vector<uint>> make_connectivities(
//     const int post_order) const = 0;

//// tag to parametric function basis matrix values
// mutable std::map<int, std::vector<double>> tag_to_SFM_vals_;
//// tag to parametric function bases
// mutable std::map<int, std::vector<Polynomial>> tag_to_PF_bases_;

// mutable std::map<int, std::vector<std::vector<uint>>>
// post_order_to_connectivities_;

// template <typename Function>
// Vector_Function<Function> operator*(
//     const Matrix& matrix, const Vector_Function<Function>& vector_function) {
//   const auto [num_row, num_column] = matrix.size();
//
//   std::vector<Function> functions(num_row);
//
//   for (size_t i = 0; i < num_row; ++i)
//     for (size_t j = 0; j < num_column; ++j)
//       functions[i] += matrix.at(i, j) * vector_function[j];
//
//   return functions;
// }

// const std::vector<Euclidean_Vector>&
// Reference_Geometry::get_post_points(const ushort post_order) const
//{
//	if (!this->post_order_to_post_points_.contains(post_order))
//	{
//		this->post_order_to_post_points_.emplace(post_order,
// this->make_post_points(post_order));
//	}
//
//	return this->post_order_to_post_points_.at(post_order);
// }
//
// const std::vector<std::vector<uint>>&
// Reference_Geometry::get_connectivities(const ushort post_order) const
//{
//	if (!this->post_order_to_connectivities_.contains(post_order))
//	{
//		this->post_order_to_connectivities_.emplace(post_order,
// this->make_connectivities(post_order));
//	}
//
//	return this->post_order_to_connectivities_.at(post_order);
// }
//

// std::vector<ushort> Reference_Geometry::vertex_node_index_sequneces(void)
// const
//{
//	const auto num_vertices = this->num_vertices();
//	std::vector<ushort> vnode_index_orders(num_vertices);
//
//	for (ushort i = 0; i < num_vertices; ++i)
//		vnode_index_orders[i] = i;
//
//	return vnode_index_orders;
// }

//
// std::vector<std::vector<uint>>
// Reference_Geometry::quadrilateral_connectivities(const std::array<uint, 4>&
// node_indexes) const
//{
//	//   2式式式式式3
//	//   弛     弛
//	//   0式式式式式1
//	return { { node_indexes[0],node_indexes[1],node_indexes[2] },
//			 { node_indexes[1],node_indexes[3],node_indexes[2] } };
// };