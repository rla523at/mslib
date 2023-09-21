#pragma once
#include "Reference_Geometry.h"

#include <array>
#include <map>

// forward declaration
namespace ms::geo
{
class Reference_Geometry_Container;
}

/*





*/

// class declaration
namespace ms::geo
{

class Reference_Point : public Reference_Geometry
{
  // overriding methods from Reference_Geometry
public:
  ms::sym::Polynomials          cal_normal_functions(const ms::sym::Polynomials& parametric_functions) const override;
  ms::sym::Polynomials          cal_parametric_functions(const std::vector<Node_View>& consisting_nodes) const override;
  int                           cal_parameter_order(const int num_points) const override;
  ms::sym::Symbol               cal_scale_function(const ms::sym::Polynomials& parametric_functions) const override;
  Node_View            center_point(void) const override;
  int                           dimension(void) const override;
  Figure                        face_figure(const int face_index) const override;
  std::vector<std::vector<int>> face_index_to_face_vnode_indexes(void) const override;
  const std::vector<double>&    get_quadrature_weights(const int integrand_degree) const override;
  const Partition_Data&         get_partition_data(const int partition_order) const override;
  bool                          is_valid_num_points(const int num_points) const override;
  bool                          is_point(void) const override;
  bool                          is_line(void) const override;
  std::vector<int>              node_indexes(const int parameter_order) const override;
  int                           num_faces(void) const override;
  int                           num_vertices(void) const override;
  Nodes_View           quadrature_points(const int integrand_degree) const override;

private:
  static constexpr std::array<double, 1> _center_coords = {0.0};

private:
  friend Reference_Geometry_Container;

private:
  // All Reference Geometry objects are created only once in the Container class.
  // To prevent unnecessary creation, acess to the constructor has been blocked
  Reference_Point(void) = default;
};

class Reference_Geometry_Common : public Reference_Geometry
{
  // overriding methods from Reference_Geometry
public:
  ms::sym::Polynomials       cal_parametric_functions(const std::vector<Node_View>& consisting_nodes) const override;
  const std::vector<double>& get_quadrature_weights(const int integrand_degree) const override;
  const Partition_Data&      get_partition_data(const int partition_order) const override;
  bool                       is_point(void) const override;
  Nodes_View        quadrature_points(const int integrand_degree) const override;

protected:
  virtual int                  cal_quadrature_rule_tag(const int integrand_degree) const        = 0;
  virtual void                 create_and_store_partition_data(const int partition_order) const = 0;
  virtual void                 create_and_store_quadrature_points(const int tag) const          = 0;
  virtual void                 create_and_store_quadrature_weights(const int tag) const         = 0;
  virtual int                  num_quadrature_points(const int tag) const                       = 0;
  virtual int                  num_parametric_function_reference_points(const int tag) const    = 0;
  virtual ms::sym::Polynomials make_parametric_function_bases(const int tag) const              = 0;
  virtual std::vector<double>  make_parametric_functions_reference_coords(const int tag) const  = 0;

private:
  const ms::sym::Polynomials& get_shape_functions(const int parameter_order) const;
  ms::sym::Polynomials        make_shape_functions(const int parameter_order) const;

protected:
  mutable std::map<int, Partition_Data>       _order_to_partition_data;
  mutable std::map<int, Nodes>                _tag_to_quadrature_points;
  mutable std::map<int, std::vector<double>>  _tag_to_quadrature_weights;
  mutable std::map<int, ms::sym::Polynomials> _parameter_order_to_shape_functions;

protected:
  virtual ~Reference_Geometry_Common(void) = default;
};

class Reference_Line : public Reference_Geometry_Common
{
  // overriding methods from Reference_Geometry
public:
  Node_View            center_point(void) const override;
  ms::sym::Polynomials          cal_normal_functions(const ms::sym::Polynomials& curve) const override;
  int                           cal_parameter_order(const int num_points) const override;
  ms::sym::Symbol               cal_scale_function(const ms::sym::Polynomials& parametric_functions) const override;
  int                           dimension(void) const override;
  Figure                        face_figure(const int face_index) const override;
  std::vector<std::vector<int>> face_index_to_face_vnode_indexes(void) const override;
  bool                          is_valid_num_points(const int num_points) const override;
  bool                          is_line(void) const override;
  std::vector<int>              node_indexes(const int parameter_order) const override;
  int                           num_faces(void) const override;
  int                           num_vertices(void) const override;

  // overriding methods from Reference_Geometry_Common
protected:
  int                  cal_quadrature_rule_tag(const int integrand_degree) const override;
  void                 create_and_store_partition_data(const int partition_order) const override;
  void                 create_and_store_quadrature_points(const int tag) const override;
  void                 create_and_store_quadrature_weights(const int tag) const override;
  int                  num_quadrature_points(const int tag) const override;
  int                  num_parametric_function_reference_points(const int param_order) const override;
  ms::sym::Polynomials make_parametric_function_bases(const int tag) const override;
  std::vector<double>  make_parametric_functions_reference_coords(const int tag) const override;

private:
  static constexpr std::array<double, 1> center_coords_ = {0.0};

private:
  friend Reference_Geometry_Container;

private:
  // All Reference Geometry objects are created only once in the Container class.
  // To prevent unnecessary creation, acess to the constructor has been blocked
  Reference_Line(void) = default;
};

} // namespace ms::geo

// std::vector<std::shared_ptr<const Reference_Geometry>> face_reference_geometries(void) const override;
// Figure figure(void) const override;
// bool is_simplex(void) const override;
// bool is_line(void) const override;
// int num_vertices(void) const override;
// int num_post_nodes(const int post_order) const override;
// int num_post_elements(const int post_order) const override;
// Vector_Function<Polynomial> make_normal_vector_function(
//     const Vector_Function<Polynomial>& mapping_function) const override;
// Euclidean_Vector random_point(void) const override;
// std::vector<std::vector<int>> set_of_face_vertex_index_sequences(
//     void) const override;
// std::vector<std::vector<int>> set_of_face_node_index_sequences(
//     void) const override;
// std::vector<std::shared_ptr<const Reference_Geometry>>
// sub_simplex_reference_geometries(void) const override;
// std::vector<std::vector<int>> set_of_sub_simplex_vertex_index_sequences(
//     void) const override;
// Irrational_Function scale_function(
//     const Vector_Function<Polynomial>& mapping_function) const override;
// int scale_function_order(void) const override;

// protected:
//  Quadrature_Rule make_quadrature_rule(
//      const int integrand_order) const override;
//  std::vector<Euclidean_Vector> make_mapping_points(void) const override;
//  std::vector<Euclidean_Vector> make_post_points(
//      const int post_order) const override;
//  std::vector<std::vector<uint>> make_connectivities(
//      const int post_order) const override;
//  Vector_Function<Polynomial> make_mapping_monomial_vector_function(
//      void) const override;
//};

//
// class Reference_Triangle : public Reference_Geometry {
// public:
//  Reference_Triangle(const int order);
//
// public:  // Query
//  Node_UPtr center_point(void) const override;
//  std::vector<std::shared_ptr<const Reference_Geometry>>
//  face_reference_geometries(void) const override;
//  Figure figure(void) const override;
//  bool is_simplex(void) const override;
//  bool is_line(void) const override;
//  int num_vertices(void) const override;
//  int num_post_nodes(const int post_order) const override;
//  int num_post_elements(const int post_order) const override;
//  Vector_Function<Polynomial> make_normal_vector_function(
//      const Vector_Function<Polynomial>& mapping_function)
//      const override;  // 2D Element ����
//  Euclidean_Vector random_point(void) const override;
//  std::vector<std::vector<int>> set_of_face_vertex_index_sequences(
//      void) const override;
//  std::vector<std::vector<int>> set_of_face_node_index_sequences(
//      void) const override;
//  std::vector<std::shared_ptr<const Reference_Geometry>>
//  sub_simplex_reference_geometries(void) const override;
//  std::vector<std::vector<int>> set_of_sub_simplex_vertex_index_sequences(
//      void) const override;
//  Irrational_Function scale_function(
//      const Vector_Function<Polynomial>& mapping_function)
//      const override;  // 2D Element ����
//  int scale_function_order(void) const override;
//
// protected:
//  std::vector<Euclidean_Vector> make_mapping_points(void) const override;
//  Quadrature_Rule make_quadrature_rule(
//      const int integrand_order) const override;
//  std::vector<Euclidean_Vector> make_post_points(
//      const int post_order) const override;
//  std::vector<std::vector<uint>> make_connectivities(
//      const int post_order) const override;
//  Vector_Function<Polynomial> make_mapping_monomial_vector_function(
//      void) const override;
//};
//
// class Reference_Quadrilateral : public Reference_Geometry {
// public:
//  Reference_Quadrilateral(const int order);
//
// public:  // Query
//  Node_UPtr center_point(void) const override;
//  std::vector<std::shared_ptr<const Reference_Geometry>>
//  face_reference_geometries(void) const override;
//  Figure figure(void) const override;
//  bool is_simplex(void) const override;
//  bool is_line(void) const override;
//  int num_vertices(void) const override;
//  int num_post_nodes(const int post_order) const override;
//  int num_post_elements(const int post_order) const override;
//  Vector_Function<Polynomial> make_normal_vector_function(
//      const Vector_Function<Polynomial>& mapping_function)
//      const override;  // 2D Element ����
//  Euclidean_Vector random_point(void) const override;
//  std::vector<std::vector<int>> set_of_face_vertex_index_sequences(
//      void) const override;
//  std::vector<std::vector<int>> set_of_face_node_index_sequences(
//      void) const override;
//  std::vector<std::shared_ptr<const Reference_Geometry>>
//  sub_simplex_reference_geometries(void) const override;
//  std::vector<std::vector<int>> set_of_sub_simplex_vertex_index_sequences(
//      void) const override;
//  Irrational_Function scale_function(
//      const Vector_Function<Polynomial>& mapping_function)
//      const override;  // 2D Element ����
//  int scale_function_order(void) const override;
//
// protected:
//  std::vector<Euclidean_Vector> make_mapping_points(void) const override;
//  Quadrature_Rule make_quadrature_rule(
//      const int integrand_order) const override;
//  std::vector<Euclidean_Vector> make_post_points(
//      const int post_order) const override;
//  std::vector<std::vector<uint>> make_connectivities(
//      const int post_order) const override;
//  Vector_Function<Polynomial> make_mapping_monomial_vector_function(
//      void) const override;
//};
//
// class Reference_Geometry_Container {
// public:
//  static const std::shared_ptr<const Reference_Geometry>& get_shared_ptr(
//      const Figure figure, const int order);
//
// private:
//  static void store_line(const int order);
//  static void store_triangle(const int order);
//  static void store_quadrilateral(const int order);
//
// private:
//  Reference_Geometry_Container(void) = delete;
//
// private:
//  static inline std::map<int, std::shared_ptr<const Reference_Geometry>>
//      order_to_reference_line_;
//  static inline std::map<int, std::shared_ptr<const Reference_Geometry>>
//      order_to_reference_triangle_;
//  static inline std::map<int, std::shared_ptr<const Reference_Geometry>>
//      order_to_reference_quadrilateral_;
//  static inline std::shared_ptr<const Reference_Geometry> error_ = nullptr;
//};

// Euclidean_Vector<space_dimension>
// ReferenceGeometry<space_dimension>::center_node(void) const {
// case
// Figure::tetrahedral:	return { -0.5, -0.5, -0.5 }; 		case
// Figure::hexahedral: return { 0, 0, 0 };
// }

// int ReferenceGeometry<space_dimension>::num_vertex(void) const {
//	case Figure::tetrahedral:	return 4;
//	case Figure::hexahedral:	return 8;
// }

// std::vector<std::vector<int>>
// ReferenceGeometry<space_dimension>::set_of_face_vertex_node_index_orders(void)
// const {
///	case Figure::tetrahedral: {
//		//   3
//		//   ��
//		//   0��������2
//		//  /
//		// 1
//		//index ������� transformation �Ǿ��� ��, cell �߽ɺθ�
// �ٶ󺸰� 		std::vector<int> face0_node_index = { 0,1,2 };
// std::vector<int> face1_node_index = { 0,2,3 }; 		std::vector<int>
// face2_node_index = { 0,3,1
//}; 		std::vector<int> face3_node_index = { 1,3,2 }; 		return {
// face0_node_index,face1_node_index, face2_node_index,face3_node_index };
//	}
//	case Figure::hexahedral: {
//		//    4����������7
//		//   /��    /��
//		//  5����������6 ��
//		//  �� 0����������3
//		//  ��/    ��/
//		//  1����������2
//		//index ������� transformation �Ǿ��� ��, cell �߽ɺθ�
// �ٶ󺸰� 		std::vector<int> face0_node_index = { 0,1,2,3 };
//		std::vector<int> face1_node_index = { 0,4,5,1 };
//		std::vector<int> face2_node_index = { 1,5,6,2 };
//		std::vector<int> face3_node_index = { 2,6,7,3 };
//		std::vector<int> face4_node_index = { 0,3,7,4 };
//		std::vector<int> face5_node_index = { 4,7,6,5 };
//
//		return { face0_node_index, face1_node_index, face2_node_index,
// face3_node_index, face4_node_index, face5_node_index };
//	}
//	default:
//		throw std::runtime_error("not supported element figure");
//		return std::vector<std::vector<int>>();
//	}
//}

// std::vector<std::vector<int>>
// ReferenceGeometry<space_dimension>::set_of_face_node_index_orders(void) const
// { 	case Figure::tetrahedral: {
//		//   3
//		//   ��
//		//   0��������2
//		//  /
//		// 1
//		//index ������� transformation �Ǿ��� ��, cell �߽ɺθ�
// �ٶ󺸰� 		dynamic_require(this->figure_order_ == 1, "this figure
// does not support high order mesh yet"); 		std::vector<int>
// face0_node_index = { 0,1,2
//}; 		std::vector<int> face1_node_index = { 0,2,3 };
// std::vector<int> face2_node_index = { 0,3,1 }; 		std::vector<int>
// face3_node_index = { 1,3,2
//};
//
//		return { face0_node_index,face1_node_index,
// face2_node_index,face3_node_index };
//	}
//	case Figure::hexahedral: {
//		//    4����������7
//		//   /��    /��
//		//  5����������6 ��
//		//  �� 0����������3
//		//  ��/    ��/
//		//  1����������2
//		//index ������� transformation �Ǿ��� ��, cell �߽ɺθ�
// �ٶ󺸰� 		dynamic_require(this->figure_order_ == 1, "this figure
// does not support high order mesh yet"); 		std::vector<int>
// face0_node_index = { 0,1,2,3 }; 		std::vector<int>
// face1_node_index = { 0,4,5,1 }; 		std::vector<int>
// face2_node_index = { 1,5,6,2
//}; 		std::vector<int> face3_node_index = { 2,6,7,3 };
// std::vector<int> face4_node_index = { 0,3,7,4 }; 		std::vector<int>
// face5_node_index = { 4,7,6,5
//};
//
//		return { face0_node_index, face1_node_index, face2_node_index,
// face3_node_index, face4_node_index, face5_node_index };
//	}
// }

// std::vector<ReferenceGeometry<space_dimension>>
// ReferenceGeometry<space_dimension>::face_reference_geometries(void) const {
//	case Figure::tetrahedral: {
//		//   3
//		//   ��
//		//   0��������2
//		//  /
//		// 1
//		const ReferenceGeometry face0_refrence_geometry = {
// Figure::triangle, this->figure_order_ }; 		const ReferenceGeometry
// face1_refrence_geometry = { Figure::triangle, this->figure_order_ };
// const ReferenceGeometry face2_refrence_geometry = { Figure::triangle,
// this->figure_order_ }; 		const ReferenceGeometry
// face3_refrence_geometry = {
// Figure::triangle, this->figure_order_ }; 		return {
// face0_refrence_geometry,face1_refrence_geometry,
// face2_refrence_geometry,face3_refrence_geometry };
//	}
//	case Figure::hexahedral: {
//		//    4����������7
//		//   /��    /��
//		//  5����������6 ��
//		//  �� 0����������3
//		//  ��/    ��/
//		//  1����������2
//		const ReferenceGeometry face0_reference_geometry = {
// Figure::quadrilateral, this->figure_order_ }; 		const
// ReferenceGeometry face1_reference_geometry = { Figure::quadrilateral,
// this->figure_order_ }; 		const ReferenceGeometry
// face2_reference_geometry = { Figure::quadrilateral, this->figure_order_ };
// const ReferenceGeometry face3_reference_geometry = { Figure::quadrilateral,
// this->figure_order_ }; 		const ReferenceGeometry
// face4_reference_geometry = { Figure::quadrilateral, this->figure_order_ };
// const ReferenceGeometry face5_reference_geometry = { Figure::quadrilateral,
// this->figure_order_ };
//
//		return { face0_reference_geometry, face1_reference_geometry,
// face2_reference_geometry, face3_reference_geometry, face4_reference_geometry,
// face5_reference_geometry };
//	}
// }

// Irrational_Function<space_dimension>
// ReferenceGeometry<space_dimension>::scale_function(const
// Vector_Function<Polynomial<space_dimension>, space_dimension>&
// mapping_function) const { 			case Figure::tetrahedral:
// case Figure::hexahedral:
//{ 				constexpr int r = 0;
// constexpr int s = 1; 				constexpr int t = 2;
// const auto mf_r = mapping_function.differentiate(r);
// const
// auto mf_s =
// mapping_function.differentiate(s); 				const auto mf_t
// = mapping_function.differentiate(t);
//
//				return ms::scalar_triple_product(mf_r, mf_s,
// mf_t).be_absolute();
// }

// std::vector<ReferenceGeometry<space_dimension>>
// ReferenceGeometry<space_dimension>::sub_simplex_reference_geometries(void)
// const { 	dynamic_require(!this->is_simplex(), "This is routine for
// non-simplex
// figure"); 	dynamic_require(this->figure_order_ == 1, "This is only valid
// for P1 figure");
//
//		case Figure::hexahedral: {
//			//    4����������7
//			//   /��    /��
//			//  5����������6 ��
//			//  �� 0����������3
//			//  ��/    ��/
//			//  1����������2
//
//			ReferenceGeometry<space_dimension>
// simplex1(Figure::tetrahedral, this->figure_order_);
//			ReferenceGeometry<space_dimension>
// simplex2(Figure::tetrahedral, this->figure_order_);
//			ReferenceGeometry<space_dimension>
// simplex3(Figure::tetrahedral, this->figure_order_);
//			ReferenceGeometry<space_dimension>
// simplex4(Figure::tetrahedral, this->figure_order_);
//			ReferenceGeometry<space_dimension>
// simplex5(Figure::tetrahedral, this->figure_order_);
//			ReferenceGeometry<space_dimension>
// simplex6(Figure::tetrahedral, this->figure_order_);
//			ReferenceGeometry<space_dimension>
// simplex7(Figure::tetrahedral, this->figure_order_);
//			ReferenceGeometry<space_dimension>
// simplex8(Figure::tetrahedral, this->figure_order_);
//
//			return { simplex1, simplex2, simplex3, simplex4,
// simplex5, simplex6, simplex7, simplex8 };
//		}
//	default:
//		throw std::runtime_error("not supported figure");
//	}
//
// }

// std::vector<std::vector<int>>
// ReferenceGeometry<space_dimension>::set_of_sub_simplex_vertex_node_index_orders(void)
// const { 	dynamic_require(!this->is_simplex(), "This is routine for
// non-simplex
// figure"); 	dynamic_require(this->figure_order_ == 1, "This is only valid
// for P1 figure");
//
//		case Figure::hexahedral: {
//			//    4����������7
//			//   /��    /��
//			//  5����������6 ��
//			//  �� 0����������3
//			//  ��/    ��/
//			//  1����������2
//
//			const std::vector<int> simplex1 = { 0,1,3,4 };
//			const std::vector<int> simplex2 = { 1,2,0,5 };
//			const std::vector<int> simplex3 = { 2,3,1,6 };
//			const std::vector<int> simplex4 = { 3,0,2,7 };
//			const std::vector<int> simplex5 = { 4,5,7,0 };
//			const std::vector<int> simplex6 = { 5,6,4,1 };
//			const std::vector<int> simplex7 = { 6,7,5,2 };
//			const std::vector<int> simplex8 = { 7,4,6,3 };
//
//			return { simplex1, simplex2, simplex3, simplex4,
// simplex5, simplex6, simplex7, simplex8 };
//		}
//	default:
//		throw std::runtime_error("not supported figure");
//	}
//
// }

// std::vector<Euclidean_Vector<space_dimension>>
// ReferenceGeometry<space_dimension>::mapping_nodes(void) const {
// case
// Figure::tetrahedral: {
//				//   3
//				//   ��
//				//   0��������2
//				//  /
//				// 1
//				switch (this->figure_order_) {
//					case 1:		return { { -1, -1, -1 },
//{ 1, -1, -1 }, { -1, 1, -1 }, { -1, -1, 1 } };
//					case 2:		return { { -1, -1, -1 },
//{ 1, -1, -1 }, { -1, 1, -1 }, { -1, -1, 1 }, { 0, -1, -1 }, { -1, 0, -1 }, {
// -1, -1, 0 }, { 0, 0, -1 }, { 0, -1, 0 }, { -1, 0, 0 } }; default:	throw
// std::runtime_error("unsuported figure order");
//				}
//				break;
//			}
//			case Figure::hexahedral: {
//				//	     4����������7
//				//      /��    /��
//				//     5����������6 ��
//				//     �� 0����������3
//				//	   ��/    ��/
//				//     1����������2
//
//				switch (this->figure_order_) {
//					case 1:		return { { -1, -1, -1 },
//{ 1, -1, -1 }, { 1, 1, -1 }, { -1, 1, -1 }, {-1, -1, 1}, {1, -1, 1}, {1, 1,
// 1},
//{-1, 1, 1} }; 					default:	throw
// std::runtime_error("unsuported figure order");
//				}
//				break;
//			}
// }

// Dynamic_Vector_Function<Polynomial<space_dimension>>
// ReferenceGeometry<space_dimension>::mapping_monomial_vector_function(void)
// const {

//			case Figure::tetrahedral: {
//				const auto n = this->figure_order_;
//				const auto num_monomial = static_cast<int>(1.5 *
//(n
//+ 2) * (n - 1) + 4);
// std::vector<Polynomial<space_dimension>>
// mapping_monomial_vector(num_monomial);
//
//				int index = 0;
//				for (int a = 0; a <= n; ++a) {
//					for (int b = 0; b <= a; ++b)
//						for (int c = 0; c <= b; ++c)
//							mapping_monomial_vector[index++]
//= (r ^ (a - b)) * (s ^ (b - c)) * (t ^ c);
//				}
//
//				return mapping_monomial_vector; // 1 r s t r^2
// rs rt s^2 st t^2 ...
//			}
//			case Figure::hexahedral: {
//				const auto n = this->figure_order_;
//				const auto num_monomial = (n + 1) * (n + 1) * (n
//+ 1); 				std::vector<Polynomial<space_dimension>>
// mapping_monomial_vector(num_monomial);
//
//				int index = 0;
//
//				//make bottom monomials
//				//same process with quadrilateral monomials
//				for (int a = 0; a <= n; ++a) {
//					for (int b = 0; b <= a; ++b)
//						mapping_monomial_vector[index++]
//= (r ^ a)
//* (s ^ b);
//
//					if (a == 0)
//						continue;
//
//					for (int c = static_cast<int>(a - 1); 0
//<= c;
//--c) 						mapping_monomial_vector[index++]
//= (r ^ c) * (s ^ a);
//				}
//
//				//make higher floor monomials by multiplying t
//				const auto num_bottom_monomial = (n + 1) * (n +
// 1); 				const auto num_floor = n;
//
//				for (int i = 1; i <= num_floor; ++i) {
//					for (int j = 0; j < num_bottom_monomial;
//++j) 						mapping_monomial_vector[index++]
//= mapping_monomial_vector[j] * (t ^ i);
//				}
//
//				return mapping_monomial_vector; // {1 r s rs}
//{...}
//*
// t
//												//
//{1 r rs s r^2 r^2s r^2s^2 rs^2 s^2} {...} * t  {...} * t^2
//			}
//}

// Quadrature_Rule<space_dimension>
// ReferenceGeometry<space_dimension>::reference_quadrature_rule(const int
// integrand_order) const {

//			case Figure::tetrahedral: {
//				switch (integrand_order){
//					case 0:
//					case 1:		return { { { -0.5, -0.5,
//-0.5 } } , { 1.3333333333333333333333333333333333333 } };
// case 2: return { { { -0.72360679774997896964091736687312762354,
//-0.72360679774997896964091736687312762354,
// 0.17082039324993690892275210061938287063 }, {
//-0.72360679774997896964091736687312762354,
// 0.17082039324993690892275210061938287063,
//-0.72360679774997896964091736687312762354 }, {
// 0.17082039324993690892275210061938287063,
//-0.72360679774997896964091736687312762354,
//-0.72360679774997896964091736687312762354 }, {
//-0.72360679774997896964091736687312762354,
//-0.72360679774997896964091736687312762354,
//-0.72360679774997896964091736687312762354 } } , {
// 0.33333333333333333333333333333333333333,
// 0.33333333333333333333333333333333333333,
// 0.33333333333333333333333333333333333333,
// 0.33333333333333333333333333333333333333 } };
// case 3:		return { { { -0.34367339496723662642072827083693243093,
//-0.34367339496723662642072827083693243093,
//-0.9689798150982901207378151874892027072 }, {
//-0.34367339496723662642072827083693243093,
//-0.9689798150982901207378151874892027072,
//-0.34367339496723662642072827083693243093 }, {
//-0.9689798150982901207378151874892027072,
//-0.34367339496723662642072827083693243093,
//-0.34367339496723662642072827083693243093 }, {
//-0.34367339496723662642072827083693243093,
//-0.34367339496723662642072827083693243093,
//-0.34367339496723662642072827083693243093 }, {
//-0.78390550020314279176487322158837338344,
//-0.78390550020314279176487322158837338344,
// 0.35171650060942837529461966476512015033 }, {
//-0.78390550020314279176487322158837338344,
// 0.35171650060942837529461966476512015033,
//-0.78390550020314279176487322158837338344 }, {
// 0.35171650060942837529461966476512015033,
//-0.78390550020314279176487322158837338344,
//-0.78390550020314279176487322158837338344 }, {
//-0.78390550020314279176487322158837338344,
//-0.78390550020314279176487322158837338344,
//-0.78390550020314279176487322158837338344 } } , {
// 0.18162379004944980942342872025562069427,
// 0.18162379004944980942342872025562069427,
// 0.18162379004944980942342872025562069427,
// 0.18162379004944980942342872025562069427,
// 0.15170954328388352390990461307771263906,
// 0.15170954328388352390990461307771263906,
// 0.15170954328388352390990461307771263906,
// 0.15170954328388352390990461307771263906 } };
// case 4: 					case 5: return { { {
// -0.37822816147339878040530853247308433401,
//-0.37822816147339878040530853247308433401,
//-0.86531551557980365878407440258074699796 }, {
//-0.37822816147339878040530853247308433401,
//-0.86531551557980365878407440258074699796,
//-0.37822816147339878040530853247308433401 }, {
//-0.86531551557980365878407440258074699796,
//-0.37822816147339878040530853247308433401,
//-0.37822816147339878040530853247308433401 }, {
//-0.37822816147339878040530853247308433401,
//-0.37822816147339878040530853247308433401,
//-0.37822816147339878040530853247308433401 }, {
//-0.81452949937821754719535217252593878951,
//-0.81452949937821754719535217252593878951,
// 0.44358849813465264158605651757781636853 }, {
//-0.81452949937821754719535217252593878951,
// 0.44358849813465264158605651757781636853,
//-0.81452949937821754719535217252593878951 }, {
// 0.44358849813465264158605651757781636853,
//-0.81452949937821754719535217252593878951,
//-0.81452949937821754719535217252593878951 }, {
//-0.81452949937821754719535217252593878951,
//-0.81452949937821754719535217252593878951,
//-0.81452949937821754719535217252593878951 }, {
//-0.90899259174870070101623894744132112187,
//-0.091007408251299298983761052558678878131,
//-0.091007408251299298983761052558678878131 }, {
//-0.091007408251299298983761052558678878131,
//-0.90899259174870070101623894744132112187,
//-0.091007408251299298983761052558678878131 }, {
//-0.90899259174870070101623894744132112187,
//-0.90899259174870070101623894744132112187,
//-0.091007408251299298983761052558678878131 }, {
//-0.90899259174870070101623894744132112187,
//-0.091007408251299298983761052558678878131,
//-0.90899259174870070101623894744132112187 }, {
//-0.091007408251299298983761052558678878131,
//-0.90899259174870070101623894744132112187,
//-0.90899259174870070101623894744132112187 }, {
//-0.091007408251299298983761052558678878131,
//-0.091007408251299298983761052558678878131,
//-0.90899259174870070101623894744132112187 } } , {
// 0.15025056762402113439891420311104844508,
// 0.15025056762402113439891420311104844508,
// 0.15025056762402113439891420311104844508,
// 0.15025056762402113439891420311104844508,
// 0.097990724155149266058280273981770004697,
// 0.097990724155149266058280273981770004697,
// 0.097990724155149266058280273981770004697,
// 0.097990724155149266058280273981770004697,
// 0.05672802770277528858409257082700992237,
// 0.05672802770277528858409257082700992237,
// 0.05672802770277528858409257082700992237,
// 0.05672802770277528858409257082700992237,
// 0.05672802770277528858409257082700992237,
// 0.05672802770277528858409257082700992237 } };
// case 6:		return { { { -0.91865208293077729376884110208717988144,
//-0.91865208293077729376884110208717988144,
// 0.75595624879233188130652330626153964433 }, {
//-0.91865208293077729376884110208717988144,
// 0.75595624879233188130652330626153964433,
//-0.91865208293077729376884110208717988144 }, {
// 0.75595624879233188130652330626153964433,
//-0.91865208293077729376884110208717988144,
//-0.91865208293077729376884110208717988144 }, {
//-0.91865208293077729376884110208717988144,
//-0.91865208293077729376884110208717988144,
//-0.91865208293077729376884110208717988144 }, {
//-0.35532421971544897931201105847501574951,
//-0.35532421971544897931201105847501574951,
//-0.93402734085365306206396682457495275146 }, {
//-0.35532421971544897931201105847501574951,
//-0.93402734085365306206396682457495275146,
//-0.35532421971544897931201105847501574951 }, {
//-0.93402734085365306206396682457495275146,
//-0.35532421971544897931201105847501574951,
//-0.35532421971544897931201105847501574951 }, {
//-0.35532421971544897931201105847501574951,
//-0.35532421971544897931201105847501574951,
//-0.35532421971544897931201105847501574951 }, {
//-0.5707942574816959414223215612274300172,
//-0.5707942574816959414223215612274300172,
//-0.28761722755491217573303531631770994839 }, {
//-0.5707942574816959414223215612274300172,
//-0.28761722755491217573303531631770994839,
//-0.5707942574816959414223215612274300172 }, {
//-0.28761722755491217573303531631770994839,
//-0.5707942574816959414223215612274300172,
//-0.5707942574816959414223215612274300172 }, {
//-0.5707942574816959414223215612274300172,
//-0.5707942574816959414223215612274300172,
//-0.5707942574816959414223215612274300172 }, {
// 0.20601132958329828273486227812187937257,
//-0.87267799624996494940152894478854603924,
//-0.46065533708336838393180438854478729409 }, {
// 0.20601132958329828273486227812187937257,
//-0.87267799624996494940152894478854603924,
//-0.87267799624996494940152894478854603924 }, {
//-0.87267799624996494940152894478854603924,
//-0.87267799624996494940152894478854603924,
// 0.20601132958329828273486227812187937257 }, {
//-0.46065533708336838393180438854478729409,
// 0.20601132958329828273486227812187937257,
//-0.87267799624996494940152894478854603924 }, {
//-0.87267799624996494940152894478854603924,
//-0.46065533708336838393180438854478729409,
// 0.20601132958329828273486227812187937257 }, {
//-0.87267799624996494940152894478854603924,
// 0.20601132958329828273486227812187937257,
//-0.87267799624996494940152894478854603924 }, {
//-0.46065533708336838393180438854478729409,
//-0.87267799624996494940152894478854603924,
// 0.20601132958329828273486227812187937257 }, {
//-0.87267799624996494940152894478854603924,
//-0.46065533708336838393180438854478729409,
//-0.87267799624996494940152894478854603924 }, {
//-0.87267799624996494940152894478854603924,
//-0.87267799624996494940152894478854603924,
//-0.46065533708336838393180438854478729409 }, {
//-0.87267799624996494940152894478854603924,
// 0.20601132958329828273486227812187937257,
//-0.46065533708336838393180438854478729409 }, {
//-0.46065533708336838393180438854478729409,
//-0.87267799624996494940152894478854603924,
//-0.87267799624996494940152894478854603924 }, {
// 0.20601132958329828273486227812187937257,
//-0.46065533708336838393180438854478729409,
//-0.87267799624996494940152894478854603924 } } , {
// 0.013436281407094190597350983261249151023,
// 0.013436281407094190597350983261249151023,
// 0.013436281407094190597350983261249151023,
// 0.013436281407094190597350983261249151023,
// 0.073809575391539629460204370471634688549,
// 0.073809575391539629460204370471634688549,
// 0.073809575391539629460204370471634688549,
// 0.073809575391539629460204370471634688549,
// 0.053230333677556656132920836743306636618,
// 0.053230333677556656132920836743306636618,
// 0.053230333677556656132920836743306636618,
// 0.053230333677556656132920836743306636618,
// 0.064285714285714285714285714285714285714,
// 0.064285714285714285714285714285714285714,
// 0.064285714285714285714285714285714285714,
// 0.064285714285714285714285714285714285714,
// 0.064285714285714285714285714285714285714,
// 0.064285714285714285714285714285714285714,
// 0.064285714285714285714285714285714285714,
// 0.064285714285714285714285714285714285714,
// 0.064285714285714285714285714285714285714,
// 0.064285714285714285714285714285714285714,
// 0.064285714285714285714285714285714285714,
// 0.064285714285714285714285714285714285714 } };
//					case 7:		return { { { -0.5, -0.5,
//-0.5
//}, { -0.36859770044359440115314000081337702315,
//-0.36859770044359440115314000081337702315,
//-0.89420689866921679654057999755986893056 }, {
//-0.36859770044359440115314000081337702315,
//-0.89420689866921679654057999755986893056,
//-0.36859770044359440115314000081337702315 }, {
//-0.89420689866921679654057999755986893056,
//-0.36859770044359440115314000081337702315,
//-0.36859770044359440115314000081337702315 }, {
//-0.36859770044359440115314000081337702315,
//-0.36859770044359440115314000081337702315,
//-0.36859770044359440115314000081337702315 }, {
//-0.89902035480320726247389235402687506805,
//-0.10097964519679273752610764597312493195,
//-0.10097964519679273752610764597312493195 }, {
//-0.10097964519679273752610764597312493195,
//-0.89902035480320726247389235402687506805,
//-0.10097964519679273752610764597312493195 }, {
//-0.89902035480320726247389235402687506805,
//-0.89902035480320726247389235402687506805,
//-0.10097964519679273752610764597312493195 }, {
//-0.89902035480320726247389235402687506805,
//-0.10097964519679273752610764597312493195,
//-0.89902035480320726247389235402687506805 }, {
//-0.10097964519679273752610764597312493195,
//-0.89902035480320726247389235402687506805,
//-0.89902035480320726247389235402687506805 }, {
//-0.10097964519679273752610764597312493195,
//-0.10097964519679273752610764597312493195,
//-0.89902035480320726247389235402687506805 }, {
// 0.15034327517400004696648315404461503933,
//-0.62233233794799790452713779229082848999,
//-0.90567859927800423791220756946295805935 }, {
// 0.15034327517400004696648315404461503933,
//-0.62233233794799790452713779229082848999,
//-0.62233233794799790452713779229082848999 }, {
//-0.62233233794799790452713779229082848999,
//-0.62233233794799790452713779229082848999,
// 0.15034327517400004696648315404461503933 }, {
//-0.90567859927800423791220756946295805935,
// 0.15034327517400004696648315404461503933,
//-0.62233233794799790452713779229082848999 }, {
//-0.62233233794799790452713779229082848999,
//-0.90567859927800423791220756946295805935,
// 0.15034327517400004696648315404461503933 }, {
//-0.62233233794799790452713779229082848999,
// 0.15034327517400004696648315404461503933,
//-0.62233233794799790452713779229082848999 }, {
//-0.90567859927800423791220756946295805935,
//-0.62233233794799790452713779229082848999,
// 0.15034327517400004696648315404461503933 }, {
//-0.62233233794799790452713779229082848999,
//-0.90567859927800423791220756946295805935,
//-0.62233233794799790452713779229082848999 }, {
//-0.62233233794799790452713779229082848999,
//-0.62233233794799790452713779229082848999,
//-0.90567859927800423791220756946295805935 }, {
//-0.62233233794799790452713779229082848999,
// 0.15034327517400004696648315404461503933,
//-0.90567859927800423791220756946295805935 }, {
//-0.90567859927800423791220756946295805935,
//-0.62233233794799790452713779229082848999,
//-0.62233233794799790452713779229082848999 }, {
// 0.15034327517400004696648315404461503933,
//-0.90567859927800423791220756946295805935,
//-0.62233233794799790452713779229082848999 }, {
// 0.62166048219709712223621075969646478759,
//-0.95746905491703350802232779700036011785,
//-0.70672237236303010619155516569574455189 }, {
// 0.62166048219709712223621075969646478759,
//-0.95746905491703350802232779700036011785,
//-0.95746905491703350802232779700036011785 }, {
//-0.95746905491703350802232779700036011785,
//-0.95746905491703350802232779700036011785,
// 0.62166048219709712223621075969646478759 }, {
//-0.70672237236303010619155516569574455189,
// 0.62166048219709712223621075969646478759,
//-0.95746905491703350802232779700036011785 }, {
//-0.95746905491703350802232779700036011785,
//-0.70672237236303010619155516569574455189,
// 0.62166048219709712223621075969646478759 }, {
//-0.95746905491703350802232779700036011785,
// 0.62166048219709712223621075969646478759,
//-0.95746905491703350802232779700036011785 }, {
//-0.70672237236303010619155516569574455189,
//-0.95746905491703350802232779700036011785,
// 0.62166048219709712223621075969646478759 }, {
//-0.95746905491703350802232779700036011785,
//-0.70672237236303010619155516569574455189,
//-0.95746905491703350802232779700036011785 }, {
//-0.95746905491703350802232779700036011785,
//-0.95746905491703350802232779700036011785,
//-0.70672237236303010619155516569574455189 }, {
//-0.95746905491703350802232779700036011785,
// 0.62166048219709712223621075969646478759,
//-0.70672237236303010619155516569574455189 }, {
//-0.70672237236303010619155516569574455189,
//-0.95746905491703350802232779700036011785,
//-0.95746905491703350802232779700036011785 }, {
// 0.62166048219709712223621075969646478759,
//-0.70672237236303010619155516569574455189,
//-0.95746905491703350802232779700036011785 } } , {
// 0.12731371928550779848077124815630184245,
// 0.05643944161328937210171489439806231665,
// 0.05643944161328937210171489439806231665,
// 0.05643944161328937210171489439806231665,
// 0.05643944161328937210171489439806231665,
// 0.04252923711047677324569976544392327613,
// 0.04252923711047677324569976544392327613,
// 0.04252923711047677324569976544392327613,
// 0.04252923711047677324569976544392327613,
// 0.04252923711047677324569976544392327613,
// 0.04252923711047677324569976544392327613,
// 0.049609507637779495159487414921974827082,
// 0.049609507637779495159487414921974827082,
// 0.049609507637779495159487414921974827082,
// 0.049609507637779495159487414921974827082,
// 0.049609507637779495159487414921974827082,
// 0.049609507637779495159487414921974827082,
// 0.049609507637779495159487414921974827082,
// 0.049609507637779495159487414921974827082,
// 0.049609507637779495159487414921974827082,
// 0.049609507637779495159487414921974827082,
// 0.049609507637779495159487414921974827082,
// 0.049609507637779495159487414921974827082,
// 0.010814361106537788754804577988128720209,
// 0.010814361106537788754804577988128720209,
// 0.010814361106537788754804577988128720209,
// 0.010814361106537788754804577988128720209,
// 0.010814361106537788754804577988128720209,
// 0.010814361106537788754804577988128720209,
// 0.010814361106537788754804577988128720209,
// 0.010814361106537788754804577988128720209,
// 0.010814361106537788754804577988128720209,
// 0.010814361106537788754804577988128720209,
// 0.010814361106537788754804577988128720209,
// 0.010814361106537788754804577988128720209 } };
//					case 8:		return { { {
//-0.78409455007557830303561343854585276851,
//-0.78409455007557830303561343854585276851,
// 0.35228365022673490910684031563755830553 }, {
//-0.78409455007557830303561343854585276851,
// 0.35228365022673490910684031563755830553,
//-0.78409455007557830303561343854585276851 }, {
// 0.35228365022673490910684031563755830553,
//-0.78409455007557830303561343854585276851,
//-0.78409455007557830303561343854585276851 }, {
//-0.78409455007557830303561343854585276851,
//-0.78409455007557830303561343854585276851,
//-0.78409455007557830303561343854585276851 }, {
//-0.629781024434826859446183999496503797,
//-0.629781024434826859446183999496503797,
//-0.110656926695519421661448001510488609 }, {
//-0.629781024434826859446183999496503797,
//-0.110656926695519421661448001510488609,
//-0.629781024434826859446183999496503797 }, {
//-0.110656926695519421661448001510488609,
//-0.629781024434826859446183999496503797,
//-0.629781024434826859446183999496503797 }, {
//-0.629781024434826859446183999496503797,
//-0.629781024434826859446183999496503797,
//-0.629781024434826859446183999496503797 }, {
//-0.91536691263046543675183591431855057189,
//-0.91536691263046543675183591431855057189,
// 0.74610073789139631025550774295565171568 }, {
//-0.91536691263046543675183591431855057189,
// 0.74610073789139631025550774295565171568,
//-0.91536691263046543675183591431855057189 }, {
// 0.74610073789139631025550774295565171568,
//-0.91536691263046543675183591431855057189,
//-0.91536691263046543675183591431855057189 }, {
//-0.91536691263046543675183591431855057189,
//-0.91536691263046543675183591431855057189,
//-0.91536691263046543675183591431855057189 }, {
//-0.37163658175192200927334245312373662702,
//-0.37163658175192200927334245312373662702,
//-0.88509025474423397217997264062879011895 }, {
//-0.37163658175192200927334245312373662702,
//-0.88509025474423397217997264062879011895,
//-0.37163658175192200927334245312373662702 }, {
//-0.88509025474423397217997264062879011895,
//-0.37163658175192200927334245312373662702,
//-0.37163658175192200927334245312373662702 }, {
//-0.37163658175192200927334245312373662702,
//-0.37163658175192200927334245312373662702,
//-0.37163658175192200927334245312373662702 }, {
//-0.12881734283233958954373367640680593526,
//-0.87118265716766041045626632359319406474,
//-0.87118265716766041045626632359319406474 }, {
//-0.87118265716766041045626632359319406474,
//-0.12881734283233958954373367640680593526,
//-0.87118265716766041045626632359319406474 }, {
//-0.12881734283233958954373367640680593526,
//-0.12881734283233958954373367640680593526,
//-0.87118265716766041045626632359319406474 }, {
//-0.12881734283233958954373367640680593526,
//-0.87118265716766041045626632359319406474,
//-0.12881734283233958954373367640680593526 }, {
//-0.87118265716766041045626632359319406474,
//-0.12881734283233958954373367640680593526,
//-0.12881734283233958954373367640680593526 }, {
//-0.87118265716766041045626632359319406474,
//-0.87118265716766041045626632359319406474,
//-0.12881734283233958954373367640680593526 }, {
// 0.43492812685261664658005711425116290088,
//-0.95713213974573885031100182115353118963,
//-0.52066384736113894595805347194410052161 }, {
// 0.43492812685261664658005711425116290088,
//-0.95713213974573885031100182115353118963,
//-0.95713213974573885031100182115353118963 }, {
//-0.95713213974573885031100182115353118963,
//-0.95713213974573885031100182115353118963,
// 0.43492812685261664658005711425116290088 }, {
//-0.52066384736113894595805347194410052161,
// 0.43492812685261664658005711425116290088,
//-0.95713213974573885031100182115353118963 }, {
//-0.95713213974573885031100182115353118963,
//-0.52066384736113894595805347194410052161,
// 0.43492812685261664658005711425116290088 }, {
//-0.95713213974573885031100182115353118963,
// 0.43492812685261664658005711425116290088,
//-0.95713213974573885031100182115353118963 }, {
//-0.52066384736113894595805347194410052161,
//-0.95713213974573885031100182115353118963,
// 0.43492812685261664658005711425116290088 }, {
//-0.95713213974573885031100182115353118963,
//-0.52066384736113894595805347194410052161,
//-0.95713213974573885031100182115353118963 }, {
//-0.95713213974573885031100182115353118963,
//-0.95713213974573885031100182115353118963,
//-0.52066384736113894595805347194410052161 }, {
//-0.95713213974573885031100182115353118963,
// 0.43492812685261664658005711425116290088,
//-0.52066384736113894595805347194410052161 }, {
//-0.52066384736113894595805347194410052161,
//-0.95713213974573885031100182115353118963,
//-0.95713213974573885031100182115353118963 }, {
// 0.43492812685261664658005711425116290088,
//-0.52066384736113894595805347194410052161,
//-0.95713213974573885031100182115353118963 }, {
// 0.16759475660428881185333950841466005283,
//-0.59172133224794175913941940243963088242,
//-0.98415209210840529357450070353539828799 }, {
// 0.16759475660428881185333950841466005283,
//-0.59172133224794175913941940243963088242,
//-0.59172133224794175913941940243963088242 }, {
//-0.59172133224794175913941940243963088242,
//-0.59172133224794175913941940243963088242,
// 0.16759475660428881185333950841466005283 }, {
//-0.98415209210840529357450070353539828799,
// 0.16759475660428881185333950841466005283,
//-0.59172133224794175913941940243963088242 }, {
//-0.59172133224794175913941940243963088242,
//-0.98415209210840529357450070353539828799,
// 0.16759475660428881185333950841466005283 }, {
//-0.59172133224794175913941940243963088242,
// 0.16759475660428881185333950841466005283,
//-0.59172133224794175913941940243963088242 }, {
//-0.98415209210840529357450070353539828799,
//-0.59172133224794175913941940243963088242,
// 0.16759475660428881185333950841466005283 }, {
//-0.59172133224794175913941940243963088242,
//-0.98415209210840529357450070353539828799,
//-0.59172133224794175913941940243963088242 }, {
//-0.59172133224794175913941940243963088242,
//-0.59172133224794175913941940243963088242,
//-0.98415209210840529357450070353539828799 }, {
//-0.59172133224794175913941940243963088242,
// 0.16759475660428881185333950841466005283,
//-0.98415209210840529357450070353539828799 }, {
//-0.98415209210840529357450070353539828799,
//-0.59172133224794175913941940243963088242,
//-0.59172133224794175913941940243963088242 }, {
// 0.16759475660428881185333950841466005283,
//-0.98415209210840529357450070353539828799,
//-0.59172133224794175913941940243963088242 } } , {
// 0.035235534544545108539097313565645827902,
// 0.035235534544545108539097313565645827902,
// 0.035235534544545108539097313565645827902,
// 0.035235534544545108539097313565645827902,
// 0.069375663418318043619861410028472677606,
// 0.069375663418318043619861410028472677606,
// 0.069375663418318043619861410028472677606,
// 0.069375663418318043619861410028472677606,
// 0.010033674871386933040967405110292766211,
// 0.010033674871386933040967405110292766211,
// 0.010033674871386933040967405110292766211,
// 0.010033674871386933040967405110292766211,
// 0.055685043809246530275971813171748834223,
// 0.055685043809246530275971813171748834223,
// 0.055685043809246530275971813171748834223,
// 0.055685043809246530275971813171748834223,
// 0.048374573681745097334959362864457783412,
// 0.048374573681745097334959362864457783412,
// 0.048374573681745097334959362864457783412,
// 0.048374573681745097334959362864457783412,
// 0.048374573681745097334959362864457783412,
// 0.048374573681745097334959362864457783412,
// 0.0095425371877925775336250150191510735196,
// 0.0095425371877925775336250150191510735196,
// 0.0095425371877925775336250150191510735196,
// 0.0095425371877925775336250150191510735196,
// 0.0095425371877925775336250150191510735196,
// 0.0095425371877925775336250150191510735196,
// 0.0095425371877925775336250150191510735196,
// 0.0095425371877925775336250150191510735196,
// 0.0095425371877925775336250150191510735196,
// 0.0095425371877925775336250150191510735196,
// 0.0095425371877925775336250150191510735196,
// 0.0095425371877925775336250150191510735196,
// 0.020604648201280446418040434034344443905,
// 0.020604648201280446418040434034344443905,
// 0.020604648201280446418040434034344443905,
// 0.020604648201280446418040434034344443905,
// 0.020604648201280446418040434034344443905,
// 0.020604648201280446418040434034344443905,
// 0.020604648201280446418040434034344443905,
// 0.020604648201280446418040434034344443905,
// 0.020604648201280446418040434034344443905,
// 0.020604648201280446418040434034344443905,
// 0.020604648201280446418040434034344443905,
// 0.020604648201280446418040434034344443905 } };
//					case 9:		return { { { -0.5, -0.5,
//-0.5
//}, { -0.99999999876036601109069831529161985423,
//-0.99999999876036601109069831529161985423,
// 0.9999999962810980332720949458748595627 }, {
//-0.99999999876036601109069831529161985423,
// 0.9999999962810980332720949458748595627,
//-0.99999999876036601109069831529161985423 }, {
// 0.9999999962810980332720949458748595627,
//-0.99999999876036601109069831529161985423,
//-0.99999999876036601109069831529161985423 }, {
//-0.99999999876036601109069831529161985423,
//-0.99999999876036601109069831529161985423,
//-0.99999999876036601109069831529161985423 }, {
//-0.67845092920947681169601712222921603853,
//-0.67845092920947681169601712222921603853,
// 0.035352787628430435088051366687648115583 }, {
//-0.67845092920947681169601712222921603853,
// 0.035352787628430435088051366687648115583,
//-0.67845092920947681169601712222921603853 }, {
// 0.035352787628430435088051366687648115583,
//-0.67845092920947681169601712222921603853,
//-0.67845092920947681169601712222921603853 }, {
//-0.67845092920947681169601712222921603853,
//-0.67845092920947681169601712222921603853,
//-0.67845092920947681169601712222921603853 }, {
//-0.35544695635715805159186380319222134027,
//-0.35544695635715805159186380319222134027,
//-0.93365913092852584522440859042333597919 }, {
//-0.35544695635715805159186380319222134027,
//-0.93365913092852584522440859042333597919,
//-0.35544695635715805159186380319222134027 }, {
//-0.93365913092852584522440859042333597919,
//-0.35544695635715805159186380319222134027,
//-0.35544695635715805159186380319222134027 }, {
//-0.35544695635715805159186380319222134027,
//-0.35544695635715805159186380319222134027,
//-0.35544695635715805159186380319222134027 }, {
//-0.90978216330917283407807385501200614331,
//-0.90978216330917283407807385501200614331,
// 0.72934648992751850223422156503601842994 }, {
//-0.90978216330917283407807385501200614331,
// 0.72934648992751850223422156503601842994,
//-0.90978216330917283407807385501200614331 }, {
// 0.72934648992751850223422156503601842994,
//-0.90978216330917283407807385501200614331,
//-0.90978216330917283407807385501200614331 }, {
//-0.90978216330917283407807385501200614331,
//-0.90978216330917283407807385501200614331,
//-0.90978216330917283407807385501200614331 }, {
//-0.77540690799124790934402070093968875729,
//-0.22459309200875209065597929906031124271,
//-0.22459309200875209065597929906031124271 }, {
//-0.22459309200875209065597929906031124271,
//-0.77540690799124790934402070093968875729,
//-0.22459309200875209065597929906031124271 }, {
//-0.77540690799124790934402070093968875729,
//-0.77540690799124790934402070093968875729,
//-0.22459309200875209065597929906031124271 }, {
//-0.77540690799124790934402070093968875729,
//-0.22459309200875209065597929906031124271,
//-0.77540690799124790934402070093968875729 }, {
//-0.22459309200875209065597929906031124271,
//-0.77540690799124790934402070093968875729,
//-0.77540690799124790934402070093968875729 }, {
//-0.22459309200875209065597929906031124271,
//-0.22459309200875209065597929906031124271,
//-0.77540690799124790934402070093968875729 }, {
//-0.994890841533917338064711949103255315,
//-0.082257102495081453456435781598997046246,
//-0.84059495347591975502241648769875059251 }, {
//-0.994890841533917338064711949103255315,
//-0.082257102495081453456435781598997046246,
//-0.082257102495081453456435781598997046246 }, {
//-0.082257102495081453456435781598997046246,
//-0.082257102495081453456435781598997046246,
//-0.994890841533917338064711949103255315 }, {
//-0.84059495347591975502241648769875059251,
//-0.994890841533917338064711949103255315,
//-0.082257102495081453456435781598997046246 }, {
//-0.082257102495081453456435781598997046246,
//-0.84059495347591975502241648769875059251,
//-0.994890841533917338064711949103255315 }, {
//-0.082257102495081453456435781598997046246,
//-0.994890841533917338064711949103255315,
//-0.082257102495081453456435781598997046246 }, {
//-0.84059495347591975502241648769875059251,
//-0.082257102495081453456435781598997046246,
//-0.994890841533917338064711949103255315 }, {
//-0.082257102495081453456435781598997046246,
//-0.84059495347591975502241648769875059251,
//-0.082257102495081453456435781598997046246 }, {
//-0.082257102495081453456435781598997046246,
//-0.082257102495081453456435781598997046246,
//-0.84059495347591975502241648769875059251 }, {
//-0.082257102495081453456435781598997046246,
//-0.994890841533917338064711949103255315,
//-0.84059495347591975502241648769875059251 }, {
//-0.84059495347591975502241648769875059251,
//-0.082257102495081453456435781598997046246,
//-0.082257102495081453456435781598997046246 }, {
//-0.994890841533917338064711949103255315,
//-0.84059495347591975502241648769875059251,
//-0.082257102495081453456435781598997046246 }, {
// 0.43670065288414901811354448912928750012,
//-0.93244825862932284418902737202159413719,
//-0.57180413562550332973548974508609922574 }, {
// 0.43670065288414901811354448912928750012,
//-0.93244825862932284418902737202159413719,
//-0.93244825862932284418902737202159413719 }, {
//-0.93244825862932284418902737202159413719,
//-0.93244825862932284418902737202159413719,
// 0.43670065288414901811354448912928750012 }, {
//-0.57180413562550332973548974508609922574,
// 0.43670065288414901811354448912928750012,
//-0.93244825862932284418902737202159413719 }, {
//-0.93244825862932284418902737202159413719,
//-0.57180413562550332973548974508609922574,
// 0.43670065288414901811354448912928750012 }, {
//-0.93244825862932284418902737202159413719,
// 0.43670065288414901811354448912928750012,
//-0.93244825862932284418902737202159413719 }, {
//-0.57180413562550332973548974508609922574,
//-0.93244825862932284418902737202159413719,
// 0.43670065288414901811354448912928750012 }, {
//-0.93244825862932284418902737202159413719,
//-0.57180413562550332973548974508609922574,
//-0.93244825862932284418902737202159413719 }, {
//-0.93244825862932284418902737202159413719,
//-0.93244825862932284418902737202159413719,
//-0.57180413562550332973548974508609922574 }, {
//-0.93244825862932284418902737202159413719,
// 0.43670065288414901811354448912928750012,
//-0.57180413562550332973548974508609922574 }, {
//-0.57180413562550332973548974508609922574,
//-0.93244825862932284418902737202159413719,
//-0.93244825862932284418902737202159413719 }, {
// 0.43670065288414901811354448912928750012,
//-0.57180413562550332973548974508609922574,
//-0.93244825862932284418902737202159413719 }, {
//-0.93116817884364945982158957577713628669,
//-0.63271726038014422042206258028726154683,
// 0.19660269960393790066571473635165938035 }, {
//-0.93116817884364945982158957577713628669,
//-0.63271726038014422042206258028726154683,
//-0.63271726038014422042206258028726154683 }, {
//-0.63271726038014422042206258028726154683,
//-0.63271726038014422042206258028726154683,
//-0.93116817884364945982158957577713628669 }, {
// 0.19660269960393790066571473635165938035,
//-0.93116817884364945982158957577713628669,
//-0.63271726038014422042206258028726154683 }, {
//-0.63271726038014422042206258028726154683,
// 0.19660269960393790066571473635165938035,
//-0.93116817884364945982158957577713628669 }, {
//-0.63271726038014422042206258028726154683,
//-0.93116817884364945982158957577713628669,
//-0.63271726038014422042206258028726154683 }, {
// 0.19660269960393790066571473635165938035,
//-0.63271726038014422042206258028726154683,
//-0.93116817884364945982158957577713628669 }, {
//-0.63271726038014422042206258028726154683,
// 0.19660269960393790066571473635165938035,
//-0.63271726038014422042206258028726154683 }, {
//-0.63271726038014422042206258028726154683,
//-0.63271726038014422042206258028726154683,
// 0.19660269960393790066571473635165938035 }, {
//-0.63271726038014422042206258028726154683,
//-0.93116817884364945982158957577713628669,
// 0.19660269960393790066571473635165938035 }, {
// 0.19660269960393790066571473635165938035,
//-0.63271726038014422042206258028726154683,
//-0.63271726038014422042206258028726154683 }, {
//-0.93116817884364945982158957577713628669,
// 0.19660269960393790066571473635165938035,
//-0.63271726038014422042206258028726154683 } } , {
// 0.07734739854997367644741906257553729918, 8.5759042345675186085903030258742437834e-05,
// 8.5759042345675186085903030258742437834e-05, 8.5759042345675186085903030258742437834e-05,
// 8.5759042345675186085903030258742437834e-05,
// 0.030897784616567277310532826869162311937,
// 0.030897784616567277310532826869162311937,
// 0.030897784616567277310532826869162311937,
// 0.030897784616567277310532826869162311937,
// 0.039417216447239047523644877985518239107,
// 0.039417216447239047523644877985518239107,
// 0.039417216447239047523644877985518239107,
// 0.039417216447239047523644877985518239107,
// 0.010751973306154910354337519688619045729,
// 0.010751973306154910354337519688619045729,
// 0.010751973306154910354337519688619045729,
// 0.010751973306154910354337519688619045729,
// 0.05084544013826995406889748131136091475,
// 0.05084544013826995406889748131136091475,
// 0.05084544013826995406889748131136091475,
// 0.05084544013826995406889748131136091475,
// 0.05084544013826995406889748131136091475,
// 0.05084544013826995406889748131136091475,
// 0.011179229597731402927583520512290878612,
// 0.011179229597731402927583520512290878612,
// 0.011179229597731402927583520512290878612,
// 0.011179229597731402927583520512290878612,
// 0.011179229597731402927583520512290878612,
// 0.011179229597731402927583520512290878612,
// 0.011179229597731402927583520512290878612,
// 0.011179229597731402927583520512290878612,
// 0.011179229597731402927583520512290878612,
// 0.011179229597731402927583520512290878612,
// 0.011179229597731402927583520512290878612,
// 0.011179229597731402927583520512290878612,
// 0.013646079136993770600501763121325612648,
// 0.013646079136993770600501763121325612648,
// 0.013646079136993770600501763121325612648,
// 0.013646079136993770600501763121325612648,
// 0.013646079136993770600501763121325612648,
// 0.013646079136993770600501763121325612648,
// 0.013646079136993770600501763121325612648,
// 0.013646079136993770600501763121325612648,
// 0.013646079136993770600501763121325612648,
// 0.013646079136993770600501763121325612648,
// 0.013646079136993770600501763121325612648,
// 0.013646079136993770600501763121325612648,
// 0.027366554623984184053091789082666607808,
// 0.027366554623984184053091789082666607808,
// 0.027366554623984184053091789082666607808,
// 0.027366554623984184053091789082666607808,
// 0.027366554623984184053091789082666607808,
// 0.027366554623984184053091789082666607808,
// 0.027366554623984184053091789082666607808,
// 0.027366554623984184053091789082666607808,
// 0.027366554623984184053091789082666607808,
// 0.027366554623984184053091789082666607808,
// 0.027366554623984184053091789082666607808,
// 0.027366554623984184053091789082666607808 } };
// case 10:	return { { { -0.5, -0.5, -0.5 }, {
//-0.3754998626096227045403833626263450896,
//-0.3754998626096227045403833626263450896,
//-0.87350041217113188637884991212096473119 }, {
//-0.3754998626096227045403833626263450896,
//-0.87350041217113188637884991212096473119,
//-0.3754998626096227045403833626263450896 }, {
//-0.87350041217113188637884991212096473119,
//-0.3754998626096227045403833626263450896,
//-0.3754998626096227045403833626263450896 }, {
//-0.3754998626096227045403833626263450896,
//-0.3754998626096227045403833626263450896,
//-0.3754998626096227045403833626263450896 }, {
//-0.77138069228530769882525760469269910942,
//-0.77138069228530769882525760469269910942,
// 0.31414207685592309647577281407809732827 }, {
//-0.77138069228530769882525760469269910942,
// 0.31414207685592309647577281407809732827,
//-0.77138069228530769882525760469269910942 }, {
// 0.31414207685592309647577281407809732827,
//-0.77138069228530769882525760469269910942,
//-0.77138069228530769882525760469269910942 }, {
//-0.77138069228530769882525760469269910942,
//-0.77138069228530769882525760469269910942,
//-0.77138069228530769882525760469269910942 }, {
//-0.66902794876077789679101975111094717095,
//-0.17913852156206901142420431149697662502,
//-0.97269500811508408036057162589509957901 }, {
//-0.66902794876077789679101975111094717095,
//-0.17913852156206901142420431149697662502,
//-0.17913852156206901142420431149697662502 }, {
//-0.17913852156206901142420431149697662502,
//-0.17913852156206901142420431149697662502,
//-0.66902794876077789679101975111094717095 }, {
//-0.97269500811508408036057162589509957901,
//-0.66902794876077789679101975111094717095,
//-0.17913852156206901142420431149697662502 }, {
//-0.17913852156206901142420431149697662502,
//-0.97269500811508408036057162589509957901,
//-0.66902794876077789679101975111094717095 }, {
//-0.17913852156206901142420431149697662502,
//-0.66902794876077789679101975111094717095,
//-0.17913852156206901142420431149697662502 }, {
//-0.97269500811508408036057162589509957901,
//-0.17913852156206901142420431149697662502,
//-0.66902794876077789679101975111094717095 }, {
//-0.17913852156206901142420431149697662502,
//-0.97269500811508408036057162589509957901,
//-0.17913852156206901142420431149697662502 }, {
//-0.17913852156206901142420431149697662502,
//-0.17913852156206901142420431149697662502,
//-0.97269500811508408036057162589509957901 }, {
//-0.17913852156206901142420431149697662502,
//-0.66902794876077789679101975111094717095,
//-0.97269500811508408036057162589509957901 }, {
//-0.97269500811508408036057162589509957901,
//-0.17913852156206901142420431149697662502,
//-0.17913852156206901142420431149697662502 }, {
//-0.66902794876077789679101975111094717095,
//-0.97269500811508408036057162589509957901,
//-0.17913852156206901142420431149697662502 }, {
// 0.8859775346904097323952611738365015259,
//-0.98772398235041850430481257350316929785,
//-0.9105295699895727237856360268301629302 }, {
// 0.8859775346904097323952611738365015259,
//-0.98772398235041850430481257350316929785,
//-0.98772398235041850430481257350316929785 }, {
//-0.98772398235041850430481257350316929785,
//-0.98772398235041850430481257350316929785,
// 0.8859775346904097323952611738365015259 }, {
//-0.9105295699895727237856360268301629302,
// 0.8859775346904097323952611738365015259,
//-0.98772398235041850430481257350316929785 }, {
//-0.98772398235041850430481257350316929785,
//-0.9105295699895727237856360268301629302,
// 0.8859775346904097323952611738365015259 }, {
//-0.98772398235041850430481257350316929785,
// 0.8859775346904097323952611738365015259,
//-0.98772398235041850430481257350316929785 }, {
//-0.9105295699895727237856360268301629302,
//-0.98772398235041850430481257350316929785,
// 0.8859775346904097323952611738365015259 }, {
//-0.98772398235041850430481257350316929785,
//-0.9105295699895727237856360268301629302,
//-0.98772398235041850430481257350316929785 }, {
//-0.98772398235041850430481257350316929785,
//-0.98772398235041850430481257350316929785,
//-0.9105295699895727237856360268301629302 }, {
//-0.98772398235041850430481257350316929785,
// 0.8859775346904097323952611738365015259,
//-0.9105295699895727237856360268301629302 }, {
//-0.9105295699895727237856360268301629302,
//-0.98772398235041850430481257350316929785,
//-0.98772398235041850430481257350316929785 }, {
// 0.8859775346904097323952611738365015259,
//-0.9105295699895727237856360268301629302,
//-0.98772398235041850430481257350316929785 }, {
//-0.045619240191439298911787183406185558751,
//-0.75789963770882114801220999680989894474,
//-0.43858148439091840506379282297401655176 }, {
//-0.045619240191439298911787183406185558751,
//-0.75789963770882114801220999680989894474,
//-0.75789963770882114801220999680989894474 }, {
//-0.75789963770882114801220999680989894474,
//-0.75789963770882114801220999680989894474,
//-0.045619240191439298911787183406185558751 }, {
//-0.43858148439091840506379282297401655176,
//-0.045619240191439298911787183406185558751,
//-0.75789963770882114801220999680989894474 }, {
//-0.75789963770882114801220999680989894474,
//-0.43858148439091840506379282297401655176,
//-0.045619240191439298911787183406185558751 }, {
//-0.75789963770882114801220999680989894474,
//-0.045619240191439298911787183406185558751,
//-0.75789963770882114801220999680989894474 }, {
//-0.43858148439091840506379282297401655176,
//-0.75789963770882114801220999680989894474,
//-0.045619240191439298911787183406185558751 }, {
//-0.75789963770882114801220999680989894474,
//-0.43858148439091840506379282297401655176,
//-0.75789963770882114801220999680989894474 }, {
//-0.75789963770882114801220999680989894474,
//-0.75789963770882114801220999680989894474,
//-0.43858148439091840506379282297401655176 }, {
//-0.75789963770882114801220999680989894474,
//-0.045619240191439298911787183406185558751,
//-0.43858148439091840506379282297401655176 }, {
//-0.43858148439091840506379282297401655176,
//-0.75789963770882114801220999680989894474,
//-0.75789963770882114801220999680989894474 }, {
//-0.045619240191439298911787183406185558751,
//-0.43858148439091840506379282297401655176,
//-0.75789963770882114801220999680989894474 }, {
// 0.18851253896001405132314006877136059652,
//-0.93444106356711465845055795933535161918,
//-0.31963041182578473442202415010065735815 }, {
// 0.18851253896001405132314006877136059652,
//-0.93444106356711465845055795933535161918,
//-0.93444106356711465845055795933535161918 }, {
//-0.93444106356711465845055795933535161918,
//-0.93444106356711465845055795933535161918,
// 0.18851253896001405132314006877136059652 }, {
//-0.31963041182578473442202415010065735815,
// 0.18851253896001405132314006877136059652,
//-0.93444106356711465845055795933535161918 }, {
//-0.93444106356711465845055795933535161918,
//-0.31963041182578473442202415010065735815,
// 0.18851253896001405132314006877136059652 }, {
//-0.93444106356711465845055795933535161918,
// 0.18851253896001405132314006877136059652,
//-0.93444106356711465845055795933535161918 }, {
//-0.31963041182578473442202415010065735815,
//-0.93444106356711465845055795933535161918,
// 0.18851253896001405132314006877136059652 }, {
//-0.93444106356711465845055795933535161918,
//-0.31963041182578473442202415010065735815,
//-0.93444106356711465845055795933535161918 }, {
//-0.93444106356711465845055795933535161918,
//-0.93444106356711465845055795933535161918,
//-0.31963041182578473442202415010065735815 }, {
//-0.93444106356711465845055795933535161918,
// 0.18851253896001405132314006877136059652,
//-0.31963041182578473442202415010065735815 }, {
//-0.31963041182578473442202415010065735815,
//-0.93444106356711465845055795933535161918,
//-0.93444106356711465845055795933535161918 }, {
// 0.18851253896001405132314006877136059652,
//-0.31963041182578473442202415010065735815,
//-0.93444106356711465845055795933535161918 }, {
// 0.60235456931668878246228336815716196359,
//-0.93502943687035390432897012004314760659,
//-0.73229569557598097380434312807086675041 }, {
// 0.60235456931668878246228336815716196359,
//-0.93502943687035390432897012004314760659,
//-0.93502943687035390432897012004314760659 }, {
//-0.93502943687035390432897012004314760659,
//-0.93502943687035390432897012004314760659,
// 0.60235456931668878246228336815716196359 }, {
//-0.73229569557598097380434312807086675041,
// 0.60235456931668878246228336815716196359,
//-0.93502943687035390432897012004314760659 }, {
//-0.93502943687035390432897012004314760659,
//-0.73229569557598097380434312807086675041,
// 0.60235456931668878246228336815716196359 }, {
//-0.93502943687035390432897012004314760659,
// 0.60235456931668878246228336815716196359,
//-0.93502943687035390432897012004314760659 }, {
//-0.73229569557598097380434312807086675041,
//-0.93502943687035390432897012004314760659,
// 0.60235456931668878246228336815716196359 }, {
//-0.93502943687035390432897012004314760659,
//-0.73229569557598097380434312807086675041,
//-0.93502943687035390432897012004314760659 }, {
//-0.93502943687035390432897012004314760659,
//-0.93502943687035390432897012004314760659,
//-0.73229569557598097380434312807086675041 }, {
//-0.93502943687035390432897012004314760659,
// 0.60235456931668878246228336815716196359,
//-0.73229569557598097380434312807086675041 }, {
//-0.73229569557598097380434312807086675041,
//-0.93502943687035390432897012004314760659,
//-0.93502943687035390432897012004314760659 }, {
// 0.60235456931668878246228336815716196359,
//-0.73229569557598097380434312807086675041,
//-0.93502943687035390432897012004314760659 }, {
// 0.2561436909507320213865521444358193313,
//-0.65004131563212195143010154694337920556,
//-0.95606105968648811852634905054906092018 }, {
// 0.2561436909507320213865521444358193313,
//-0.65004131563212195143010154694337920556,
//-0.65004131563212195143010154694337920556 }, {
//-0.65004131563212195143010154694337920556,
//-0.65004131563212195143010154694337920556,
// 0.2561436909507320213865521444358193313 }, {
//-0.95606105968648811852634905054906092018,
// 0.2561436909507320213865521444358193313,
//-0.65004131563212195143010154694337920556 }, {
//-0.65004131563212195143010154694337920556,
//-0.95606105968648811852634905054906092018,
// 0.2561436909507320213865521444358193313 }, {
//-0.65004131563212195143010154694337920556,
// 0.2561436909507320213865521444358193313,
//-0.65004131563212195143010154694337920556 }, {
//-0.95606105968648811852634905054906092018,
//-0.65004131563212195143010154694337920556,
// 0.2561436909507320213865521444358193313 }, {
//-0.65004131563212195143010154694337920556,
//-0.95606105968648811852634905054906092018,
//-0.65004131563212195143010154694337920556 }, {
//-0.65004131563212195143010154694337920556,
//-0.65004131563212195143010154694337920556,
//-0.95606105968648811852634905054906092018 }, {
//-0.65004131563212195143010154694337920556,
// 0.2561436909507320213865521444358193313,
//-0.95606105968648811852634905054906092018 }, {
//-0.95606105968648811852634905054906092018,
//-0.65004131563212195143010154694337920556,
//-0.65004131563212195143010154694337920556 }, {
// 0.2561436909507320213865521444358193313,
//-0.95606105968648811852634905054906092018,
//-0.65004131563212195143010154694337920556 } } , {
// 0.063199698074694317846318428237401466134,
// 0.035916079989691599737018881339842778615,
// 0.035916079989691599737018881339842778615,
// 0.035916079989691599737018881339842778615,
// 0.035916079989691599737018881339842778615,
// 0.013158879622391177646076980573564101693,
// 0.013158879622391177646076980573564101693,
// 0.013158879622391177646076980573564101693,
// 0.013158879622391177646076980573564101693,
// 0.015191841626926975498161246507619114195,
// 0.015191841626926975498161246507619114195,
// 0.015191841626926975498161246507619114195,
// 0.015191841626926975498161246507619114195,
// 0.015191841626926975498161246507619114195,
// 0.015191841626926975498161246507619114195,
// 0.015191841626926975498161246507619114195,
// 0.015191841626926975498161246507619114195,
// 0.015191841626926975498161246507619114195,
// 0.015191841626926975498161246507619114195,
// 0.015191841626926975498161246507619114195,
// 0.015191841626926975498161246507619114195,
// 0.00048259245911900483231983784641134908616,
// 0.00048259245911900483231983784641134908616,
// 0.00048259245911900483231983784641134908616,
// 0.00048259245911900483231983784641134908616,
// 0.00048259245911900483231983784641134908616,
// 0.00048259245911900483231983784641134908616,
// 0.00048259245911900483231983784641134908616,
// 0.00048259245911900483231983784641134908616,
// 0.00048259245911900483231983784641134908616,
// 0.00048259245911900483231983784641134908616,
// 0.00048259245911900483231983784641134908616,
// 0.00048259245911900483231983784641134908616,
// 0.034319642640608095038714683012872940214,
// 0.034319642640608095038714683012872940214,
// 0.034319642640608095038714683012872940214,
// 0.034319642640608095038714683012872940214,
// 0.034319642640608095038714683012872940214,
// 0.034319642640608095038714683012872940214,
// 0.034319642640608095038714683012872940214,
// 0.034319642640608095038714683012872940214,
// 0.034319642640608095038714683012872940214,
// 0.034319642640608095038714683012872940214,
// 0.034319642640608095038714683012872940214,
// 0.034319642640608095038714683012872940214,
// 0.013514495573007723718021960153355702928,
// 0.013514495573007723718021960153355702928,
// 0.013514495573007723718021960153355702928,
// 0.013514495573007723718021960153355702928,
// 0.013514495573007723718021960153355702928,
// 0.013514495573007723718021960153355702928,
// 0.013514495573007723718021960153355702928,
// 0.013514495573007723718021960153355702928,
// 0.013514495573007723718021960153355702928,
// 0.013514495573007723718021960153355702928,
// 0.013514495573007723718021960153355702928,
// 0.013514495573007723718021960153355702928,
// 0.0087681963693812055566076536006009347954,
// 0.0087681963693812055566076536006009347954,
// 0.0087681963693812055566076536006009347954,
// 0.0087681963693812055566076536006009347954,
// 0.0087681963693812055566076536006009347954,
// 0.0087681963693812055566076536006009347954,
// 0.0087681963693812055566076536006009347954,
// 0.0087681963693812055566076536006009347954,
// 0.0087681963693812055566076536006009347954,
// 0.0087681963693812055566076536006009347954,
// 0.0087681963693812055566076536006009347954,
// 0.0087681963693812055566076536006009347954,
// 0.017209381065149320852393906999331987612,
// 0.017209381065149320852393906999331987612,
// 0.017209381065149320852393906999331987612,
// 0.017209381065149320852393906999331987612,
// 0.017209381065149320852393906999331987612,
// 0.017209381065149320852393906999331987612,
// 0.017209381065149320852393906999331987612,
// 0.017209381065149320852393906999331987612,
// 0.017209381065149320852393906999331987612,
// 0.017209381065149320852393906999331987612,
// 0.017209381065149320852393906999331987612,
// 0.017209381065149320852393906999331987612 } };
// default: 						throw
// std::runtime_error("not supported order");
//				}
//			}
//			case Figure::hexahedral: {
//				switch (integrand_order) {
//					case 0:
//					case 1:		return { { { 0, 0, 0 } }
//,
//{
// 8 } }; 					case 2:
// case 3: return { { { -0.577350269189626, -0.577350269189626,
// -0.577350269189626 },
// {
//-0.577350269189626, -0.577350269189626, 0.577350269189626 }, {
//-0.577350269189626, 0.577350269189626, -0.577350269189626 }, {
// -0.577350269189626, 0.577350269189626, 0.577350269189626 }, {
// 0.577350269189626, -0.577350269189626, -0.577350269189626 }, {
// 0.577350269189626, -0.577350269189626, 0.577350269189626 }, {
// 0.577350269189626, 0.577350269189626, -0.577350269189626 }, {
// 0.577350269189626, 0.577350269189626, 0.577350269189626 } } , { 1, 1, 1, 1,
// 1, 1, 1, 1 } }; 					case 4:
// case 5:		return { { { -0.774596669241483, -0.774596669241483,
//-0.774596669241483 }, { -0.774596669241483, -0.774596669241483, 0 }, {
//-0.774596669241483, -0.774596669241483, 0.774596669241483 }, {
//-0.774596669241483, 0, -0.774596669241483 }, { -0.774596669241483, 0, 0 }, {
//-0.774596669241483, 0, 0.774596669241483 }, { -0.774596669241483,
// 0.774596669241483, -0.774596669241483 }, { -0.774596669241483,
// 0.774596669241483, 0 }, { -0.774596669241483, 0.774596669241483,
// 0.774596669241483 }, { 0, -0.774596669241483, -0.774596669241483 }, { 0,
//-0.774596669241483, 0 }, { 0, -0.774596669241483, 0.774596669241483 }, { 0, 0,
//-0.774596669241483 }, { 0, 0, 0 }, { 0, 0, 0.774596669241483 }, { 0,
// 0.774596669241483, -0.774596669241483 }, { 0, 0.774596669241483, 0 }, { 0,
// 0.774596669241483, 0.774596669241483 }, { 0.774596669241483,
//-0.774596669241483, -0.774596669241483 }, { 0.774596669241483,
//-0.774596669241483, 0 }, { 0.774596669241483, -0.774596669241483,
// 0.774596669241483 }, { 0.774596669241483, 0, -0.774596669241483 }, {
// 0.774596669241483, 0, 0 }, { 0.774596669241483, 0, 0.774596669241483 }, {
// 0.774596669241483, 0.774596669241483, -0.774596669241483 }, {
// 0.774596669241483, 0.774596669241483, 0 }, { 0.774596669241483,
// 0.774596669241483, 0.774596669241483 } } , { 0.1714677640603552,
// 0.2743484224965692, 0.1714677640603552, 0.2743484224965692,
// 0.4389574759945119, 0.2743484224965692, 0.1714677640603552,
// 0.2743484224965692, 0.1714677640603552, 0.2743484224965692,
// 0.4389574759945119, 0.2743484224965692, 0.4389574759945118,
// 0.7023319615912209, 0.4389574759945118, 0.2743484224965692,
// 0.4389574759945119, 0.2743484224965692, 0.1714677640603552,
// 0.2743484224965692, 0.1714677640603552, 0.2743484224965692,
// 0.4389574759945119, 0.2743484224965692, 0.1714677640603552,
// 0.2743484224965692, 0.1714677640603552 } }; case
// 6: 					case 7:		return { { {
// -0.861136311594052, -0.861136311594052,
//-0.861136311594052 }, { -0.861136311594052, -0.861136311594052,
//-0.339981043584856 }, { -0.861136311594052, -0.861136311594052,
// 0.339981043584856 }, { -0.861136311594052, -0.861136311594052,
// 0.861136311594052 }, { -0.861136311594052, -0.339981043584856,
//-0.861136311594052 }, { -0.861136311594052, -0.339981043584856,
//-0.339981043584856 }, { -0.861136311594052, -0.339981043584856,
// 0.339981043584856 }, { -0.861136311594052, -0.339981043584856,
// 0.861136311594052 }, { -0.861136311594052, 0.339981043584856,
//-0.861136311594052 }, { -0.861136311594052, 0.339981043584856,
//-0.339981043584856 }, { -0.861136311594052, 0.339981043584856,
// 0.339981043584856 }, { -0.861136311594052, 0.339981043584856,
// 0.861136311594052 }, { -0.861136311594052, 0.861136311594052,
//-0.861136311594052 }, { -0.861136311594052, 0.861136311594052,
//-0.339981043584856 }, { -0.861136311594052, 0.861136311594052,
// 0.339981043584856 }, { -0.861136311594052, 0.861136311594052,
// 0.861136311594052 }, { -0.339981043584856, -0.861136311594052,
//-0.861136311594052 }, { -0.339981043584856, -0.861136311594052,
//-0.339981043584856 }, { -0.339981043584856, -0.861136311594052,
// 0.339981043584856 }, { -0.339981043584856, -0.861136311594052,
// 0.861136311594052 }, { -0.339981043584856, -0.339981043584856,
//-0.861136311594052 }, { -0.339981043584856, -0.339981043584856,
//-0.339981043584856 }, { -0.339981043584856, -0.339981043584856,
// 0.339981043584856 }, { -0.339981043584856, -0.339981043584856,
// 0.861136311594052 }, { -0.339981043584856, 0.339981043584856,
//-0.861136311594052 }, { -0.339981043584856, 0.339981043584856,
//-0.339981043584856 }, { -0.339981043584856, 0.339981043584856,
// 0.339981043584856 }, { -0.339981043584856, 0.339981043584856,
// 0.861136311594052 }, { -0.339981043584856, 0.861136311594052,
//-0.861136311594052 }, { -0.339981043584856, 0.861136311594052,
//-0.339981043584856 }, { -0.339981043584856, 0.861136311594052,
// 0.339981043584856 }, { -0.339981043584856, 0.861136311594052,
// 0.861136311594052 }, { 0.339981043584856, -0.861136311594052,
// -0.861136311594052 }, { 0.339981043584856, -0.861136311594052,
// -0.339981043584856 }, { 0.339981043584856, -0.861136311594052,
// 0.339981043584856 }, { 0.339981043584856, -0.861136311594052,
// 0.861136311594052 }, { 0.339981043584856, -0.339981043584856,
// -0.861136311594052 }, { 0.339981043584856, -0.339981043584856,
// -0.339981043584856 }, { 0.339981043584856, -0.339981043584856,
// 0.339981043584856 }, { 0.339981043584856, -0.339981043584856,
// 0.861136311594052 }, { 0.339981043584856, 0.339981043584856,
// -0.861136311594052 }, { 0.339981043584856, 0.339981043584856,
// -0.339981043584856 }, { 0.339981043584856, 0.339981043584856,
// 0.339981043584856 }, { 0.339981043584856, 0.339981043584856,
// 0.861136311594052 }, { 0.339981043584856, 0.861136311594052,
// -0.861136311594052 }, { 0.339981043584856, 0.861136311594052,
// -0.339981043584856 }, { 0.339981043584856, 0.861136311594052,
// 0.339981043584856 }, { 0.339981043584856, 0.861136311594052,
// 0.861136311594052 }, { 0.861136311594052, -0.861136311594052,
// -0.861136311594052 }, { 0.861136311594052, -0.861136311594052,
// -0.339981043584856 }, { 0.861136311594052, -0.861136311594052,
// 0.339981043584856 }, { 0.861136311594052, -0.861136311594052,
// 0.861136311594052 }, { 0.861136311594052, -0.339981043584856,
// -0.861136311594052 }, { 0.861136311594052, -0.339981043584856,
// -0.339981043584856 }, { 0.861136311594052, -0.339981043584856,
// 0.339981043584856 }, { 0.861136311594052, -0.339981043584856,
// 0.861136311594052 }, { 0.861136311594052, 0.339981043584856,
// -0.861136311594052 }, { 0.861136311594052, 0.339981043584856,
// -0.339981043584856 }, { 0.861136311594052, 0.339981043584856,
// 0.339981043584856 }, { 0.861136311594052, 0.339981043584856,
// 0.861136311594052 }, { 0.861136311594052, 0.861136311594052,
// -0.861136311594052 }, { 0.861136311594052, 0.861136311594052,
// -0.339981043584856 }, { 0.861136311594052, 0.861136311594052,
// 0.339981043584856 }, { 0.861136311594052, 0.861136311594052,
// 0.861136311594052 } } , { 0.04209147749053151, 0.07891151579507061,
// 0.07891151579507061, 0.04209147749053151, 0.07891151579507059,
// 0.1479403360567813, 0.1479403360567813, 0.07891151579507059,
// 0.07891151579507059, 0.1479403360567813, 0.1479403360567813,
// 0.07891151579507059, 0.04209147749053151, 0.07891151579507061,
// 0.07891151579507061, 0.04209147749053151, 0.07891151579507059,
// 0.1479403360567813, 0.1479403360567813, 0.07891151579507059,
// 0.1479403360567813, 0.2773529669539128, 0.2773529669539128,
// 0.1479403360567813, 0.1479403360567813, 0.2773529669539128,
// 0.2773529669539128, 0.1479403360567813, 0.07891151579507059,
// 0.1479403360567813, 0.1479403360567813, 0.07891151579507059,
// 0.07891151579507059, 0.1479403360567813, 0.1479403360567813,
// 0.07891151579507059, 0.1479403360567813, 0.2773529669539128,
// 0.2773529669539128, 0.1479403360567813, 0.1479403360567813,
// 0.2773529669539128, 0.2773529669539128, 0.1479403360567813,
// 0.07891151579507059, 0.1479403360567813, 0.1479403360567813,
// 0.07891151579507059, 0.04209147749053151, 0.07891151579507061,
// 0.07891151579507061, 0.04209147749053151, 0.07891151579507059,
// 0.1479403360567813, 0.1479403360567813, 0.07891151579507059,
// 0.07891151579507059, 0.1479403360567813, 0.1479403360567813,
// 0.07891151579507059, 0.04209147749053151, 0.07891151579507061,
// 0.07891151579507061, 0.04209147749053151 } };
// case 8: 					case 9:		return { { {
//-0.906179845938664, -0.906179845938664, -0.906179845938664 }, {
//-0.906179845938664, -0.906179845938664, -0.538469310105683 }, {
//-0.906179845938664, -0.906179845938664, 0 }, { -0.906179845938664,
//-0.906179845938664, 0.538469310105683 }, { -0.906179845938664,
//-0.906179845938664, 0.906179845938664 }, { -0.906179845938664,
//-0.538469310105683, -0.906179845938664 }, { -0.906179845938664,
//-0.538469310105683, -0.538469310105683 }, { -0.906179845938664,
//-0.538469310105683, 0 }, { -0.906179845938664, -0.538469310105683,
// 0.538469310105683 }, { -0.906179845938664, -0.538469310105683,
// 0.906179845938664 }, { -0.906179845938664, 0, -0.906179845938664 }, {
//-0.906179845938664, 0, -0.538469310105683 }, { -0.906179845938664, 0, 0 }, {
//-0.906179845938664, 0, 0.538469310105683 }, { -0.906179845938664, 0,
// 0.906179845938664 }, { -0.906179845938664, 0.538469310105683,
//-0.906179845938664 }, { -0.906179845938664, 0.538469310105683,
//-0.538469310105683 }, { -0.906179845938664, 0.538469310105683, 0 }, {
//-0.906179845938664, 0.538469310105683, 0.538469310105683 }, {
//-0.906179845938664, 0.538469310105683, 0.906179845938664 }, {
//-0.906179845938664, 0.906179845938664, -0.906179845938664 }, {
//-0.906179845938664, 0.906179845938664, -0.538469310105683 }, {
//-0.906179845938664, 0.906179845938664, 0 }, { -0.906179845938664,
// 0.906179845938664, 0.538469310105683 }, { -0.906179845938664,
// 0.906179845938664, 0.906179845938664 }, { -0.538469310105683,
//-0.906179845938664, -0.906179845938664 }, { -0.538469310105683,
//-0.906179845938664, -0.538469310105683 }, { -0.538469310105683,
//-0.906179845938664, 0 }, { -0.538469310105683, -0.906179845938664,
// 0.538469310105683 }, { -0.538469310105683, -0.906179845938664,
// 0.906179845938664 }, { -0.538469310105683, -0.538469310105683,
//-0.906179845938664 }, { -0.538469310105683, -0.538469310105683,
//-0.538469310105683 }, { -0.538469310105683, -0.538469310105683, 0 }, {
//-0.538469310105683, -0.538469310105683, 0.538469310105683 }, {
//-0.538469310105683, -0.538469310105683, 0.906179845938664 }, {
//-0.538469310105683, 0, -0.906179845938664 }, { -0.538469310105683, 0,
//-0.538469310105683 }, { -0.538469310105683, 0, 0 }, { -0.538469310105683, 0,
// 0.538469310105683 }, { -0.538469310105683, 0, 0.906179845938664 }, {
//-0.538469310105683, 0.538469310105683, -0.906179845938664 }, {
//-0.538469310105683, 0.538469310105683, -0.538469310105683 }, {
//-0.538469310105683, 0.538469310105683, 0 }, { -0.538469310105683,
// 0.538469310105683, 0.538469310105683 }, { -0.538469310105683,
// 0.538469310105683, 0.906179845938664 }, { -0.538469310105683,
// 0.906179845938664, -0.906179845938664 }, { -0.538469310105683,
// 0.906179845938664, -0.538469310105683 }, { -0.538469310105683,
// 0.906179845938664, 0 }, { -0.538469310105683, 0.906179845938664,
// 0.538469310105683 }, { -0.538469310105683, 0.906179845938664,
// 0.906179845938664 }, { 0, -0.906179845938664, -0.906179845938664 }, { 0,
//-0.906179845938664, -0.538469310105683 }, { 0, -0.906179845938664, 0 }, { 0,
//-0.906179845938664, 0.538469310105683 }, { 0, -0.906179845938664,
// 0.906179845938664 }, { 0, -0.538469310105683, -0.906179845938664 }, { 0,
//-0.538469310105683, -0.538469310105683 }, { 0, -0.538469310105683, 0 }, { 0,
//-0.538469310105683, 0.538469310105683 }, { 0, -0.538469310105683,
// 0.906179845938664 }, { 0, 0, -0.906179845938664 }, { 0, 0, -0.538469310105683
//}, { 0, 0, 0 }, { 0, 0, 0.538469310105683 }, { 0, 0, 0.906179845938664 }, { 0,
// 0.538469310105683, -0.906179845938664 }, { 0, 0.538469310105683,
//-0.538469310105683 }, { 0, 0.538469310105683, 0 }, { 0, 0.538469310105683,
// 0.538469310105683 }, { 0, 0.538469310105683, 0.906179845938664 }, { 0,
// 0.906179845938664, -0.906179845938664 }, { 0, 0.906179845938664,
//-0.538469310105683 }, { 0, 0.906179845938664, 0 }, { 0, 0.906179845938664,
// 0.538469310105683 }, { 0, 0.906179845938664, 0.906179845938664 }, {
// 0.538469310105683, -0.906179845938664, -0.906179845938664 }, {
// 0.538469310105683, -0.906179845938664, -0.538469310105683 }, {
// 0.538469310105683, -0.906179845938664, 0 }, { 0.538469310105683,
//-0.906179845938664, 0.538469310105683 }, { 0.538469310105683,
//-0.906179845938664, 0.906179845938664 }, { 0.538469310105683,
//-0.538469310105683, -0.906179845938664 }, { 0.538469310105683,
//-0.538469310105683, -0.538469310105683 }, { 0.538469310105683,
//-0.538469310105683, 0 }, { 0.538469310105683, -0.538469310105683,
// 0.538469310105683 }, { 0.538469310105683, -0.538469310105683,
// 0.906179845938664 }, { 0.538469310105683, 0, -0.906179845938664 }, {
// 0.538469310105683, 0, -0.538469310105683 }, { 0.538469310105683, 0, 0 }, {
// 0.538469310105683, 0, 0.538469310105683 }, { 0.538469310105683, 0,
// 0.906179845938664 }, { 0.538469310105683, 0.538469310105683,
//-0.906179845938664 }, { 0.538469310105683, 0.538469310105683,
//-0.538469310105683 }, { 0.538469310105683, 0.538469310105683, 0 }, {
// 0.538469310105683, 0.538469310105683, 0.538469310105683 }, {
// 0.538469310105683, 0.538469310105683, 0.906179845938664 }, {
// 0.538469310105683, 0.906179845938664, -0.906179845938664 }, {
// 0.538469310105683, 0.906179845938664, -0.538469310105683 }, {
// 0.538469310105683, 0.906179845938664, 0 }, { 0.538469310105683,
// 0.906179845938664, 0.538469310105683 }, { 0.538469310105683,
// 0.906179845938664, 0.906179845938664 }, { 0.906179845938664,
//-0.906179845938664, -0.906179845938664 }, { 0.906179845938664,
//-0.906179845938664, -0.538469310105683 }, { 0.906179845938664,
//-0.906179845938664, 0 }, { 0.906179845938664, -0.906179845938664,
// 0.538469310105683 }, { 0.906179845938664, -0.906179845938664,
// 0.906179845938664 }, { 0.906179845938664, -0.538469310105683,
//-0.906179845938664 }, { 0.906179845938664, -0.538469310105683,
//-0.538469310105683 }, { 0.906179845938664, -0.538469310105683, 0 }, {
// 0.906179845938664, -0.538469310105683, 0.538469310105683 }, {
// 0.906179845938664, -0.538469310105683, 0.906179845938664 }, {
// 0.906179845938664, 0, -0.906179845938664 }, { 0.906179845938664, 0,
//-0.538469310105683 }, { 0.906179845938664, 0, 0 }, { 0.906179845938664, 0,
// 0.538469310105683 }, { 0.906179845938664, 0, 0.906179845938664 }, {
// 0.906179845938664, 0.538469310105683, -0.906179845938664 }, {
// 0.906179845938664, 0.538469310105683, -0.538469310105683 }, {
// 0.906179845938664, 0.538469310105683, 0 }, { 0.906179845938664,
// 0.538469310105683, 0.538469310105683 }, { 0.906179845938664,
// 0.538469310105683, 0.906179845938664 }, { 0.906179845938664,
// 0.906179845938664, -0.906179845938664 }, { 0.906179845938664,
// 0.906179845938664, -0.538469310105683 }, { 0.906179845938664,
// 0.906179845938664, 0 }, { 0.906179845938664, 0.906179845938664,
// 0.538469310105683 }, { 0.906179845938664, 0.906179845938664,
// 0.906179845938664 } } , { 0.01329973642063263, 0.0268675087653718,
// 0.03193420735284827, 0.0268675087653718, 0.01329973642063263,
// 0.0268675087653718, 0.05427649123462804, 0.06451199999999993,
// 0.05427649123462804, 0.0268675087653718, 0.03193420735284828,
// 0.06451199999999993, 0.07667773006934522, 0.06451199999999993,
// 0.03193420735284828, 0.0268675087653718, 0.05427649123462804,
// 0.06451199999999993, 0.05427649123462804, 0.0268675087653718,
// 0.01329973642063263, 0.0268675087653718, 0.03193420735284827,
// 0.0268675087653718, 0.01329973642063263, 0.0268675087653718,
// 0.05427649123462804, 0.06451199999999993, 0.05427649123462804,
// 0.0268675087653718, 0.05427649123462804, 0.1096468424545385,
// 0.1303241410696481, 0.1096468424545385, 0.05427649123462804,
// 0.06451199999999993, 0.1303241410696481, 0.1549007829622048,
// 0.1303241410696481, 0.06451199999999993, 0.05427649123462804,
// 0.1096468424545385, 0.1303241410696481, 0.1096468424545385,
// 0.05427649123462804, 0.0268675087653718, 0.05427649123462804,
// 0.06451199999999993, 0.05427649123462804, 0.0268675087653718,
// 0.03193420735284828, 0.06451199999999993, 0.07667773006934522,
// 0.06451199999999993, 0.03193420735284828, 0.06451199999999993,
// 0.1303241410696481, 0.1549007829622048, 0.1303241410696481,
// 0.06451199999999993, 0.07667773006934522, 0.1549007829622047,
// 0.1841121097393691, 0.1549007829622047, 0.07667773006934522,
// 0.06451199999999993, 0.1303241410696481, 0.1549007829622048,
// 0.1303241410696481, 0.06451199999999993, 0.03193420735284828,
// 0.06451199999999993, 0.07667773006934522, 0.06451199999999993,
// 0.03193420735284828, 0.0268675087653718, 0.05427649123462804,
// 0.06451199999999993, 0.05427649123462804, 0.0268675087653718,
// 0.05427649123462804, 0.1096468424545385, 0.1303241410696481,
// 0.1096468424545385, 0.05427649123462804, 0.06451199999999993,
// 0.1303241410696481, 0.1549007829622048, 0.1303241410696481,
// 0.06451199999999993, 0.05427649123462804, 0.1096468424545385,
// 0.1303241410696481, 0.1096468424545385, 0.05427649123462804,
// 0.0268675087653718, 0.05427649123462804, 0.06451199999999993,
// 0.05427649123462804, 0.0268675087653718, 0.01329973642063263,
// 0.0268675087653718, 0.03193420735284827, 0.0268675087653718,
// 0.01329973642063263, 0.0268675087653718, 0.05427649123462804,
// 0.06451199999999993, 0.05427649123462804, 0.0268675087653718,
// 0.03193420735284828, 0.06451199999999993, 0.07667773006934522,
// 0.06451199999999993, 0.03193420735284828, 0.0268675087653718,
// 0.05427649123462804, 0.06451199999999993, 0.05427649123462804,
// 0.0268675087653718, 0.01329973642063263, 0.0268675087653718,
// 0.03193420735284827, 0.0268675087653718, 0.01329973642063263 } };
// case 10: 					case 11:	return { { {
//-0.9324695142031521, -0.9324695142031521, -0.9324695142031521 }, {
//-0.9324695142031521, -0.9324695142031521, -0.661209386466264 }, {
//-0.9324695142031521, -0.9324695142031521, -0.238619186083197 }, {
//-0.9324695142031521, -0.9324695142031521, 0.238619186083197 }, {
//-0.9324695142031521, -0.9324695142031521, 0.661209386466264 }, {
//-0.9324695142031521, -0.9324695142031521, 0.9324695142031521 }, {
//-0.9324695142031521, -0.661209386466264, -0.9324695142031521 }, {
//-0.9324695142031521, -0.661209386466264, -0.661209386466264 }, {
//-0.9324695142031521, -0.661209386466264, -0.238619186083197 }, {
//-0.9324695142031521, -0.661209386466264, 0.238619186083197 }, {
//-0.9324695142031521, -0.661209386466264, 0.661209386466264 }, {
//-0.9324695142031521, -0.661209386466264, 0.9324695142031521 }, {
//-0.9324695142031521, -0.238619186083197, -0.9324695142031521 }, {
//-0.9324695142031521, -0.238619186083197, -0.661209386466264 }, {
//-0.9324695142031521, -0.238619186083197, -0.238619186083197 }, {
//-0.9324695142031521, -0.238619186083197, 0.238619186083197 }, {
//-0.9324695142031521, -0.238619186083197, 0.661209386466264 }, {
//-0.9324695142031521, -0.238619186083197, 0.9324695142031521 }, {
//-0.9324695142031521, 0.238619186083197, -0.9324695142031521 }, {
//-0.9324695142031521, 0.238619186083197, -0.661209386466264 }, {
//-0.9324695142031521, 0.238619186083197, -0.238619186083197 }, {
//-0.9324695142031521, 0.238619186083197, 0.238619186083197 }, {
//-0.9324695142031521, 0.238619186083197, 0.661209386466264 }, {
//-0.9324695142031521, 0.238619186083197, 0.9324695142031521 }, {
//-0.9324695142031521, 0.661209386466264, -0.9324695142031521 }, {
//-0.9324695142031521, 0.661209386466264, -0.661209386466264 }, {
//-0.9324695142031521, 0.661209386466264, -0.238619186083197 }, {
//-0.9324695142031521, 0.661209386466264, 0.238619186083197 }, {
//-0.9324695142031521, 0.661209386466264, 0.661209386466264 }, {
//-0.9324695142031521, 0.661209386466264, 0.9324695142031521 }, {
//-0.9324695142031521, 0.9324695142031521, -0.9324695142031521 }, {
//-0.9324695142031521, 0.9324695142031521, -0.661209386466264 }, {
//-0.9324695142031521, 0.9324695142031521, -0.238619186083197 }, {
//-0.9324695142031521, 0.9324695142031521, 0.238619186083197 }, {
//-0.9324695142031521, 0.9324695142031521, 0.661209386466264 }, {
//-0.9324695142031521, 0.9324695142031521, 0.9324695142031521 }, {
//-0.661209386466264, -0.9324695142031521, -0.9324695142031521 }, {
//-0.661209386466264, -0.9324695142031521, -0.661209386466264 }, {
//-0.661209386466264, -0.9324695142031521, -0.238619186083197 }, {
//-0.661209386466264, -0.9324695142031521, 0.238619186083197 }, {
//-0.661209386466264, -0.9324695142031521, 0.661209386466264 }, {
//-0.661209386466264, -0.9324695142031521, 0.9324695142031521 }, {
//-0.661209386466264, -0.661209386466264, -0.9324695142031521 }, {
//-0.661209386466264, -0.661209386466264, -0.661209386466264 }, {
//-0.661209386466264, -0.661209386466264, -0.238619186083197 }, {
//-0.661209386466264, -0.661209386466264, 0.238619186083197 }, {
//-0.661209386466264, -0.661209386466264, 0.661209386466264 }, {
//-0.661209386466264, -0.661209386466264, 0.9324695142031521 }, {
//-0.661209386466264, -0.238619186083197, -0.9324695142031521 }, {
//-0.661209386466264, -0.238619186083197, -0.661209386466264 }, {
//-0.661209386466264, -0.238619186083197, -0.238619186083197 }, {
//-0.661209386466264, -0.238619186083197, 0.238619186083197 }, {
//-0.661209386466264, -0.238619186083197, 0.661209386466264 }, {
//-0.661209386466264, -0.238619186083197, 0.9324695142031521 }, {
//-0.661209386466264, 0.238619186083197, -0.9324695142031521 }, {
//-0.661209386466264, 0.238619186083197, -0.661209386466264 }, {
//-0.661209386466264, 0.238619186083197, -0.238619186083197 }, {
//-0.661209386466264, 0.238619186083197, 0.238619186083197 }, {
//-0.661209386466264, 0.238619186083197, 0.661209386466264 }, {
//-0.661209386466264, 0.238619186083197, 0.9324695142031521 }, {
//-0.661209386466264, 0.661209386466264, -0.9324695142031521 }, {
//-0.661209386466264, 0.661209386466264, -0.661209386466264 }, {
//-0.661209386466264, 0.661209386466264, -0.238619186083197 }, {
//-0.661209386466264, 0.661209386466264, 0.238619186083197 }, {
//-0.661209386466264, 0.661209386466264, 0.661209386466264 }, {
//-0.661209386466264, 0.661209386466264, 0.9324695142031521 }, {
//-0.661209386466264, 0.9324695142031521, -0.9324695142031521 }, {
//-0.661209386466264, 0.9324695142031521, -0.661209386466264 }, {
//-0.661209386466264, 0.9324695142031521, -0.238619186083197 }, {
//-0.661209386466264, 0.9324695142031521, 0.238619186083197 }, {
//-0.661209386466264, 0.9324695142031521, 0.661209386466264 }, {
//-0.661209386466264, 0.9324695142031521, 0.9324695142031521 }, {
//-0.238619186083197, -0.9324695142031521, -0.9324695142031521 }, {
//-0.238619186083197, -0.9324695142031521, -0.661209386466264 }, {
//-0.238619186083197, -0.9324695142031521, -0.238619186083197 }, {
//-0.238619186083197, -0.9324695142031521, 0.238619186083197 }, {
//-0.238619186083197, -0.9324695142031521, 0.661209386466264 }, {
//-0.238619186083197, -0.9324695142031521, 0.9324695142031521 }, {
//-0.238619186083197, -0.661209386466264, -0.9324695142031521 }, {
//-0.238619186083197, -0.661209386466264, -0.661209386466264 }, {
//-0.238619186083197, -0.661209386466264, -0.238619186083197 }, {
//-0.238619186083197, -0.661209386466264, 0.238619186083197 }, {
//-0.238619186083197, -0.661209386466264, 0.661209386466264 }, {
//-0.238619186083197, -0.661209386466264, 0.9324695142031521 }, {
//-0.238619186083197, -0.238619186083197, -0.9324695142031521 }, {
//-0.238619186083197, -0.238619186083197, -0.661209386466264 }, {
//-0.238619186083197, -0.238619186083197, -0.238619186083197 }, {
//-0.238619186083197, -0.238619186083197, 0.238619186083197 }, {
//-0.238619186083197, -0.238619186083197, 0.661209386466264 }, {
//-0.238619186083197, -0.238619186083197, 0.9324695142031521 }, {
//-0.238619186083197, 0.238619186083197, -0.9324695142031521 }, {
//-0.238619186083197, 0.238619186083197, -0.661209386466264 }, {
//-0.238619186083197, 0.238619186083197, -0.238619186083197 }, {
//-0.238619186083197, 0.238619186083197, 0.238619186083197 }, {
//-0.238619186083197, 0.238619186083197, 0.661209386466264 }, {
//-0.238619186083197, 0.238619186083197, 0.9324695142031521 }, {
//-0.238619186083197, 0.661209386466264, -0.9324695142031521 }, {
//-0.238619186083197, 0.661209386466264, -0.661209386466264 }, {
//-0.238619186083197, 0.661209386466264, -0.238619186083197 }, {
//-0.238619186083197, 0.661209386466264, 0.238619186083197 }, {
//-0.238619186083197, 0.661209386466264, 0.661209386466264 }, {
//-0.238619186083197, 0.661209386466264, 0.9324695142031521 }, {
//-0.238619186083197, 0.9324695142031521, -0.9324695142031521 }, {
//-0.238619186083197, 0.9324695142031521, -0.661209386466264 }, {
//-0.238619186083197, 0.9324695142031521, -0.238619186083197 }, {
//-0.238619186083197, 0.9324695142031521, 0.238619186083197 }, {
//-0.238619186083197, 0.9324695142031521, 0.661209386466264 }, {
//-0.238619186083197, 0.9324695142031521, 0.9324695142031521 }, {
// 0.238619186083197, -0.9324695142031521, -0.9324695142031521 }, {
// 0.238619186083197, -0.9324695142031521, -0.661209386466264 }, {
// 0.238619186083197, -0.9324695142031521, -0.238619186083197 }, {
// 0.238619186083197, -0.9324695142031521, 0.238619186083197 }, {
// 0.238619186083197, -0.9324695142031521, 0.661209386466264 }, {
// 0.238619186083197, -0.9324695142031521, 0.9324695142031521 }, {
// 0.238619186083197, -0.661209386466264, -0.9324695142031521 }, {
// 0.238619186083197, -0.661209386466264, -0.661209386466264 }, {
// 0.238619186083197, -0.661209386466264, -0.238619186083197 }, {
// 0.238619186083197, -0.661209386466264, 0.238619186083197 }, {
// 0.238619186083197, -0.661209386466264, 0.661209386466264 }, {
// 0.238619186083197, -0.661209386466264, 0.9324695142031521 }, {
// 0.238619186083197, -0.238619186083197, -0.9324695142031521 }, {
// 0.238619186083197, -0.238619186083197, -0.661209386466264 }, {
// 0.238619186083197, -0.238619186083197, -0.238619186083197 }, {
// 0.238619186083197, -0.238619186083197, 0.238619186083197 }, {
// 0.238619186083197, -0.238619186083197, 0.661209386466264 }, {
// 0.238619186083197, -0.238619186083197, 0.9324695142031521 }, {
// 0.238619186083197, 0.238619186083197, -0.9324695142031521 }, {
// 0.238619186083197, 0.238619186083197, -0.661209386466264 }, {
// 0.238619186083197, 0.238619186083197, -0.238619186083197 }, {
// 0.238619186083197, 0.238619186083197, 0.238619186083197 }, {
// 0.238619186083197, 0.238619186083197, 0.661209386466264 }, {
// 0.238619186083197, 0.238619186083197, 0.9324695142031521 }, {
// 0.238619186083197, 0.661209386466264, -0.9324695142031521 }, {
// 0.238619186083197, 0.661209386466264, -0.661209386466264 }, {
// 0.238619186083197, 0.661209386466264, -0.238619186083197 }, {
// 0.238619186083197, 0.661209386466264, 0.238619186083197 }, {
// 0.238619186083197, 0.661209386466264, 0.661209386466264 }, {
// 0.238619186083197, 0.661209386466264, 0.9324695142031521 }, {
// 0.238619186083197, 0.9324695142031521, -0.9324695142031521 }, {
// 0.238619186083197, 0.9324695142031521, -0.661209386466264 }, {
// 0.238619186083197, 0.9324695142031521, -0.238619186083197 }, {
// 0.238619186083197, 0.9324695142031521, 0.238619186083197 }, {
// 0.238619186083197, 0.9324695142031521, 0.661209386466264 }, {
// 0.238619186083197, 0.9324695142031521, 0.9324695142031521 }, {
// 0.661209386466264, -0.9324695142031521, -0.9324695142031521 }, {
// 0.661209386466264, -0.9324695142031521, -0.661209386466264 }, {
// 0.661209386466264, -0.9324695142031521, -0.238619186083197 }, {
// 0.661209386466264, -0.9324695142031521, 0.238619186083197 }, {
// 0.661209386466264, -0.9324695142031521, 0.661209386466264 }, {
// 0.661209386466264, -0.9324695142031521, 0.9324695142031521 }, {
// 0.661209386466264, -0.661209386466264, -0.9324695142031521 }, {
// 0.661209386466264, -0.661209386466264, -0.661209386466264 }, {
// 0.661209386466264, -0.661209386466264, -0.238619186083197 }, {
// 0.661209386466264, -0.661209386466264, 0.238619186083197 }, {
// 0.661209386466264, -0.661209386466264, 0.661209386466264 }, {
// 0.661209386466264, -0.661209386466264, 0.9324695142031521 }, {
// 0.661209386466264, -0.238619186083197, -0.9324695142031521 }, {
// 0.661209386466264, -0.238619186083197, -0.661209386466264 }, {
// 0.661209386466264, -0.238619186083197, -0.238619186083197 }, {
// 0.661209386466264, -0.238619186083197, 0.238619186083197 }, {
// 0.661209386466264, -0.238619186083197, 0.661209386466264 }, {
// 0.661209386466264, -0.238619186083197, 0.9324695142031521 }, {
// 0.661209386466264, 0.238619186083197, -0.9324695142031521 }, {
// 0.661209386466264, 0.238619186083197, -0.661209386466264 }, {
// 0.661209386466264, 0.238619186083197, -0.238619186083197 }, {
// 0.661209386466264, 0.238619186083197, 0.238619186083197 }, {
// 0.661209386466264, 0.238619186083197, 0.661209386466264 }, {
// 0.661209386466264, 0.238619186083197, 0.9324695142031521 }, {
// 0.661209386466264, 0.661209386466264, -0.9324695142031521 }, {
// 0.661209386466264, 0.661209386466264, -0.661209386466264 }, {
// 0.661209386466264, 0.661209386466264, -0.238619186083197 }, {
// 0.661209386466264, 0.661209386466264, 0.238619186083197 }, {
// 0.661209386466264, 0.661209386466264, 0.661209386466264 }, {
// 0.661209386466264, 0.661209386466264, 0.9324695142031521 }, {
// 0.661209386466264, 0.9324695142031521, -0.9324695142031521 }, {
// 0.661209386466264, 0.9324695142031521, -0.661209386466264 }, {
// 0.661209386466264, 0.9324695142031521, -0.238619186083197 }, {
// 0.661209386466264, 0.9324695142031521, 0.238619186083197 }, {
// 0.661209386466264, 0.9324695142031521, 0.661209386466264 }, {
// 0.661209386466264, 0.9324695142031521, 0.9324695142031521 }, {
// 0.9324695142031521, -0.9324695142031521, -0.9324695142031521 }, {
// 0.9324695142031521, -0.9324695142031521, -0.661209386466264 }, {
// 0.9324695142031521, -0.9324695142031521, -0.238619186083197 }, {
// 0.9324695142031521, -0.9324695142031521, 0.238619186083197 }, {
// 0.9324695142031521, -0.9324695142031521, 0.661209386466264 }, {
// 0.9324695142031521, -0.9324695142031521, 0.9324695142031521 }, {
// 0.9324695142031521, -0.661209386466264, -0.9324695142031521 }, {
// 0.9324695142031521, -0.661209386466264, -0.661209386466264 }, {
// 0.9324695142031521, -0.661209386466264, -0.238619186083197 }, {
// 0.9324695142031521, -0.661209386466264, 0.238619186083197 }, {
// 0.9324695142031521, -0.661209386466264, 0.661209386466264 }, {
// 0.9324695142031521, -0.661209386466264, 0.9324695142031521 }, {
// 0.9324695142031521, -0.238619186083197, -0.9324695142031521 }, {
// 0.9324695142031521, -0.238619186083197, -0.661209386466264 }, {
// 0.9324695142031521, -0.238619186083197, -0.238619186083197 }, {
// 0.9324695142031521, -0.238619186083197, 0.238619186083197 }, {
// 0.9324695142031521, -0.238619186083197, 0.661209386466264 }, {
// 0.9324695142031521, -0.238619186083197, 0.9324695142031521 }, {
// 0.9324695142031521, 0.238619186083197, -0.9324695142031521 }, {
// 0.9324695142031521, 0.238619186083197, -0.661209386466264 }, {
// 0.9324695142031521, 0.238619186083197, -0.238619186083197 }, {
// 0.9324695142031521, 0.238619186083197, 0.238619186083197 }, {
// 0.9324695142031521, 0.238619186083197, 0.661209386466264 }, {
// 0.9324695142031521, 0.238619186083197, 0.9324695142031521 }, {
// 0.9324695142031521, 0.661209386466264, -0.9324695142031521 }, {
// 0.9324695142031521, 0.661209386466264, -0.661209386466264 }, {
// 0.9324695142031521, 0.661209386466264, -0.238619186083197 }, {
// 0.9324695142031521, 0.661209386466264, 0.238619186083197 }, {
// 0.9324695142031521, 0.661209386466264, 0.661209386466264 }, {
// 0.9324695142031521, 0.661209386466264, 0.9324695142031521 }, {
// 0.9324695142031521, 0.9324695142031521, -0.9324695142031521 }, {
// 0.9324695142031521, 0.9324695142031521, -0.661209386466264 }, {
// 0.9324695142031521, 0.9324695142031521, -0.238619186083197 }, {
// 0.9324695142031521, 0.9324695142031521, 0.238619186083197 }, {
// 0.9324695142031521, 0.9324695142031521, 0.661209386466264 }, {
// 0.9324695142031521, 0.9324695142031521, 0.9324695142031521 } } , {
// 0.005028730495636565, 0.01058910316235413, 0.01373424803098996,
// 0.01373424803098996, 0.01058910316235413, 0.005028730495636565,
// 0.01058910316235413, 0.02229769638286893, 0.02892049382716063,
// 0.02892049382716063, 0.02229769638286893, 0.01058910316235413,
// 0.01373424803098996, 0.02892049382716063, 0.0375103754596564,
// 0.0375103754596564, 0.02892049382716063, 0.01373424803098996,
// 0.01373424803098996, 0.02892049382716063, 0.0375103754596564,
// 0.0375103754596564, 0.02892049382716063, 0.01373424803098996,
// 0.01058910316235413, 0.02229769638286893, 0.02892049382716063,
// 0.02892049382716063, 0.02229769638286893, 0.01058910316235413,
// 0.005028730495636565, 0.01058910316235413, 0.01373424803098996,
// 0.01373424803098996, 0.01058910316235413, 0.005028730495636565,
// 0.01058910316235413, 0.02229769638286893, 0.02892049382716063,
// 0.02892049382716063, 0.02229769638286893, 0.01058910316235413,
// 0.02229769638286893, 0.04695272643581212, 0.06089848976948679,
// 0.06089848976948679, 0.04695272643581212, 0.02229769638286893,
// 0.02892049382716063, 0.06089848976948679, 0.07898638349094073,
// 0.07898638349094073, 0.06089848976948679, 0.02892049382716063,
// 0.02892049382716063, 0.06089848976948679, 0.07898638349094073,
// 0.07898638349094073, 0.06089848976948679, 0.02892049382716063,
// 0.02229769638286893, 0.04695272643581212, 0.06089848976948679,
// 0.06089848976948679, 0.04695272643581212, 0.02229769638286893,
// 0.01058910316235413, 0.02229769638286893, 0.02892049382716063,
// 0.02892049382716063, 0.02229769638286893, 0.01058910316235413,
// 0.01373424803098996, 0.02892049382716063, 0.0375103754596564,
// 0.0375103754596564, 0.02892049382716063, 0.01373424803098996,
// 0.02892049382716063, 0.06089848976948679, 0.07898638349094073,
// 0.07898638349094073, 0.06089848976948679, 0.02892049382716063,
// 0.0375103754596564, 0.07898638349094074, 0.1024466912166996,
// 0.1024466912166996, 0.07898638349094074, 0.0375103754596564,
// 0.0375103754596564, 0.07898638349094074, 0.1024466912166996,
// 0.1024466912166996, 0.07898638349094074, 0.0375103754596564,
// 0.02892049382716063, 0.06089848976948679, 0.07898638349094073,
// 0.07898638349094073, 0.06089848976948679, 0.02892049382716063,
// 0.01373424803098996, 0.02892049382716063, 0.0375103754596564,
// 0.0375103754596564, 0.02892049382716063, 0.01373424803098996,
// 0.01373424803098996, 0.02892049382716063, 0.0375103754596564,
// 0.0375103754596564, 0.02892049382716063, 0.01373424803098996,
// 0.02892049382716063, 0.06089848976948679, 0.07898638349094073,
// 0.07898638349094073, 0.06089848976948679, 0.02892049382716063,
// 0.0375103754596564, 0.07898638349094074, 0.1024466912166996,
// 0.1024466912166996, 0.07898638349094074, 0.0375103754596564,
// 0.0375103754596564, 0.07898638349094074, 0.1024466912166996,
// 0.1024466912166996, 0.07898638349094074, 0.0375103754596564,
// 0.02892049382716063, 0.06089848976948679, 0.07898638349094073,
// 0.07898638349094073, 0.06089848976948679, 0.02892049382716063,
// 0.01373424803098996, 0.02892049382716063, 0.0375103754596564,
// 0.0375103754596564, 0.02892049382716063, 0.01373424803098996,
// 0.01058910316235413, 0.02229769638286893, 0.02892049382716063,
// 0.02892049382716063, 0.02229769638286893, 0.01058910316235413,
// 0.02229769638286893, 0.04695272643581212, 0.06089848976948679,
// 0.06089848976948679, 0.04695272643581212, 0.02229769638286893,
// 0.02892049382716063, 0.06089848976948679, 0.07898638349094073,
// 0.07898638349094073, 0.06089848976948679, 0.02892049382716063,
// 0.02892049382716063, 0.06089848976948679, 0.07898638349094073,
// 0.07898638349094073, 0.06089848976948679, 0.02892049382716063,
// 0.02229769638286893, 0.04695272643581212, 0.06089848976948679,
// 0.06089848976948679, 0.04695272643581212, 0.02229769638286893,
// 0.01058910316235413, 0.02229769638286893, 0.02892049382716063,
// 0.02892049382716063, 0.02229769638286893, 0.01058910316235413,
// 0.005028730495636565, 0.01058910316235413, 0.01373424803098996,
// 0.01373424803098996, 0.01058910316235413, 0.005028730495636565,
// 0.01058910316235413, 0.02229769638286893, 0.02892049382716063,
// 0.02892049382716063, 0.02229769638286893, 0.01058910316235413,
// 0.01373424803098996, 0.02892049382716063, 0.0375103754596564,
// 0.0375103754596564, 0.02892049382716063, 0.01373424803098996,
// 0.01373424803098996, 0.02892049382716063, 0.0375103754596564,
// 0.0375103754596564, 0.02892049382716063, 0.01373424803098996,
// 0.01058910316235413, 0.02229769638286893, 0.02892049382716063,
// 0.02892049382716063, 0.02229769638286893, 0.01058910316235413,
// 0.005028730495636565, 0.01058910316235413, 0.01373424803098996,
// 0.01373424803098996, 0.01058910316235413, 0.005028730495636565 } };
// case 12: 					case 13:	return { { {
// -0.949107912342758, -0.949107912342758,
//-0.949107912342758 }, { -0.949107912342758, -0.949107912342758,
//-0.741531185599394 }, { -0.949107912342758, -0.949107912342758,
//-0.405845151377397 }, { -0.949107912342758, -0.949107912342758, 0 }, {
//-0.949107912342758, -0.949107912342758, 0.405845151377397 }, {
// -0.949107912342758, -0.949107912342758, 0.741531185599394 }, {
// -0.949107912342758, -0.949107912342758, 0.949107912342758 }, {
// -0.949107912342758, -0.741531185599394, -0.949107912342758 }, {
//-0.949107912342758, -0.741531185599394, -0.741531185599394 }, {
//-0.949107912342758, -0.741531185599394, -0.405845151377397 }, {
//-0.949107912342758, -0.741531185599394, 0 }, { -0.949107912342758,
//-0.741531185599394, 0.405845151377397 }, { -0.949107912342758,
//-0.741531185599394, 0.741531185599394 }, { -0.949107912342758,
//-0.741531185599394, 0.949107912342758 }, { -0.949107912342758,
//-0.405845151377397, -0.949107912342758 }, { -0.949107912342758,
//-0.405845151377397, -0.741531185599394 }, { -0.949107912342758,
//-0.405845151377397, -0.405845151377397 }, { -0.949107912342758,
//-0.405845151377397, 0 }, { -0.949107912342758, -0.405845151377397,
// 0.405845151377397 }, { -0.949107912342758, -0.405845151377397,
// 0.741531185599394 }, { -0.949107912342758, -0.405845151377397,
// 0.949107912342758 }, { -0.949107912342758, 0, -0.949107912342758 }, {
//-0.949107912342758, 0, -0.741531185599394 }, { -0.949107912342758, 0,
//-0.405845151377397 }, { -0.949107912342758, 0, 0 }, { -0.949107912342758, 0,
// 0.405845151377397 }, { -0.949107912342758, 0, 0.741531185599394 }, {
// -0.949107912342758, 0, 0.949107912342758 }, { -0.949107912342758,
// 0.405845151377397, -0.949107912342758 }, { -0.949107912342758,
// 0.405845151377397, -0.741531185599394 }, { -0.949107912342758,
// 0.405845151377397, -0.405845151377397 }, { -0.949107912342758,
// 0.405845151377397, 0 }, { -0.949107912342758, 0.405845151377397,
// 0.405845151377397 }, { -0.949107912342758, 0.405845151377397,
// 0.741531185599394 }, { -0.949107912342758, 0.405845151377397,
// 0.949107912342758 }, { -0.949107912342758, 0.741531185599394,
// -0.949107912342758 }, { -0.949107912342758, 0.741531185599394,
// -0.741531185599394 }, { -0.949107912342758, 0.741531185599394,
// -0.405845151377397 }, { -0.949107912342758, 0.741531185599394, 0 }, {
// -0.949107912342758, 0.741531185599394, 0.405845151377397 }, {
// -0.949107912342758, 0.741531185599394, 0.741531185599394 }, {
// -0.949107912342758, 0.741531185599394, 0.949107912342758 }, {
// -0.949107912342758, 0.949107912342758, -0.949107912342758 }, {
//-0.949107912342758, 0.949107912342758, -0.741531185599394 }, {
//-0.949107912342758, 0.949107912342758, -0.405845151377397 }, {
//-0.949107912342758, 0.949107912342758, 0 }, { -0.949107912342758,
// 0.949107912342758, 0.405845151377397 }, { -0.949107912342758,
// 0.949107912342758, 0.741531185599394 }, { -0.949107912342758,
// 0.949107912342758, 0.949107912342758 }, { -0.741531185599394,
//-0.949107912342758, -0.949107912342758 }, { -0.741531185599394,
//-0.949107912342758, -0.741531185599394 }, { -0.741531185599394,
//-0.949107912342758, -0.405845151377397 }, { -0.741531185599394,
//-0.949107912342758, 0 }, { -0.741531185599394, -0.949107912342758,
// 0.405845151377397 }, { -0.741531185599394, -0.949107912342758,
// 0.741531185599394 }, { -0.741531185599394, -0.949107912342758,
// 0.949107912342758 }, { -0.741531185599394, -0.741531185599394,
//-0.949107912342758 }, { -0.741531185599394, -0.741531185599394,
//-0.741531185599394 }, { -0.741531185599394, -0.741531185599394,
//-0.405845151377397 }, { -0.741531185599394, -0.741531185599394, 0 }, {
//-0.741531185599394, -0.741531185599394, 0.405845151377397 }, {
// -0.741531185599394, -0.741531185599394, 0.741531185599394 }, {
// -0.741531185599394, -0.741531185599394, 0.949107912342758 }, {
// -0.741531185599394, -0.405845151377397, -0.949107912342758 }, {
//-0.741531185599394, -0.405845151377397, -0.741531185599394 }, {
//-0.741531185599394, -0.405845151377397, -0.405845151377397 }, {
//-0.741531185599394, -0.405845151377397, 0 }, { -0.741531185599394,
//-0.405845151377397, 0.405845151377397 }, { -0.741531185599394,
//-0.405845151377397, 0.741531185599394 }, { -0.741531185599394,
//-0.405845151377397, 0.949107912342758 }, { -0.741531185599394, 0,
//-0.949107912342758 }, { -0.741531185599394, 0, -0.741531185599394 }, {
//-0.741531185599394, 0, -0.405845151377397 }, { -0.741531185599394, 0, 0 }, {
//-0.741531185599394, 0, 0.405845151377397 }, { -0.741531185599394, 0,
// 0.741531185599394 }, { -0.741531185599394, 0, 0.949107912342758 }, {
// -0.741531185599394, 0.405845151377397, -0.949107912342758 }, {
//-0.741531185599394, 0.405845151377397, -0.741531185599394 }, {
//-0.741531185599394, 0.405845151377397, -0.405845151377397 }, {
//-0.741531185599394, 0.405845151377397, 0 }, { -0.741531185599394,
// 0.405845151377397, 0.405845151377397 }, { -0.741531185599394,
// 0.405845151377397, 0.741531185599394 }, { -0.741531185599394,
// 0.405845151377397, 0.949107912342758 }, { -0.741531185599394,
// 0.741531185599394, -0.949107912342758 }, { -0.741531185599394,
// 0.741531185599394, -0.741531185599394 }, { -0.741531185599394,
// 0.741531185599394, -0.405845151377397 }, { -0.741531185599394,
// 0.741531185599394, 0 }, { -0.741531185599394, 0.741531185599394,
// 0.405845151377397 }, { -0.741531185599394, 0.741531185599394,
// 0.741531185599394 }, { -0.741531185599394, 0.741531185599394,
// 0.949107912342758 }, { -0.741531185599394, 0.949107912342758,
// -0.949107912342758 }, { -0.741531185599394, 0.949107912342758,
// -0.741531185599394 }, { -0.741531185599394, 0.949107912342758,
// -0.405845151377397 }, { -0.741531185599394, 0.949107912342758, 0 }, {
// -0.741531185599394, 0.949107912342758, 0.405845151377397 }, {
// -0.741531185599394, 0.949107912342758, 0.741531185599394 }, {
// -0.741531185599394, 0.949107912342758, 0.949107912342758 }, {
// -0.405845151377397, -0.949107912342758, -0.949107912342758 }, {
//-0.405845151377397, -0.949107912342758, -0.741531185599394 }, {
//-0.405845151377397, -0.949107912342758, -0.405845151377397 }, {
//-0.405845151377397, -0.949107912342758, 0 }, { -0.405845151377397,
//-0.949107912342758, 0.405845151377397 }, { -0.405845151377397,
//-0.949107912342758, 0.741531185599394 }, { -0.405845151377397,
//-0.949107912342758, 0.949107912342758 }, { -0.405845151377397,
//-0.741531185599394, -0.949107912342758 }, { -0.405845151377397,
//-0.741531185599394, -0.741531185599394 }, { -0.405845151377397,
//-0.741531185599394, -0.405845151377397 }, { -0.405845151377397,
//-0.741531185599394, 0 }, { -0.405845151377397, -0.741531185599394,
// 0.405845151377397 }, { -0.405845151377397, -0.741531185599394,
// 0.741531185599394 }, { -0.405845151377397, -0.741531185599394,
// 0.949107912342758 }, { -0.405845151377397, -0.405845151377397,
//-0.949107912342758 }, { -0.405845151377397, -0.405845151377397,
//-0.741531185599394 }, { -0.405845151377397, -0.405845151377397,
//-0.405845151377397 }, { -0.405845151377397, -0.405845151377397, 0 }, {
//-0.405845151377397, -0.405845151377397, 0.405845151377397 }, {
// -0.405845151377397, -0.405845151377397, 0.741531185599394 }, {
// -0.405845151377397, -0.405845151377397, 0.949107912342758 }, {
// -0.405845151377397, 0, -0.949107912342758 }, { -0.405845151377397, 0,
//-0.741531185599394 }, { -0.405845151377397, 0, -0.405845151377397 }, {
//-0.405845151377397, 0, 0 }, { -0.405845151377397, 0, 0.405845151377397 }, {
// -0.405845151377397, 0, 0.741531185599394 }, { -0.405845151377397, 0,
// 0.949107912342758 }, { -0.405845151377397, 0.405845151377397,
// -0.949107912342758 }, { -0.405845151377397, 0.405845151377397,
// -0.741531185599394 }, { -0.405845151377397, 0.405845151377397,
// -0.405845151377397 }, { -0.405845151377397, 0.405845151377397, 0 }, {
// -0.405845151377397, 0.405845151377397, 0.405845151377397 }, {
// -0.405845151377397, 0.405845151377397, 0.741531185599394 }, {
// -0.405845151377397, 0.405845151377397, 0.949107912342758 }, {
// -0.405845151377397, 0.741531185599394, -0.949107912342758 }, {
//-0.405845151377397, 0.741531185599394, -0.741531185599394 }, {
//-0.405845151377397, 0.741531185599394, -0.405845151377397 }, {
//-0.405845151377397, 0.741531185599394, 0 }, { -0.405845151377397,
// 0.741531185599394, 0.405845151377397 }, { -0.405845151377397,
// 0.741531185599394, 0.741531185599394 }, { -0.405845151377397,
// 0.741531185599394, 0.949107912342758 }, { -0.405845151377397,
// 0.949107912342758, -0.949107912342758 }, { -0.405845151377397,
// 0.949107912342758, -0.741531185599394 }, { -0.405845151377397,
// 0.949107912342758, -0.405845151377397 }, { -0.405845151377397,
// 0.949107912342758, 0 }, { -0.405845151377397, 0.949107912342758,
// 0.405845151377397 }, { -0.405845151377397, 0.949107912342758,
// 0.741531185599394 }, { -0.405845151377397, 0.949107912342758,
// 0.949107912342758 }, { 0, -0.949107912342758, -0.949107912342758 }, { 0,
//-0.949107912342758, -0.741531185599394 }, { 0, -0.949107912342758,
//-0.405845151377397 }, { 0, -0.949107912342758, 0 }, { 0, -0.949107912342758,
// 0.405845151377397 }, { 0, -0.949107912342758, 0.741531185599394 }, { 0,
// -0.949107912342758, 0.949107912342758 }, { 0, -0.741531185599394,
//-0.949107912342758 }, { 0, -0.741531185599394, -0.741531185599394 }, { 0,
//-0.741531185599394, -0.405845151377397 }, { 0, -0.741531185599394, 0 }, { 0,
//-0.741531185599394, 0.405845151377397 }, { 0, -0.741531185599394,
// 0.741531185599394 }, { 0, -0.741531185599394, 0.949107912342758 }, { 0,
// -0.405845151377397, -0.949107912342758 }, { 0, -0.405845151377397,
//-0.741531185599394 }, { 0, -0.405845151377397, -0.405845151377397 }, { 0,
//-0.405845151377397, 0 }, { 0, -0.405845151377397, 0.405845151377397 }, { 0,
// -0.405845151377397, 0.741531185599394 }, { 0, -0.405845151377397,
// 0.949107912342758 }, { 0, 0, -0.949107912342758 }, { 0, 0, -0.741531185599394
//}, { 0, 0, -0.405845151377397 }, { 0, 0, 0 }, { 0, 0, 0.405845151377397 }, {
// 0, 0, 0.741531185599394 }, { 0, 0, 0.949107912342758 }, { 0,
// 0.405845151377397, -0.949107912342758 }, { 0, 0.405845151377397,
//-0.741531185599394 }, { 0, 0.405845151377397, -0.405845151377397 }, { 0,
// 0.405845151377397, 0 }, { 0, 0.405845151377397, 0.405845151377397 }, { 0,
// 0.405845151377397, 0.741531185599394 }, { 0, 0.405845151377397,
// 0.949107912342758 }, { 0, 0.741531185599394, -0.949107912342758 }, { 0,
// 0.741531185599394, -0.741531185599394 }, { 0, 0.741531185599394,
//-0.405845151377397 }, { 0, 0.741531185599394, 0 }, { 0, 0.741531185599394,
// 0.405845151377397 }, { 0, 0.741531185599394, 0.741531185599394 }, { 0,
// 0.741531185599394, 0.949107912342758 }, { 0, 0.949107912342758,
//-0.949107912342758 }, { 0, 0.949107912342758, -0.741531185599394 }, { 0,
// 0.949107912342758, -0.405845151377397 }, { 0, 0.949107912342758, 0 }, { 0,
// 0.949107912342758, 0.405845151377397 }, { 0, 0.949107912342758,
// 0.741531185599394 }, { 0, 0.949107912342758, 0.949107912342758 }, {
// 0.405845151377397, -0.949107912342758, -0.949107912342758 }, {
// 0.405845151377397, -0.949107912342758, -0.741531185599394 }, {
// 0.405845151377397, -0.949107912342758, -0.405845151377397 }, {
// 0.405845151377397, -0.949107912342758, 0 }, { 0.405845151377397,
//-0.949107912342758, 0.405845151377397 }, { 0.405845151377397,
//-0.949107912342758, 0.741531185599394 }, { 0.405845151377397,
//-0.949107912342758, 0.949107912342758 }, { 0.405845151377397,
//-0.741531185599394, -0.949107912342758 }, { 0.405845151377397,
//-0.741531185599394, -0.741531185599394 }, { 0.405845151377397,
//-0.741531185599394, -0.405845151377397 }, { 0.405845151377397,
//-0.741531185599394, 0 }, { 0.405845151377397, -0.741531185599394,
// 0.405845151377397 }, { 0.405845151377397, -0.741531185599394,
// 0.741531185599394 }, { 0.405845151377397, -0.741531185599394,
// 0.949107912342758 }, { 0.405845151377397, -0.405845151377397,
//-0.949107912342758 }, { 0.405845151377397, -0.405845151377397,
//-0.741531185599394 }, { 0.405845151377397, -0.405845151377397,
//-0.405845151377397 }, { 0.405845151377397, -0.405845151377397, 0 }, {
// 0.405845151377397, -0.405845151377397, 0.405845151377397 }, {
// 0.405845151377397, -0.405845151377397, 0.741531185599394 }, {
// 0.405845151377397, -0.405845151377397, 0.949107912342758 }, {
// 0.405845151377397, 0, -0.949107912342758 }, { 0.405845151377397, 0,
//-0.741531185599394 }, { 0.405845151377397, 0, -0.405845151377397 }, {
// 0.405845151377397, 0, 0 }, { 0.405845151377397, 0, 0.405845151377397 }, {
// 0.405845151377397, 0, 0.741531185599394 }, { 0.405845151377397, 0,
// 0.949107912342758 }, { 0.405845151377397, 0.405845151377397,
//-0.949107912342758 }, { 0.405845151377397, 0.405845151377397,
//-0.741531185599394 }, { 0.405845151377397, 0.405845151377397,
//-0.405845151377397 }, { 0.405845151377397, 0.405845151377397, 0 }, {
// 0.405845151377397, 0.405845151377397, 0.405845151377397 }, {
// 0.405845151377397, 0.405845151377397, 0.741531185599394 }, {
// 0.405845151377397, 0.405845151377397, 0.949107912342758 }, {
// 0.405845151377397, 0.741531185599394, -0.949107912342758 }, {
// 0.405845151377397, 0.741531185599394, -0.741531185599394 }, {
// 0.405845151377397, 0.741531185599394, -0.405845151377397 }, {
// 0.405845151377397, 0.741531185599394, 0 }, { 0.405845151377397,
// 0.741531185599394, 0.405845151377397 }, { 0.405845151377397,
// 0.741531185599394, 0.741531185599394 }, { 0.405845151377397,
// 0.741531185599394, 0.949107912342758 }, { 0.405845151377397,
// 0.949107912342758, -0.949107912342758 }, { 0.405845151377397,
// 0.949107912342758, -0.741531185599394 }, { 0.405845151377397,
// 0.949107912342758, -0.405845151377397 }, { 0.405845151377397,
// 0.949107912342758, 0 }, { 0.405845151377397, 0.949107912342758,
// 0.405845151377397 }, { 0.405845151377397, 0.949107912342758,
// 0.741531185599394
//}, { 0.405845151377397, 0.949107912342758, 0.949107912342758 }, {
// 0.741531185599394, -0.949107912342758, -0.949107912342758 }, {
// 0.741531185599394, -0.949107912342758, -0.741531185599394 }, {
// 0.741531185599394, -0.949107912342758, -0.405845151377397 }, {
// 0.741531185599394, -0.949107912342758, 0 }, { 0.741531185599394,
//-0.949107912342758, 0.405845151377397 }, { 0.741531185599394,
//-0.949107912342758, 0.741531185599394 }, { 0.741531185599394,
//-0.949107912342758, 0.949107912342758 }, { 0.741531185599394,
//-0.741531185599394, -0.949107912342758 }, { 0.741531185599394,
//-0.741531185599394, -0.741531185599394 }, { 0.741531185599394,
//-0.741531185599394, -0.405845151377397 }, { 0.741531185599394,
//-0.741531185599394, 0 }, { 0.741531185599394, -0.741531185599394,
// 0.405845151377397 }, { 0.741531185599394, -0.741531185599394,
// 0.741531185599394 }, { 0.741531185599394, -0.741531185599394,
// 0.949107912342758 }, { 0.741531185599394, -0.405845151377397,
//-0.949107912342758 }, { 0.741531185599394, -0.405845151377397,
//-0.741531185599394 }, { 0.741531185599394, -0.405845151377397,
//-0.405845151377397 }, { 0.741531185599394, -0.405845151377397, 0 }, {
// 0.741531185599394, -0.405845151377397, 0.405845151377397 }, {
// 0.741531185599394, -0.405845151377397, 0.741531185599394 }, {
// 0.741531185599394, -0.405845151377397, 0.949107912342758 }, {
// 0.741531185599394, 0, -0.949107912342758 }, { 0.741531185599394, 0,
//-0.741531185599394 }, { 0.741531185599394, 0, -0.405845151377397 }, {
// 0.741531185599394, 0, 0 }, { 0.741531185599394, 0, 0.405845151377397 }, {
// 0.741531185599394, 0, 0.741531185599394 }, { 0.741531185599394, 0,
// 0.949107912342758 }, { 0.741531185599394, 0.405845151377397,
//-0.949107912342758 }, { 0.741531185599394, 0.405845151377397,
//-0.741531185599394 }, { 0.741531185599394, 0.405845151377397,
//-0.405845151377397 }, { 0.741531185599394, 0.405845151377397, 0 }, {
// 0.741531185599394, 0.405845151377397, 0.405845151377397 }, {
// 0.741531185599394, 0.405845151377397, 0.741531185599394 }, {
// 0.741531185599394, 0.405845151377397, 0.949107912342758 }, {
// 0.741531185599394, 0.741531185599394, -0.949107912342758 }, {
// 0.741531185599394, 0.741531185599394, -0.741531185599394 }, {
// 0.741531185599394, 0.741531185599394, -0.405845151377397 }, {
// 0.741531185599394, 0.741531185599394, 0 }, { 0.741531185599394,
// 0.741531185599394, 0.405845151377397 }, { 0.741531185599394,
// 0.741531185599394, 0.741531185599394 }, { 0.741531185599394,
// 0.741531185599394, 0.949107912342758 }, { 0.741531185599394,
// 0.949107912342758, -0.949107912342758 }, { 0.741531185599394,
// 0.949107912342758, -0.741531185599394 }, { 0.741531185599394,
// 0.949107912342758, -0.405845151377397 }, { 0.741531185599394,
// 0.949107912342758, 0 }, { 0.741531185599394, 0.949107912342758,
// 0.405845151377397 }, { 0.741531185599394, 0.949107912342758,
// 0.741531185599394
//}, { 0.741531185599394, 0.949107912342758, 0.949107912342758 }, {
// 0.949107912342758, -0.949107912342758, -0.949107912342758 }, {
// 0.949107912342758, -0.949107912342758, -0.741531185599394 }, {
// 0.949107912342758, -0.949107912342758, -0.405845151377397 }, {
// 0.949107912342758, -0.949107912342758, 0 }, { 0.949107912342758,
//-0.949107912342758, 0.405845151377397 }, { 0.949107912342758,
//-0.949107912342758, 0.741531185599394 }, { 0.949107912342758,
//-0.949107912342758, 0.949107912342758 }, { 0.949107912342758,
//-0.741531185599394, -0.949107912342758 }, { 0.949107912342758,
//-0.741531185599394, -0.741531185599394 }, { 0.949107912342758,
//-0.741531185599394, -0.405845151377397 }, { 0.949107912342758,
//-0.741531185599394, 0 }, { 0.949107912342758, -0.741531185599394,
// 0.405845151377397 }, { 0.949107912342758, -0.741531185599394,
// 0.741531185599394 }, { 0.949107912342758, -0.741531185599394,
// 0.949107912342758 }, { 0.949107912342758, -0.405845151377397,
//-0.949107912342758 }, { 0.949107912342758, -0.405845151377397,
//-0.741531185599394 }, { 0.949107912342758, -0.405845151377397,
//-0.405845151377397 }, { 0.949107912342758, -0.405845151377397, 0 }, {
// 0.949107912342758, -0.405845151377397, 0.405845151377397 }, {
// 0.949107912342758, -0.405845151377397, 0.741531185599394 }, {
// 0.949107912342758, -0.405845151377397, 0.949107912342758 }, {
// 0.949107912342758, 0, -0.949107912342758 }, { 0.949107912342758, 0,
//-0.741531185599394 }, { 0.949107912342758, 0, -0.405845151377397 }, {
// 0.949107912342758, 0, 0 }, { 0.949107912342758, 0, 0.405845151377397 }, {
// 0.949107912342758, 0, 0.741531185599394 }, { 0.949107912342758, 0,
// 0.949107912342758 }, { 0.949107912342758, 0.405845151377397,
//-0.949107912342758 }, { 0.949107912342758, 0.405845151377397,
//-0.741531185599394 }, { 0.949107912342758, 0.405845151377397,
//-0.405845151377397 }, { 0.949107912342758, 0.405845151377397, 0 }, {
// 0.949107912342758, 0.405845151377397, 0.405845151377397 }, {
// 0.949107912342758, 0.405845151377397, 0.741531185599394 }, {
// 0.949107912342758, 0.405845151377397, 0.949107912342758 }, {
// 0.949107912342758, 0.741531185599394, -0.949107912342758 }, {
// 0.949107912342758, 0.741531185599394, -0.741531185599394 }, {
// 0.949107912342758, 0.741531185599394, -0.405845151377397 }, {
// 0.949107912342758, 0.741531185599394, 0 }, { 0.949107912342758,
// 0.741531185599394, 0.405845151377397 }, { 0.949107912342758,
// 0.741531185599394, 0.741531185599394 }, { 0.949107912342758,
// 0.741531185599394, 0.949107912342758 }, { 0.949107912342758,
// 0.949107912342758, -0.949107912342758 }, { 0.949107912342758,
// 0.949107912342758, -0.741531185599394 }, { 0.949107912342758,
// 0.949107912342758, -0.405845151377397 }, { 0.949107912342758,
// 0.949107912342758, 0 }, { 0.949107912342758, 0.949107912342758,
// 0.405845151377397 }, { 0.949107912342758, 0.949107912342758,
// 0.741531185599394
//}, { 0.949107912342758, 0.949107912342758, 0.949107912342758 } } , {
// 0.002170991099484325, 0.004689640298542917, 0.006401898735341792,
// 0.007007652660768766, 0.006401898735341792, 0.004689640298542917,
// 0.002170991099484325, 0.004689640298542918, 0.01013027005727551,
// 0.01382898451475967, 0.01513749656732298, 0.01382898451475967,
// 0.01013027005727551, 0.004689640298542918, 0.006401898735341793,
// 0.01382898451475967, 0.01887815543200791, 0.02066442497960781,
// 0.01887815543200791, 0.01382898451475967, 0.006401898735341793,
// 0.007007652660768767, 0.01513749656732297, 0.02066442497960781,
// 0.02261971310045627, 0.02066442497960781, 0.01513749656732297,
// 0.007007652660768767, 0.006401898735341793, 0.01382898451475967,
// 0.01887815543200791, 0.02066442497960781, 0.01887815543200791,
// 0.01382898451475967, 0.006401898735341793, 0.004689640298542918,
// 0.01013027005727551, 0.01382898451475967, 0.01513749656732298,
// 0.01382898451475967, 0.01013027005727551, 0.004689640298542918,
// 0.002170991099484325, 0.004689640298542917, 0.006401898735341792,
// 0.007007652660768766, 0.006401898735341792, 0.004689640298542917,
// 0.002170991099484325, 0.004689640298542918, 0.01013027005727551,
// 0.01382898451475967, 0.01513749656732298, 0.01382898451475967,
// 0.01013027005727551, 0.004689640298542918, 0.01013027005727551,
// 0.0218827809598143, 0.02987251448600937, 0.0326990810501411,
// 0.02987251448600937, 0.0218827809598143, 0.01013027005727551,
// 0.01382898451475967, 0.02987251448600937, 0.04077942028280544,
// 0.04463800895066054, 0.04077942028280544, 0.02987251448600937,
// 0.01382898451475967, 0.01513749656732297, 0.0326990810501411,
// 0.04463800895066054, 0.04886170105560343, 0.04463800895066054,
// 0.0326990810501411, 0.01513749656732297, 0.01382898451475967,
// 0.02987251448600937, 0.04077942028280544, 0.04463800895066054,
// 0.04077942028280544, 0.02987251448600937, 0.01382898451475967,
// 0.01013027005727551, 0.0218827809598143, 0.02987251448600937,
// 0.0326990810501411, 0.02987251448600937, 0.0218827809598143,
// 0.01013027005727551, 0.004689640298542918, 0.01013027005727551,
// 0.01382898451475967, 0.01513749656732298, 0.01382898451475967,
// 0.01013027005727551, 0.004689640298542918, 0.006401898735341793,
// 0.01382898451475967, 0.01887815543200791, 0.02066442497960781,
// 0.01887815543200791, 0.01382898451475967, 0.006401898735341793,
// 0.01382898451475967, 0.02987251448600937, 0.04077942028280544,
// 0.04463800895066054, 0.04077942028280544, 0.02987251448600937,
// 0.01382898451475967, 0.01887815543200791, 0.04077942028280544,
// 0.05566860196454243, 0.06093601957877199, 0.05566860196454243,
// 0.04077942028280544, 0.01887815543200791, 0.02066442497960781,
// 0.04463800895066054, 0.06093601957877198, 0.06670184540415741,
// 0.06093601957877198, 0.04463800895066054, 0.02066442497960781,
// 0.01887815543200791, 0.04077942028280544, 0.05566860196454243,
// 0.06093601957877199, 0.05566860196454243, 0.04077942028280544,
// 0.01887815543200791, 0.01382898451475967, 0.02987251448600937,
// 0.04077942028280544, 0.04463800895066054, 0.04077942028280544,
// 0.02987251448600937, 0.01382898451475967, 0.006401898735341793,
// 0.01382898451475967, 0.01887815543200791, 0.02066442497960781,
// 0.01887815543200791, 0.01382898451475967, 0.006401898735341793,
// 0.007007652660768767, 0.01513749656732297, 0.02066442497960781,
// 0.02261971310045627, 0.02066442497960781, 0.01513749656732297,
// 0.007007652660768767, 0.01513749656732297, 0.0326990810501411,
// 0.04463800895066054, 0.04886170105560343, 0.04463800895066054,
// 0.0326990810501411, 0.01513749656732297, 0.02066442497960781,
// 0.04463800895066054, 0.06093601957877198, 0.06670184540415741,
// 0.06093601957877198, 0.04463800895066054, 0.02066442497960781,
// 0.02261971310045627, 0.04886170105560343, 0.06670184540415743,
// 0.07301323931355114, 0.06670184540415743, 0.04886170105560343,
// 0.02261971310045627, 0.02066442497960781, 0.04463800895066054,
// 0.06093601957877198, 0.06670184540415741, 0.06093601957877198,
// 0.04463800895066054, 0.02066442497960781, 0.01513749656732297,
// 0.0326990810501411, 0.04463800895066054, 0.04886170105560343,
// 0.04463800895066054, 0.0326990810501411, 0.01513749656732297,
// 0.007007652660768767, 0.01513749656732297, 0.02066442497960781,
// 0.02261971310045627, 0.02066442497960781, 0.01513749656732297,
// 0.007007652660768767, 0.006401898735341793, 0.01382898451475967,
// 0.01887815543200791, 0.02066442497960781, 0.01887815543200791,
// 0.01382898451475967, 0.006401898735341793, 0.01382898451475967,
// 0.02987251448600937, 0.04077942028280544, 0.04463800895066054,
// 0.04077942028280544, 0.02987251448600937, 0.01382898451475967,
// 0.01887815543200791, 0.04077942028280544, 0.05566860196454243,
// 0.06093601957877199, 0.05566860196454243, 0.04077942028280544,
// 0.01887815543200791, 0.02066442497960781, 0.04463800895066054,
// 0.06093601957877198, 0.06670184540415741, 0.06093601957877198,
// 0.04463800895066054, 0.02066442497960781, 0.01887815543200791,
// 0.04077942028280544, 0.05566860196454243, 0.06093601957877199,
// 0.05566860196454243, 0.04077942028280544, 0.01887815543200791,
// 0.01382898451475967, 0.02987251448600937, 0.04077942028280544,
// 0.04463800895066054, 0.04077942028280544, 0.02987251448600937,
// 0.01382898451475967, 0.006401898735341793, 0.01382898451475967,
// 0.01887815543200791, 0.02066442497960781, 0.01887815543200791,
// 0.01382898451475967, 0.006401898735341793, 0.004689640298542918,
// 0.01013027005727551, 0.01382898451475967, 0.01513749656732298,
// 0.01382898451475967, 0.01013027005727551, 0.004689640298542918,
// 0.01013027005727551, 0.0218827809598143, 0.02987251448600937,
// 0.0326990810501411, 0.02987251448600937, 0.0218827809598143,
// 0.01013027005727551, 0.01382898451475967, 0.02987251448600937,
// 0.04077942028280544, 0.04463800895066054, 0.04077942028280544,
// 0.02987251448600937, 0.01382898451475967, 0.01513749656732297,
// 0.0326990810501411, 0.04463800895066054, 0.04886170105560343,
// 0.04463800895066054, 0.0326990810501411, 0.01513749656732297,
// 0.01382898451475967, 0.02987251448600937, 0.04077942028280544,
// 0.04463800895066054, 0.04077942028280544, 0.02987251448600937,
// 0.01382898451475967, 0.01013027005727551, 0.0218827809598143,
// 0.02987251448600937, 0.0326990810501411, 0.02987251448600937,
// 0.0218827809598143, 0.01013027005727551, 0.004689640298542918,
// 0.01013027005727551, 0.01382898451475967, 0.01513749656732298,
// 0.01382898451475967, 0.01013027005727551, 0.004689640298542918,
// 0.002170991099484325, 0.004689640298542917, 0.006401898735341792,
// 0.007007652660768766, 0.006401898735341792, 0.004689640298542917,
// 0.002170991099484325, 0.004689640298542918, 0.01013027005727551,
// 0.01382898451475967, 0.01513749656732298, 0.01382898451475967,
// 0.01013027005727551, 0.004689640298542918, 0.006401898735341793,
// 0.01382898451475967, 0.01887815543200791, 0.02066442497960781,
// 0.01887815543200791, 0.01382898451475967, 0.006401898735341793,
// 0.007007652660768767, 0.01513749656732297, 0.02066442497960781,
// 0.02261971310045627, 0.02066442497960781, 0.01513749656732297,
// 0.007007652660768767, 0.006401898735341793, 0.01382898451475967,
// 0.01887815543200791, 0.02066442497960781, 0.01887815543200791,
// 0.01382898451475967, 0.006401898735341793, 0.004689640298542918,
// 0.01013027005727551, 0.01382898451475967, 0.01513749656732298,
// 0.01382898451475967, 0.01013027005727551, 0.004689640298542918,
// 0.002170991099484325, 0.004689640298542917, 0.006401898735341792,
// 0.007007652660768766, 0.006401898735341792, 0.004689640298542917,
// 0.002170991099484325 } }; 					case 14:
// case 15:	return { { { -0.960289856497536, -0.960289856497536,
//-0.960289856497536 }, { -0.960289856497536, -0.960289856497536,
//-0.796666477413627 }, { -0.960289856497536, -0.960289856497536,
//-0.525532409916329 }, { -0.960289856497536, -0.960289856497536,
//-0.18343464249565 }, { -0.960289856497536, -0.960289856497536,
// 0.18343464249565 }, { -0.960289856497536, -0.960289856497536,
// 0.525532409916329 }, { -0.960289856497536, -0.960289856497536,
// 0.796666477413627 }, { -0.960289856497536, -0.960289856497536,
// 0.960289856497536 }, { -0.960289856497536, -0.796666477413627,
//-0.960289856497536 }, { -0.960289856497536, -0.796666477413627,
//-0.796666477413627 }, { -0.960289856497536, -0.796666477413627,
//-0.525532409916329 }, { -0.960289856497536, -0.796666477413627,
//-0.18343464249565 }, { -0.960289856497536, -0.796666477413627,
// 0.18343464249565 }, { -0.960289856497536, -0.796666477413627,
// 0.525532409916329 }, { -0.960289856497536, -0.796666477413627,
// 0.796666477413627 }, { -0.960289856497536, -0.796666477413627,
// 0.960289856497536 }, { -0.960289856497536, -0.525532409916329,
//-0.960289856497536 }, { -0.960289856497536, -0.525532409916329,
//-0.796666477413627 }, { -0.960289856497536, -0.525532409916329,
//-0.525532409916329 }, { -0.960289856497536, -0.525532409916329,
//-0.18343464249565 }, { -0.960289856497536, -0.525532409916329,
// 0.18343464249565 }, { -0.960289856497536, -0.525532409916329,
// 0.525532409916329 }, { -0.960289856497536, -0.525532409916329,
// 0.796666477413627 }, { -0.960289856497536, -0.525532409916329,
// 0.960289856497536 }, { -0.960289856497536, -0.18343464249565,
//-0.960289856497536 }, { -0.960289856497536, -0.18343464249565,
//-0.796666477413627 }, { -0.960289856497536, -0.18343464249565,
//-0.525532409916329 }, { -0.960289856497536, -0.18343464249565,
//-0.18343464249565 }, { -0.960289856497536, -0.18343464249565, 0.18343464249565
//}, { -0.960289856497536, -0.18343464249565, 0.525532409916329 }, {
//-0.960289856497536, -0.18343464249565, 0.796666477413627 }, {
//-0.960289856497536, -0.18343464249565, 0.960289856497536 }, {
//-0.960289856497536, 0.18343464249565, -0.960289856497536 }, {
// -0.960289856497536, 0.18343464249565, -0.796666477413627 }, {
// -0.960289856497536, 0.18343464249565, -0.525532409916329 }, {
// -0.960289856497536, 0.18343464249565, -0.18343464249565 }, {
// -0.960289856497536, 0.18343464249565, 0.18343464249565 }, {
// -0.960289856497536, 0.18343464249565, 0.525532409916329
//}, { -0.960289856497536, 0.18343464249565, 0.796666477413627 }, {
//-0.960289856497536, 0.18343464249565, 0.960289856497536 }, {
//-0.960289856497536, 0.525532409916329, -0.960289856497536 }, {
//-0.960289856497536, 0.525532409916329, -0.796666477413627 }, {
//-0.960289856497536, 0.525532409916329, -0.525532409916329 }, {
//-0.960289856497536, 0.525532409916329, -0.18343464249565 }, {
//-0.960289856497536, 0.525532409916329, 0.18343464249565 }, {
//-0.960289856497536, 0.525532409916329, 0.525532409916329 }, {
//-0.960289856497536, 0.525532409916329, 0.796666477413627 }, {
//-0.960289856497536, 0.525532409916329, 0.960289856497536 }, {
//-0.960289856497536, 0.796666477413627, -0.960289856497536 }, {
//-0.960289856497536, 0.796666477413627, -0.796666477413627 }, {
//-0.960289856497536, 0.796666477413627, -0.525532409916329 }, {
//-0.960289856497536, 0.796666477413627, -0.18343464249565 }, {
//-0.960289856497536, 0.796666477413627, 0.18343464249565 }, {
//-0.960289856497536, 0.796666477413627, 0.525532409916329 }, {
//-0.960289856497536, 0.796666477413627, 0.796666477413627 }, {
//-0.960289856497536, 0.796666477413627, 0.960289856497536 }, {
//-0.960289856497536, 0.960289856497536, -0.960289856497536 }, {
//-0.960289856497536, 0.960289856497536, -0.796666477413627 }, {
//-0.960289856497536, 0.960289856497536, -0.525532409916329 }, {
//-0.960289856497536, 0.960289856497536, -0.18343464249565 }, {
//-0.960289856497536, 0.960289856497536, 0.18343464249565 }, {
//-0.960289856497536, 0.960289856497536, 0.525532409916329 }, {
//-0.960289856497536, 0.960289856497536, 0.796666477413627 }, {
//-0.960289856497536, 0.960289856497536, 0.960289856497536 }, {
//-0.796666477413627, -0.960289856497536, -0.960289856497536 }, {
//-0.796666477413627, -0.960289856497536, -0.796666477413627 }, {
//-0.796666477413627, -0.960289856497536, -0.525532409916329 }, {
//-0.796666477413627, -0.960289856497536, -0.18343464249565 }, {
//-0.796666477413627, -0.960289856497536, 0.18343464249565 }, {
//-0.796666477413627, -0.960289856497536, 0.525532409916329 }, {
//-0.796666477413627, -0.960289856497536, 0.796666477413627 }, {
//-0.796666477413627, -0.960289856497536, 0.960289856497536 }, {
//-0.796666477413627, -0.796666477413627, -0.960289856497536 }, {
//-0.796666477413627, -0.796666477413627, -0.796666477413627 }, {
//-0.796666477413627, -0.796666477413627, -0.525532409916329 }, {
//-0.796666477413627, -0.796666477413627, -0.18343464249565 }, {
//-0.796666477413627, -0.796666477413627, 0.18343464249565 }, {
//-0.796666477413627, -0.796666477413627, 0.525532409916329 }, {
//-0.796666477413627, -0.796666477413627, 0.796666477413627 }, {
//-0.796666477413627, -0.796666477413627, 0.960289856497536 }, {
//-0.796666477413627, -0.525532409916329, -0.960289856497536 }, {
//-0.796666477413627, -0.525532409916329, -0.796666477413627 }, {
//-0.796666477413627, -0.525532409916329, -0.525532409916329 }, {
//-0.796666477413627, -0.525532409916329, -0.18343464249565 }, {
//-0.796666477413627, -0.525532409916329, 0.18343464249565 }, {
//-0.796666477413627, -0.525532409916329, 0.525532409916329 }, {
//-0.796666477413627, -0.525532409916329, 0.796666477413627 }, {
//-0.796666477413627, -0.525532409916329, 0.960289856497536 }, {
//-0.796666477413627, -0.18343464249565, -0.960289856497536 }, {
//-0.796666477413627, -0.18343464249565, -0.796666477413627 }, {
//-0.796666477413627, -0.18343464249565, -0.525532409916329 }, {
//-0.796666477413627, -0.18343464249565, -0.18343464249565 }, {
//-0.796666477413627, -0.18343464249565, 0.18343464249565 }, {
//-0.796666477413627, -0.18343464249565, 0.525532409916329 }, {
//-0.796666477413627, -0.18343464249565, 0.796666477413627 }, {
//-0.796666477413627, -0.18343464249565, 0.960289856497536 }, {
//-0.796666477413627, 0.18343464249565, -0.960289856497536 }, {
//-0.796666477413627, 0.18343464249565, -0.796666477413627 }, {
//-0.796666477413627, 0.18343464249565, -0.525532409916329 }, {
//-0.796666477413627, 0.18343464249565, -0.18343464249565 }, {
//-0.796666477413627, 0.18343464249565, 0.18343464249565 }, {
//-0.796666477413627, 0.18343464249565, 0.525532409916329 }, {
//-0.796666477413627, 0.18343464249565, 0.796666477413627 }, {
//-0.796666477413627, 0.18343464249565, 0.960289856497536 }, {
//-0.796666477413627, 0.525532409916329, -0.960289856497536 }, {
//-0.796666477413627, 0.525532409916329, -0.796666477413627 }, {
//-0.796666477413627, 0.525532409916329, -0.525532409916329 }, {
//-0.796666477413627, 0.525532409916329, -0.18343464249565 }, {
//-0.796666477413627, 0.525532409916329, 0.18343464249565 }, {
//-0.796666477413627, 0.525532409916329, 0.525532409916329 }, {
//-0.796666477413627, 0.525532409916329, 0.796666477413627 }, {
//-0.796666477413627, 0.525532409916329, 0.960289856497536 }, {
//-0.796666477413627, 0.796666477413627, -0.960289856497536 }, {
//-0.796666477413627, 0.796666477413627, -0.796666477413627 }, {
//-0.796666477413627, 0.796666477413627, -0.525532409916329 }, {
//-0.796666477413627, 0.796666477413627, -0.18343464249565 }, {
//-0.796666477413627, 0.796666477413627, 0.18343464249565 }, {
//-0.796666477413627, 0.796666477413627, 0.525532409916329 }, {
//-0.796666477413627, 0.796666477413627, 0.796666477413627 }, {
//-0.796666477413627, 0.796666477413627, 0.960289856497536 }, {
//-0.796666477413627, 0.960289856497536, -0.960289856497536 }, {
//-0.796666477413627, 0.960289856497536, -0.796666477413627 }, {
//-0.796666477413627, 0.960289856497536, -0.525532409916329 }, {
//-0.796666477413627, 0.960289856497536, -0.18343464249565 }, {
//-0.796666477413627, 0.960289856497536, 0.18343464249565 }, {
//-0.796666477413627, 0.960289856497536, 0.525532409916329 }, {
//-0.796666477413627, 0.960289856497536, 0.796666477413627 }, {
//-0.796666477413627, 0.960289856497536, 0.960289856497536 }, {
//-0.525532409916329, -0.960289856497536, -0.960289856497536 }, {
//-0.525532409916329, -0.960289856497536, -0.796666477413627 }, {
//-0.525532409916329, -0.960289856497536, -0.525532409916329 }, {
//-0.525532409916329, -0.960289856497536, -0.18343464249565 }, {
//-0.525532409916329, -0.960289856497536, 0.18343464249565 }, {
//-0.525532409916329, -0.960289856497536, 0.525532409916329 }, {
//-0.525532409916329, -0.960289856497536, 0.796666477413627 }, {
//-0.525532409916329, -0.960289856497536, 0.960289856497536 }, {
//-0.525532409916329, -0.796666477413627, -0.960289856497536 }, {
//-0.525532409916329, -0.796666477413627, -0.796666477413627 }, {
//-0.525532409916329, -0.796666477413627, -0.525532409916329 }, {
//-0.525532409916329, -0.796666477413627, -0.18343464249565 }, {
//-0.525532409916329, -0.796666477413627, 0.18343464249565 }, {
//-0.525532409916329, -0.796666477413627, 0.525532409916329 }, {
//-0.525532409916329, -0.796666477413627, 0.796666477413627 }, {
//-0.525532409916329, -0.796666477413627, 0.960289856497536 }, {
//-0.525532409916329, -0.525532409916329, -0.960289856497536 }, {
//-0.525532409916329, -0.525532409916329, -0.796666477413627 }, {
//-0.525532409916329, -0.525532409916329, -0.525532409916329 }, {
//-0.525532409916329, -0.525532409916329, -0.18343464249565 }, {
//-0.525532409916329, -0.525532409916329, 0.18343464249565 }, {
//-0.525532409916329, -0.525532409916329, 0.525532409916329 }, {
//-0.525532409916329, -0.525532409916329, 0.796666477413627 }, {
//-0.525532409916329, -0.525532409916329, 0.960289856497536 }, {
//-0.525532409916329, -0.18343464249565, -0.960289856497536 }, {
//-0.525532409916329, -0.18343464249565, -0.796666477413627 }, {
//-0.525532409916329, -0.18343464249565, -0.525532409916329 }, {
//-0.525532409916329, -0.18343464249565, -0.18343464249565 }, {
//-0.525532409916329, -0.18343464249565, 0.18343464249565 }, {
//-0.525532409916329, -0.18343464249565, 0.525532409916329 }, {
//-0.525532409916329, -0.18343464249565, 0.796666477413627 }, {
//-0.525532409916329, -0.18343464249565, 0.960289856497536 }, {
//-0.525532409916329, 0.18343464249565, -0.960289856497536 }, {
//-0.525532409916329, 0.18343464249565, -0.796666477413627 }, {
//-0.525532409916329, 0.18343464249565, -0.525532409916329 }, {
//-0.525532409916329, 0.18343464249565, -0.18343464249565 }, {
//-0.525532409916329, 0.18343464249565, 0.18343464249565 }, {
//-0.525532409916329, 0.18343464249565, 0.525532409916329 }, {
//-0.525532409916329, 0.18343464249565, 0.796666477413627 }, {
//-0.525532409916329, 0.18343464249565, 0.960289856497536 }, {
//-0.525532409916329, 0.525532409916329, -0.960289856497536 }, {
//-0.525532409916329, 0.525532409916329, -0.796666477413627 }, {
//-0.525532409916329, 0.525532409916329, -0.525532409916329 }, {
//-0.525532409916329, 0.525532409916329, -0.18343464249565 }, {
//-0.525532409916329, 0.525532409916329, 0.18343464249565 }, {
//-0.525532409916329, 0.525532409916329, 0.525532409916329 }, {
//-0.525532409916329, 0.525532409916329, 0.796666477413627 }, {
//-0.525532409916329, 0.525532409916329, 0.960289856497536 }, {
//-0.525532409916329, 0.796666477413627, -0.960289856497536 }, {
//-0.525532409916329, 0.796666477413627, -0.796666477413627 }, {
//-0.525532409916329, 0.796666477413627, -0.525532409916329 }, {
//-0.525532409916329, 0.796666477413627, -0.18343464249565 }, {
//-0.525532409916329, 0.796666477413627, 0.18343464249565 }, {
//-0.525532409916329, 0.796666477413627, 0.525532409916329 }, {
//-0.525532409916329, 0.796666477413627, 0.796666477413627 }, {
//-0.525532409916329, 0.796666477413627, 0.960289856497536 }, {
//-0.525532409916329, 0.960289856497536, -0.960289856497536 }, {
//-0.525532409916329, 0.960289856497536, -0.796666477413627 }, {
//-0.525532409916329, 0.960289856497536, -0.525532409916329 }, {
//-0.525532409916329, 0.960289856497536, -0.18343464249565 }, {
//-0.525532409916329, 0.960289856497536, 0.18343464249565 }, {
//-0.525532409916329, 0.960289856497536, 0.525532409916329 }, {
//-0.525532409916329, 0.960289856497536, 0.796666477413627 }, {
//-0.525532409916329, 0.960289856497536, 0.960289856497536 }, {
//-0.18343464249565, -0.960289856497536, -0.960289856497536 }, {
//-0.18343464249565, -0.960289856497536, -0.796666477413627 }, {
//-0.18343464249565, -0.960289856497536, -0.525532409916329 }, {
//-0.18343464249565, -0.960289856497536, -0.18343464249565 }, {
//-0.18343464249565, -0.960289856497536, 0.18343464249565 }, {
//-0.18343464249565, -0.960289856497536, 0.525532409916329 }, {
//-0.18343464249565, -0.960289856497536, 0.796666477413627 }, {
//-0.18343464249565, -0.960289856497536, 0.960289856497536 }, {
//-0.18343464249565, -0.796666477413627, -0.960289856497536 }, {
//-0.18343464249565, -0.796666477413627, -0.796666477413627 }, {
//-0.18343464249565, -0.796666477413627, -0.525532409916329 }, {
//-0.18343464249565, -0.796666477413627, -0.18343464249565 }, {
//-0.18343464249565, -0.796666477413627, 0.18343464249565 }, {
//-0.18343464249565, -0.796666477413627, 0.525532409916329 }, {
//-0.18343464249565, -0.796666477413627, 0.796666477413627 }, {
//-0.18343464249565, -0.796666477413627, 0.960289856497536 }, {
//-0.18343464249565, -0.525532409916329, -0.960289856497536 }, {
//-0.18343464249565, -0.525532409916329, -0.796666477413627 }, {
//-0.18343464249565, -0.525532409916329, -0.525532409916329 }, {
//-0.18343464249565, -0.525532409916329, -0.18343464249565 }, {
//-0.18343464249565, -0.525532409916329, 0.18343464249565 }, {
//-0.18343464249565, -0.525532409916329, 0.525532409916329 }, {
//-0.18343464249565, -0.525532409916329, 0.796666477413627 }, {
//-0.18343464249565, -0.525532409916329, 0.960289856497536 }, {
//-0.18343464249565, -0.18343464249565, -0.960289856497536 }, {
//-0.18343464249565, -0.18343464249565, -0.796666477413627 }, {
//-0.18343464249565, -0.18343464249565, -0.525532409916329 }, {
//-0.18343464249565, -0.18343464249565, -0.18343464249565 }, {
//-0.18343464249565, -0.18343464249565, 0.18343464249565 }, { -0.18343464249565,
//-0.18343464249565, 0.525532409916329 }, { -0.18343464249565,
//-0.18343464249565, 0.796666477413627 }, { -0.18343464249565,
//-0.18343464249565, 0.960289856497536 }, { -0.18343464249565, 0.18343464249565,
//-0.960289856497536 }, { -0.18343464249565, 0.18343464249565,
//-0.796666477413627 }, { -0.18343464249565, 0.18343464249565,
//-0.525532409916329 }, { -0.18343464249565, 0.18343464249565, -0.18343464249565
//}, { -0.18343464249565, 0.18343464249565, 0.18343464249565 }, {
//-0.18343464249565, 0.18343464249565, 0.525532409916329 }, { -0.18343464249565,
// 0.18343464249565, 0.796666477413627 }, { -0.18343464249565, 0.18343464249565,
// 0.960289856497536 }, { -0.18343464249565, 0.525532409916329,
//-0.960289856497536 }, { -0.18343464249565, 0.525532409916329,
//-0.796666477413627 }, { -0.18343464249565, 0.525532409916329,
//-0.525532409916329 }, { -0.18343464249565, 0.525532409916329,
//-0.18343464249565 }, { -0.18343464249565, 0.525532409916329, 0.18343464249565
//}, { -0.18343464249565, 0.525532409916329, 0.525532409916329 }, {
//-0.18343464249565, 0.525532409916329, 0.796666477413627 }, {
//-0.18343464249565, 0.525532409916329, 0.960289856497536 }, {
//-0.18343464249565, 0.796666477413627, -0.960289856497536 }, {
//-0.18343464249565, 0.796666477413627, -0.796666477413627 }, {
//-0.18343464249565, 0.796666477413627, -0.525532409916329 }, {
//-0.18343464249565, 0.796666477413627, -0.18343464249565 }, {
//-0.18343464249565, 0.796666477413627, 0.18343464249565 }, { -0.18343464249565,
// 0.796666477413627, 0.525532409916329 }, { -0.18343464249565,
// 0.796666477413627, 0.796666477413627 }, { -0.18343464249565,
// 0.796666477413627, 0.960289856497536 }, { -0.18343464249565,
// 0.960289856497536, -0.960289856497536 }, { -0.18343464249565,
// 0.960289856497536, -0.796666477413627 }, { -0.18343464249565,
// 0.960289856497536, -0.525532409916329 }, { -0.18343464249565,
// 0.960289856497536, -0.18343464249565 }, { -0.18343464249565,
// 0.960289856497536, 0.18343464249565 }, { -0.18343464249565,
// 0.960289856497536, 0.525532409916329 }, { -0.18343464249565,
// 0.960289856497536, 0.796666477413627
//}, { -0.18343464249565, 0.960289856497536, 0.960289856497536 }, {
// 0.18343464249565, -0.960289856497536, -0.960289856497536 }, {
// 0.18343464249565, -0.960289856497536, -0.796666477413627 }, {
// 0.18343464249565, -0.960289856497536, -0.525532409916329 }, {
// 0.18343464249565, -0.960289856497536, -0.18343464249565 }, {
// 0.18343464249565, -0.960289856497536, 0.18343464249565 }, { 0.18343464249565,
//-0.960289856497536, 0.525532409916329 }, { 0.18343464249565,
//-0.960289856497536, 0.796666477413627 }, { 0.18343464249565,
//-0.960289856497536, 0.960289856497536 }, { 0.18343464249565,
//-0.796666477413627, -0.960289856497536 }, { 0.18343464249565,
//-0.796666477413627, -0.796666477413627 }, { 0.18343464249565,
//-0.796666477413627, -0.525532409916329 }, { 0.18343464249565,
//-0.796666477413627, -0.18343464249565 }, { 0.18343464249565,
//-0.796666477413627, 0.18343464249565 }, { 0.18343464249565,
//-0.796666477413627, 0.525532409916329 }, { 0.18343464249565,
//-0.796666477413627, 0.796666477413627 }, { 0.18343464249565,
//-0.796666477413627, 0.960289856497536 }, { 0.18343464249565,
//-0.525532409916329, -0.960289856497536 }, { 0.18343464249565,
//-0.525532409916329, -0.796666477413627 }, { 0.18343464249565,
//-0.525532409916329, -0.525532409916329 }, { 0.18343464249565,
//-0.525532409916329, -0.18343464249565 }, { 0.18343464249565,
//-0.525532409916329, 0.18343464249565 }, { 0.18343464249565,
//-0.525532409916329, 0.525532409916329 }, { 0.18343464249565,
//-0.525532409916329, 0.796666477413627 }, { 0.18343464249565,
//-0.525532409916329, 0.960289856497536 }, { 0.18343464249565,
//-0.18343464249565, -0.960289856497536 }, { 0.18343464249565,
//-0.18343464249565, -0.796666477413627 }, { 0.18343464249565,
//-0.18343464249565, -0.525532409916329 }, { 0.18343464249565,
//-0.18343464249565, -0.18343464249565 }, { 0.18343464249565, -0.18343464249565,
// 0.18343464249565 }, { 0.18343464249565, -0.18343464249565, 0.525532409916329
//}, { 0.18343464249565, -0.18343464249565, 0.796666477413627 }, {
// 0.18343464249565, -0.18343464249565, 0.960289856497536 }, { 0.18343464249565,
// 0.18343464249565, -0.960289856497536 }, { 0.18343464249565, 0.18343464249565,
//-0.796666477413627 }, { 0.18343464249565, 0.18343464249565, -0.525532409916329
//}, { 0.18343464249565, 0.18343464249565, -0.18343464249565 }, {
// 0.18343464249565, 0.18343464249565, 0.18343464249565 }, { 0.18343464249565,
// 0.18343464249565, 0.525532409916329 }, { 0.18343464249565, 0.18343464249565,
// 0.796666477413627 }, { 0.18343464249565, 0.18343464249565, 0.960289856497536
//}, { 0.18343464249565, 0.525532409916329, -0.960289856497536 }, {
// 0.18343464249565, 0.525532409916329, -0.796666477413627 }, {
// 0.18343464249565, 0.525532409916329, -0.525532409916329 }, {
// 0.18343464249565, 0.525532409916329, -0.18343464249565 }, { 0.18343464249565,
// 0.525532409916329, 0.18343464249565 }, { 0.18343464249565, 0.525532409916329,
// 0.525532409916329
//}, { 0.18343464249565, 0.525532409916329, 0.796666477413627 }, {
// 0.18343464249565, 0.525532409916329, 0.960289856497536 }, { 0.18343464249565,
// 0.796666477413627, -0.960289856497536 }, { 0.18343464249565,
// 0.796666477413627, -0.796666477413627 }, { 0.18343464249565,
// 0.796666477413627, -0.525532409916329 }, { 0.18343464249565,
// 0.796666477413627, -0.18343464249565 }, { 0.18343464249565,
// 0.796666477413627, 0.18343464249565 }, { 0.18343464249565, 0.796666477413627,
// 0.525532409916329
//}, { 0.18343464249565, 0.796666477413627, 0.796666477413627 }, {
// 0.18343464249565, 0.796666477413627, 0.960289856497536 }, { 0.18343464249565,
// 0.960289856497536, -0.960289856497536 }, { 0.18343464249565,
// 0.960289856497536, -0.796666477413627 }, { 0.18343464249565,
// 0.960289856497536, -0.525532409916329 }, { 0.18343464249565,
// 0.960289856497536, -0.18343464249565 }, { 0.18343464249565,
// 0.960289856497536, 0.18343464249565 }, { 0.18343464249565, 0.960289856497536,
// 0.525532409916329
//}, { 0.18343464249565, 0.960289856497536, 0.796666477413627 }, {
// 0.18343464249565, 0.960289856497536, 0.960289856497536 }, {
// 0.525532409916329, -0.960289856497536, -0.960289856497536 }, {
// 0.525532409916329, -0.960289856497536, -0.796666477413627 }, {
// 0.525532409916329, -0.960289856497536, -0.525532409916329 }, {
// 0.525532409916329, -0.960289856497536, -0.18343464249565 }, {
// 0.525532409916329, -0.960289856497536, 0.18343464249565 }, {
// 0.525532409916329, -0.960289856497536, 0.525532409916329 }, {
// 0.525532409916329, -0.960289856497536, 0.796666477413627 }, {
// 0.525532409916329, -0.960289856497536, 0.960289856497536 }, {
// 0.525532409916329, -0.796666477413627, -0.960289856497536 }, {
// 0.525532409916329, -0.796666477413627, -0.796666477413627 }, {
// 0.525532409916329, -0.796666477413627, -0.525532409916329 }, {
// 0.525532409916329, -0.796666477413627, -0.18343464249565 }, {
// 0.525532409916329, -0.796666477413627, 0.18343464249565 }, {
// 0.525532409916329, -0.796666477413627, 0.525532409916329 }, {
// 0.525532409916329, -0.796666477413627, 0.796666477413627 }, {
// 0.525532409916329, -0.796666477413627, 0.960289856497536 }, {
// 0.525532409916329, -0.525532409916329, -0.960289856497536 }, {
// 0.525532409916329, -0.525532409916329, -0.796666477413627 }, {
// 0.525532409916329, -0.525532409916329, -0.525532409916329 }, {
// 0.525532409916329, -0.525532409916329, -0.18343464249565 }, {
// 0.525532409916329, -0.525532409916329, 0.18343464249565 }, {
// 0.525532409916329, -0.525532409916329, 0.525532409916329 }, {
// 0.525532409916329, -0.525532409916329, 0.796666477413627 }, {
// 0.525532409916329, -0.525532409916329, 0.960289856497536 }, {
// 0.525532409916329, -0.18343464249565, -0.960289856497536 }, {
// 0.525532409916329, -0.18343464249565, -0.796666477413627 }, {
// 0.525532409916329, -0.18343464249565, -0.525532409916329 }, {
// 0.525532409916329, -0.18343464249565, -0.18343464249565 }, {
// 0.525532409916329, -0.18343464249565, 0.18343464249565 }, {
// 0.525532409916329, -0.18343464249565, 0.525532409916329 }, {
// 0.525532409916329,
// -0.18343464249565, 0.796666477413627
//}, { 0.525532409916329, -0.18343464249565, 0.960289856497536 }, {
// 0.525532409916329, 0.18343464249565, -0.960289856497536 }, {
// 0.525532409916329, 0.18343464249565, -0.796666477413627 }, {
// 0.525532409916329, 0.18343464249565, -0.525532409916329 }, {
// 0.525532409916329, 0.18343464249565, -0.18343464249565 }, {
// 0.525532409916329, 0.18343464249565, 0.18343464249565 }, { 0.525532409916329,
// 0.18343464249565, 0.525532409916329 }, { 0.525532409916329, 0.18343464249565,
// 0.796666477413627
//}, { 0.525532409916329, 0.18343464249565, 0.960289856497536 }, {
// 0.525532409916329, 0.525532409916329, -0.960289856497536 }, {
// 0.525532409916329, 0.525532409916329, -0.796666477413627 }, {
// 0.525532409916329, 0.525532409916329, -0.525532409916329 }, {
// 0.525532409916329, 0.525532409916329, -0.18343464249565 }, {
// 0.525532409916329, 0.525532409916329, 0.18343464249565 }, {
// 0.525532409916329, 0.525532409916329, 0.525532409916329 }, {
// 0.525532409916329, 0.525532409916329, 0.796666477413627 }, {
// 0.525532409916329, 0.525532409916329, 0.960289856497536 }, {
// 0.525532409916329, 0.796666477413627, -0.960289856497536 }, {
// 0.525532409916329, 0.796666477413627, -0.796666477413627 }, {
// 0.525532409916329, 0.796666477413627, -0.525532409916329 }, {
// 0.525532409916329, 0.796666477413627, -0.18343464249565 }, {
// 0.525532409916329, 0.796666477413627, 0.18343464249565 }, {
// 0.525532409916329, 0.796666477413627, 0.525532409916329 }, {
// 0.525532409916329, 0.796666477413627, 0.796666477413627
//}, { 0.525532409916329, 0.796666477413627, 0.960289856497536 }, {
// 0.525532409916329, 0.960289856497536, -0.960289856497536 }, {
// 0.525532409916329, 0.960289856497536, -0.796666477413627 }, {
// 0.525532409916329, 0.960289856497536, -0.525532409916329 }, {
// 0.525532409916329, 0.960289856497536, -0.18343464249565 }, {
// 0.525532409916329, 0.960289856497536, 0.18343464249565 }, {
// 0.525532409916329, 0.960289856497536, 0.525532409916329 }, {
// 0.525532409916329, 0.960289856497536, 0.796666477413627 }, {
// 0.525532409916329, 0.960289856497536, 0.960289856497536 }, {
// 0.796666477413627, -0.960289856497536, -0.960289856497536 }, {
// 0.796666477413627, -0.960289856497536, -0.796666477413627 }, {
// 0.796666477413627, -0.960289856497536, -0.525532409916329 }, {
// 0.796666477413627, -0.960289856497536, -0.18343464249565 }, {
// 0.796666477413627, -0.960289856497536, 0.18343464249565 }, {
// 0.796666477413627, -0.960289856497536, 0.525532409916329 }, {
// 0.796666477413627, -0.960289856497536, 0.796666477413627 }, {
// 0.796666477413627, -0.960289856497536, 0.960289856497536 }, {
// 0.796666477413627, -0.796666477413627, -0.960289856497536 }, {
// 0.796666477413627, -0.796666477413627, -0.796666477413627 }, {
// 0.796666477413627, -0.796666477413627, -0.525532409916329 }, {
// 0.796666477413627, -0.796666477413627, -0.18343464249565 }, {
// 0.796666477413627, -0.796666477413627, 0.18343464249565 }, {
// 0.796666477413627, -0.796666477413627, 0.525532409916329 }, {
// 0.796666477413627, -0.796666477413627, 0.796666477413627 }, {
// 0.796666477413627, -0.796666477413627, 0.960289856497536 }, {
// 0.796666477413627, -0.525532409916329, -0.960289856497536 }, {
// 0.796666477413627, -0.525532409916329, -0.796666477413627 }, {
// 0.796666477413627, -0.525532409916329, -0.525532409916329 }, {
// 0.796666477413627, -0.525532409916329, -0.18343464249565 }, {
// 0.796666477413627, -0.525532409916329, 0.18343464249565 }, {
// 0.796666477413627, -0.525532409916329, 0.525532409916329 }, {
// 0.796666477413627, -0.525532409916329, 0.796666477413627 }, {
// 0.796666477413627, -0.525532409916329, 0.960289856497536 }, {
// 0.796666477413627, -0.18343464249565, -0.960289856497536 }, {
// 0.796666477413627, -0.18343464249565, -0.796666477413627 }, {
// 0.796666477413627, -0.18343464249565, -0.525532409916329 }, {
// 0.796666477413627, -0.18343464249565, -0.18343464249565 }, {
// 0.796666477413627, -0.18343464249565, 0.18343464249565 }, {
// 0.796666477413627, -0.18343464249565, 0.525532409916329 }, {
// 0.796666477413627,
// -0.18343464249565, 0.796666477413627
//}, { 0.796666477413627, -0.18343464249565, 0.960289856497536 }, {
// 0.796666477413627, 0.18343464249565, -0.960289856497536 }, {
// 0.796666477413627, 0.18343464249565, -0.796666477413627 }, {
// 0.796666477413627, 0.18343464249565, -0.525532409916329 }, {
// 0.796666477413627, 0.18343464249565, -0.18343464249565 }, {
// 0.796666477413627, 0.18343464249565, 0.18343464249565 }, { 0.796666477413627,
// 0.18343464249565, 0.525532409916329 }, { 0.796666477413627, 0.18343464249565,
// 0.796666477413627
//}, { 0.796666477413627, 0.18343464249565, 0.960289856497536 }, {
// 0.796666477413627, 0.525532409916329, -0.960289856497536 }, {
// 0.796666477413627, 0.525532409916329, -0.796666477413627 }, {
// 0.796666477413627, 0.525532409916329, -0.525532409916329 }, {
// 0.796666477413627, 0.525532409916329, -0.18343464249565 }, {
// 0.796666477413627, 0.525532409916329, 0.18343464249565 }, {
// 0.796666477413627, 0.525532409916329, 0.525532409916329 }, {
// 0.796666477413627, 0.525532409916329, 0.796666477413627 }, {
// 0.796666477413627, 0.525532409916329, 0.960289856497536 }, {
// 0.796666477413627, 0.796666477413627, -0.960289856497536 }, {
// 0.796666477413627, 0.796666477413627, -0.796666477413627 }, {
// 0.796666477413627, 0.796666477413627, -0.525532409916329 }, {
// 0.796666477413627, 0.796666477413627, -0.18343464249565 }, {
// 0.796666477413627, 0.796666477413627, 0.18343464249565 }, {
// 0.796666477413627, 0.796666477413627, 0.525532409916329 }, {
// 0.796666477413627, 0.796666477413627, 0.796666477413627
//}, { 0.796666477413627, 0.796666477413627, 0.960289856497536 }, {
// 0.796666477413627, 0.960289856497536, -0.960289856497536 }, {
// 0.796666477413627, 0.960289856497536, -0.796666477413627 }, {
// 0.796666477413627, 0.960289856497536, -0.525532409916329 }, {
// 0.796666477413627, 0.960289856497536, -0.18343464249565 }, {
// 0.796666477413627, 0.960289856497536, 0.18343464249565 }, {
// 0.796666477413627, 0.960289856497536, 0.525532409916329 }, {
// 0.796666477413627, 0.960289856497536, 0.796666477413627 }, {
// 0.796666477413627, 0.960289856497536, 0.960289856497536 }, {
// 0.960289856497536, -0.960289856497536, -0.960289856497536 }, {
// 0.960289856497536, -0.960289856497536, -0.796666477413627 }, {
// 0.960289856497536, -0.960289856497536, -0.525532409916329 }, {
// 0.960289856497536, -0.960289856497536, -0.18343464249565 }, {
// 0.960289856497536, -0.960289856497536, 0.18343464249565 }, {
// 0.960289856497536, -0.960289856497536, 0.525532409916329 }, {
// 0.960289856497536, -0.960289856497536, 0.796666477413627 }, {
// 0.960289856497536, -0.960289856497536, 0.960289856497536 }, {
// 0.960289856497536, -0.796666477413627, -0.960289856497536 }, {
// 0.960289856497536, -0.796666477413627, -0.796666477413627 }, {
// 0.960289856497536, -0.796666477413627, -0.525532409916329 }, {
// 0.960289856497536, -0.796666477413627, -0.18343464249565 }, {
// 0.960289856497536, -0.796666477413627, 0.18343464249565 }, {
// 0.960289856497536, -0.796666477413627, 0.525532409916329 }, {
// 0.960289856497536, -0.796666477413627, 0.796666477413627 }, {
// 0.960289856497536, -0.796666477413627, 0.960289856497536 }, {
// 0.960289856497536, -0.525532409916329, -0.960289856497536 }, {
// 0.960289856497536, -0.525532409916329, -0.796666477413627 }, {
// 0.960289856497536, -0.525532409916329, -0.525532409916329 }, {
// 0.960289856497536, -0.525532409916329, -0.18343464249565 }, {
// 0.960289856497536, -0.525532409916329, 0.18343464249565 }, {
// 0.960289856497536, -0.525532409916329, 0.525532409916329 }, {
// 0.960289856497536, -0.525532409916329, 0.796666477413627 }, {
// 0.960289856497536, -0.525532409916329, 0.960289856497536 }, {
// 0.960289856497536, -0.18343464249565, -0.960289856497536 }, {
// 0.960289856497536, -0.18343464249565, -0.796666477413627 }, {
// 0.960289856497536, -0.18343464249565, -0.525532409916329 }, {
// 0.960289856497536, -0.18343464249565, -0.18343464249565 }, {
// 0.960289856497536, -0.18343464249565, 0.18343464249565 }, {
// 0.960289856497536, -0.18343464249565, 0.525532409916329 }, {
// 0.960289856497536,
// -0.18343464249565, 0.796666477413627
//}, { 0.960289856497536, -0.18343464249565, 0.960289856497536 }, {
// 0.960289856497536, 0.18343464249565, -0.960289856497536 }, {
// 0.960289856497536, 0.18343464249565, -0.796666477413627 }, {
// 0.960289856497536, 0.18343464249565, -0.525532409916329 }, {
// 0.960289856497536, 0.18343464249565, -0.18343464249565 }, {
// 0.960289856497536, 0.18343464249565, 0.18343464249565 }, { 0.960289856497536,
// 0.18343464249565, 0.525532409916329 }, { 0.960289856497536, 0.18343464249565,
// 0.796666477413627
//}, { 0.960289856497536, 0.18343464249565, 0.960289856497536 }, {
// 0.960289856497536, 0.525532409916329, -0.960289856497536 }, {
// 0.960289856497536, 0.525532409916329, -0.796666477413627 }, {
// 0.960289856497536, 0.525532409916329, -0.525532409916329 }, {
// 0.960289856497536, 0.525532409916329, -0.18343464249565 }, {
// 0.960289856497536, 0.525532409916329, 0.18343464249565 }, {
// 0.960289856497536, 0.525532409916329, 0.525532409916329 }, {
// 0.960289856497536, 0.525532409916329, 0.796666477413627 }, {
// 0.960289856497536, 0.525532409916329, 0.960289856497536 }, {
// 0.960289856497536, 0.796666477413627, -0.960289856497536 }, {
// 0.960289856497536, 0.796666477413627, -0.796666477413627 }, {
// 0.960289856497536, 0.796666477413627, -0.525532409916329 }, {
// 0.960289856497536, 0.796666477413627, -0.18343464249565 }, {
// 0.960289856497536, 0.796666477413627, 0.18343464249565 }, {
// 0.960289856497536, 0.796666477413627, 0.525532409916329 }, {
// 0.960289856497536, 0.796666477413627, 0.796666477413627
//}, { 0.960289856497536, 0.796666477413627, 0.960289856497536 }, {
// 0.960289856497536, 0.960289856497536, -0.960289856497536 }, {
// 0.960289856497536, 0.960289856497536, -0.796666477413627 }, {
// 0.960289856497536, 0.960289856497536, -0.525532409916329 }, {
// 0.960289856497536, 0.960289856497536, -0.18343464249565 }, {
// 0.960289856497536, 0.960289856497536, 0.18343464249565 }, {
// 0.960289856497536, 0.960289856497536, 0.525532409916329 }, {
// 0.960289856497536, 0.960289856497536, 0.796666477413627 }, {
// 0.960289856497536, 0.960289856497536, 0.960289856497536 } } , {
// 0.001037310733367905, 0.002278786618767613, 0.00321461993646263,
// 0.003716499270894022, 0.003716499270894022, 0.00321461993646263,
// 0.002278786618767613, 0.001037310733367905, 0.002278786618767614,
// 0.005006087652264334, 0.007061946492976767, 0.008164485852446232,
// 0.008164485852446232, 0.007061946492976767, 0.005006087652264334,
// 0.002278786618767614, 0.003214619936462631, 0.007061946492976767,
// 0.009962088507800969, 0.01151740964954176, 0.01151740964954176,
// 0.009962088507800969, 0.007061946492976767, 0.003214619936462631,
// 0.003716499270894023, 0.008164485852446234, 0.01151740964954176,
// 0.01331555375476573, 0.01331555375476573, 0.01151740964954176,
// 0.008164485852446234, 0.003716499270894023, 0.003716499270894023,
// 0.008164485852446234, 0.01151740964954176, 0.01331555375476573,
// 0.01331555375476573, 0.01151740964954176, 0.008164485852446234,
// 0.003716499270894023, 0.003214619936462631, 0.007061946492976767,
// 0.009962088507800969, 0.01151740964954176, 0.01151740964954176,
// 0.009962088507800969, 0.007061946492976767, 0.003214619936462631,
// 0.002278786618767614, 0.005006087652264334, 0.007061946492976767,
// 0.008164485852446232, 0.008164485852446232, 0.007061946492976767,
// 0.005006087652264334, 0.002278786618767614, 0.001037310733367905,
// 0.002278786618767613, 0.00321461993646263, 0.003716499270894022,
// 0.003716499270894022, 0.00321461993646263, 0.002278786618767613,
// 0.001037310733367905, 0.002278786618767614, 0.005006087652264334,
// 0.007061946492976767, 0.008164485852446232, 0.008164485852446232,
// 0.007061946492976767, 0.005006087652264334, 0.002278786618767614,
// 0.005006087652264334, 0.0109974814560332, 0.01551383655155982,
// 0.01793591882469536, 0.01793591882469536, 0.01551383655155982,
// 0.0109974814560332, 0.005006087652264334, 0.007061946492976767,
// 0.01551383655155982, 0.02188493115543988, 0.02530169422524632,
// 0.02530169422524632, 0.02188493115543988, 0.01551383655155982,
// 0.007061946492976767, 0.008164485852446232, 0.01793591882469536,
// 0.02530169422524631, 0.02925189602475564, 0.02925189602475564,
// 0.02530169422524631, 0.01793591882469536, 0.008164485852446232,
// 0.008164485852446232, 0.01793591882469536, 0.02530169422524631,
// 0.02925189602475564, 0.02925189602475564, 0.02530169422524631,
// 0.01793591882469536, 0.008164485852446232, 0.007061946492976767,
// 0.01551383655155982, 0.02188493115543988, 0.02530169422524632,
// 0.02530169422524632, 0.02188493115543988, 0.01551383655155982,
// 0.007061946492976767, 0.005006087652264334, 0.0109974814560332,
// 0.01551383655155982, 0.01793591882469536, 0.01793591882469536,
// 0.01551383655155982, 0.0109974814560332, 0.005006087652264334,
// 0.002278786618767614, 0.005006087652264334, 0.007061946492976767,
// 0.008164485852446232, 0.008164485852446232, 0.007061946492976767,
// 0.005006087652264334, 0.002278786618767614, 0.003214619936462631,
// 0.007061946492976767, 0.009962088507800969, 0.01151740964954176,
// 0.01151740964954176, 0.009962088507800969, 0.007061946492976767,
// 0.003214619936462631, 0.007061946492976767, 0.01551383655155982,
// 0.02188493115543988, 0.02530169422524632, 0.02530169422524632,
// 0.02188493115543988, 0.01551383655155982, 0.007061946492976767,
// 0.009962088507800969, 0.02188493115543988, 0.03087245441103915,
// 0.03569238559367399, 0.03569238559367399, 0.03087245441103915,
// 0.02188493115543988, 0.009962088507800969, 0.01151740964954176,
// 0.02530169422524632, 0.03569238559367399, 0.04126482372946601,
// 0.04126482372946601, 0.03569238559367399, 0.02530169422524632,
// 0.01151740964954176, 0.01151740964954176, 0.02530169422524632,
// 0.03569238559367399, 0.04126482372946601, 0.04126482372946601,
// 0.03569238559367399, 0.02530169422524632, 0.01151740964954176,
// 0.009962088507800969, 0.02188493115543988, 0.03087245441103915,
// 0.03569238559367399, 0.03569238559367399, 0.03087245441103915,
// 0.02188493115543988, 0.009962088507800969, 0.007061946492976767,
// 0.01551383655155982, 0.02188493115543988, 0.02530169422524632,
// 0.02530169422524632, 0.02188493115543988, 0.01551383655155982,
// 0.007061946492976767, 0.003214619936462631, 0.007061946492976767,
// 0.009962088507800969, 0.01151740964954176, 0.01151740964954176,
// 0.009962088507800969, 0.007061946492976767, 0.003214619936462631,
// 0.003716499270894023, 0.008164485852446234, 0.01151740964954176,
// 0.01331555375476573, 0.01331555375476573, 0.01151740964954176,
// 0.008164485852446234, 0.003716499270894023, 0.008164485852446232,
// 0.01793591882469536, 0.02530169422524631, 0.02925189602475564,
// 0.02925189602475564, 0.02530169422524631, 0.01793591882469536,
// 0.008164485852446232, 0.01151740964954176, 0.02530169422524632,
// 0.03569238559367399, 0.04126482372946601, 0.04126482372946601,
// 0.03569238559367399, 0.02530169422524632, 0.01151740964954176,
// 0.01331555375476573, 0.02925189602475564, 0.04126482372946601,
// 0.04770725321665521, 0.04770725321665521, 0.04126482372946601,
// 0.02925189602475564, 0.01331555375476573, 0.01331555375476573,
// 0.02925189602475564, 0.04126482372946601, 0.04770725321665521,
// 0.04770725321665521, 0.04126482372946601, 0.02925189602475564,
// 0.01331555375476573, 0.01151740964954176, 0.02530169422524632,
// 0.03569238559367399, 0.04126482372946601, 0.04126482372946601,
// 0.03569238559367399, 0.02530169422524632, 0.01151740964954176,
// 0.008164485852446232, 0.01793591882469536, 0.02530169422524631,
// 0.02925189602475564, 0.02925189602475564, 0.02530169422524631,
// 0.01793591882469536, 0.008164485852446232, 0.003716499270894023,
// 0.008164485852446234, 0.01151740964954176, 0.01331555375476573,
// 0.01331555375476573, 0.01151740964954176, 0.008164485852446234,
// 0.003716499270894023, 0.003716499270894023, 0.008164485852446234,
// 0.01151740964954176, 0.01331555375476573, 0.01331555375476573,
// 0.01151740964954176, 0.008164485852446234, 0.003716499270894023,
// 0.008164485852446232, 0.01793591882469536, 0.02530169422524631,
// 0.02925189602475564, 0.02925189602475564, 0.02530169422524631,
// 0.01793591882469536, 0.008164485852446232, 0.01151740964954176,
// 0.02530169422524632, 0.03569238559367399, 0.04126482372946601,
// 0.04126482372946601, 0.03569238559367399, 0.02530169422524632,
// 0.01151740964954176, 0.01331555375476573, 0.02925189602475564,
// 0.04126482372946601, 0.04770725321665521, 0.04770725321665521,
// 0.04126482372946601, 0.02925189602475564, 0.01331555375476573,
// 0.01331555375476573, 0.02925189602475564, 0.04126482372946601,
// 0.04770725321665521, 0.04770725321665521, 0.04126482372946601,
// 0.02925189602475564, 0.01331555375476573, 0.01151740964954176,
// 0.02530169422524632, 0.03569238559367399, 0.04126482372946601,
// 0.04126482372946601, 0.03569238559367399, 0.02530169422524632,
// 0.01151740964954176, 0.008164485852446232, 0.01793591882469536,
// 0.02530169422524631, 0.02925189602475564, 0.02925189602475564,
// 0.02530169422524631, 0.01793591882469536, 0.008164485852446232,
// 0.003716499270894023, 0.008164485852446234, 0.01151740964954176,
// 0.01331555375476573, 0.01331555375476573, 0.01151740964954176,
// 0.008164485852446234, 0.003716499270894023, 0.003214619936462631,
// 0.007061946492976767, 0.009962088507800969, 0.01151740964954176,
// 0.01151740964954176, 0.009962088507800969, 0.007061946492976767,
// 0.003214619936462631, 0.007061946492976767, 0.01551383655155982,
// 0.02188493115543988, 0.02530169422524632, 0.02530169422524632,
// 0.02188493115543988, 0.01551383655155982, 0.007061946492976767,
// 0.009962088507800969, 0.02188493115543988, 0.03087245441103915,
// 0.03569238559367399, 0.03569238559367399, 0.03087245441103915,
// 0.02188493115543988, 0.009962088507800969, 0.01151740964954176,
// 0.02530169422524632, 0.03569238559367399, 0.04126482372946601,
// 0.04126482372946601, 0.03569238559367399, 0.02530169422524632,
// 0.01151740964954176, 0.01151740964954176, 0.02530169422524632,
// 0.03569238559367399, 0.04126482372946601, 0.04126482372946601,
// 0.03569238559367399, 0.02530169422524632, 0.01151740964954176,
// 0.009962088507800969, 0.02188493115543988, 0.03087245441103915,
// 0.03569238559367399, 0.03569238559367399, 0.03087245441103915,
// 0.02188493115543988, 0.009962088507800969, 0.007061946492976767,
// 0.01551383655155982, 0.02188493115543988, 0.02530169422524632,
// 0.02530169422524632, 0.02188493115543988, 0.01551383655155982,
// 0.007061946492976767, 0.003214619936462631, 0.007061946492976767,
// 0.009962088507800969, 0.01151740964954176, 0.01151740964954176,
// 0.009962088507800969, 0.007061946492976767, 0.003214619936462631,
// 0.002278786618767614, 0.005006087652264334, 0.007061946492976767,
// 0.008164485852446232, 0.008164485852446232, 0.007061946492976767,
// 0.005006087652264334, 0.002278786618767614, 0.005006087652264334,
// 0.0109974814560332, 0.01551383655155982, 0.01793591882469536,
// 0.01793591882469536, 0.01551383655155982, 0.0109974814560332,
// 0.005006087652264334, 0.007061946492976767, 0.01551383655155982,
// 0.02188493115543988, 0.02530169422524632, 0.02530169422524632,
// 0.02188493115543988, 0.01551383655155982, 0.007061946492976767,
// 0.008164485852446232, 0.01793591882469536, 0.02530169422524631,
// 0.02925189602475564, 0.02925189602475564, 0.02530169422524631,
// 0.01793591882469536, 0.008164485852446232, 0.008164485852446232,
// 0.01793591882469536, 0.02530169422524631, 0.02925189602475564,
// 0.02925189602475564, 0.02530169422524631, 0.01793591882469536,
// 0.008164485852446232, 0.007061946492976767, 0.01551383655155982,
// 0.02188493115543988, 0.02530169422524632, 0.02530169422524632,
// 0.02188493115543988, 0.01551383655155982, 0.007061946492976767,
// 0.005006087652264334, 0.0109974814560332, 0.01551383655155982,
// 0.01793591882469536, 0.01793591882469536, 0.01551383655155982,
// 0.0109974814560332, 0.005006087652264334, 0.002278786618767614,
// 0.005006087652264334, 0.007061946492976767, 0.008164485852446232,
// 0.008164485852446232, 0.007061946492976767, 0.005006087652264334,
// 0.002278786618767614, 0.001037310733367905, 0.002278786618767613,
// 0.00321461993646263, 0.003716499270894022, 0.003716499270894022,
// 0.00321461993646263, 0.002278786618767613, 0.001037310733367905,
// 0.002278786618767614, 0.005006087652264334, 0.007061946492976767,
// 0.008164485852446232, 0.008164485852446232, 0.007061946492976767,
// 0.005006087652264334, 0.002278786618767614, 0.003214619936462631,
// 0.007061946492976767, 0.009962088507800969, 0.01151740964954176,
// 0.01151740964954176, 0.009962088507800969, 0.007061946492976767,
// 0.003214619936462631, 0.003716499270894023, 0.008164485852446234,
// 0.01151740964954176, 0.01331555375476573, 0.01331555375476573,
// 0.01151740964954176, 0.008164485852446234, 0.003716499270894023,
// 0.003716499270894023, 0.008164485852446234, 0.01151740964954176,
// 0.01331555375476573, 0.01331555375476573, 0.01151740964954176,
// 0.008164485852446234, 0.003716499270894023, 0.003214619936462631,
// 0.007061946492976767, 0.009962088507800969, 0.01151740964954176,
// 0.01151740964954176, 0.009962088507800969, 0.007061946492976767,
// 0.003214619936462631, 0.002278786618767614, 0.005006087652264334,
// 0.007061946492976767, 0.008164485852446232, 0.008164485852446232,
// 0.007061946492976767, 0.005006087652264334, 0.002278786618767614,
// 0.001037310733367905, 0.002278786618767613, 0.00321461993646263,
// 0.003716499270894022, 0.003716499270894022, 0.00321461993646263,
// 0.002278786618767613, 0.001037310733367905 } };
// case 16: 					case 17:	return { { {
// -0.968160239507626, -0.968160239507626, -0.968160239507626 }, {
//-0.968160239507626, -0.968160239507626, -0.836031107326636 }, {
//-0.968160239507626, -0.968160239507626, -0.61337143270059 }, {
//-0.968160239507626, -0.968160239507626, -0.324253423403809 }, {
//-0.968160239507626, -0.968160239507626, 0 }, { -0.968160239507626,
//-0.968160239507626, 0.324253423403809 }, { -0.968160239507626,
// -0.968160239507626, 0.61337143270059 }, { -0.968160239507626,
// -0.968160239507626, 0.836031107326636 }, { -0.968160239507626,
// -0.968160239507626, 0.968160239507626 }, { -0.968160239507626,
// -0.836031107326636, -0.968160239507626 }, { -0.968160239507626,
//-0.836031107326636, -0.836031107326636 }, { -0.968160239507626,
//-0.836031107326636, -0.61337143270059 }, { -0.968160239507626,
//-0.836031107326636, -0.324253423403809 }, { -0.968160239507626,
//-0.836031107326636, 0 }, { -0.968160239507626, -0.836031107326636,
// 0.324253423403809 }, { -0.968160239507626, -0.836031107326636,
// 0.61337143270059 }, { -0.968160239507626, -0.836031107326636,
// 0.836031107326636 }, { -0.968160239507626, -0.836031107326636,
// 0.968160239507626 }, { -0.968160239507626, -0.61337143270059,
//-0.968160239507626 }, { -0.968160239507626, -0.61337143270059,
//-0.836031107326636 }, { -0.968160239507626, -0.61337143270059,
//-0.61337143270059 }, { -0.968160239507626, -0.61337143270059,
//-0.324253423403809 }, { -0.968160239507626, -0.61337143270059, 0 }, {
//-0.968160239507626, -0.61337143270059, 0.324253423403809 }, {
//-0.968160239507626, -0.61337143270059, 0.61337143270059 }, {
//-0.968160239507626, -0.61337143270059, 0.836031107326636 }, {
//-0.968160239507626, -0.61337143270059, 0.968160239507626 }, {
//-0.968160239507626, -0.324253423403809, -0.968160239507626 }, {
//-0.968160239507626, -0.324253423403809, -0.836031107326636 }, {
//-0.968160239507626, -0.324253423403809, -0.61337143270059 }, {
//-0.968160239507626, -0.324253423403809, -0.324253423403809 }, {
//-0.968160239507626, -0.324253423403809, 0 }, { -0.968160239507626,
//-0.324253423403809, 0.324253423403809 }, { -0.968160239507626,
// -0.324253423403809, 0.61337143270059 }, { -0.968160239507626,
// -0.324253423403809, 0.836031107326636 }, { -0.968160239507626,
// -0.324253423403809, 0.968160239507626 }, { -0.968160239507626, 0,
// -0.968160239507626 }, { -0.968160239507626, 0, -0.836031107326636 }, {
//-0.968160239507626, 0, -0.61337143270059 }, { -0.968160239507626, 0,
//-0.324253423403809 }, { -0.968160239507626, 0, 0 }, { -0.968160239507626, 0,
// 0.324253423403809 }, { -0.968160239507626, 0, 0.61337143270059 }, {
//-0.968160239507626, 0, 0.836031107326636 }, { -0.968160239507626, 0,
// 0.968160239507626 }, { -0.968160239507626, 0.324253423403809,
//-0.968160239507626 }, { -0.968160239507626, 0.324253423403809,
//-0.836031107326636 }, { -0.968160239507626, 0.324253423403809,
//-0.61337143270059 }, { -0.968160239507626, 0.324253423403809,
//-0.324253423403809 }, { -0.968160239507626, 0.324253423403809, 0 }, {
//-0.968160239507626, 0.324253423403809, 0.324253423403809 }, {
// -0.968160239507626, 0.324253423403809, 0.61337143270059 }, {
// -0.968160239507626, 0.324253423403809, 0.836031107326636 }, {
// -0.968160239507626, 0.324253423403809, 0.968160239507626 }, {
// -0.968160239507626, 0.61337143270059, -0.968160239507626 }, {
// -0.968160239507626, 0.61337143270059, -0.836031107326636 }, {
// -0.968160239507626, 0.61337143270059, -0.61337143270059 }, {
// -0.968160239507626, 0.61337143270059, -0.324253423403809 }, {
//-0.968160239507626, 0.61337143270059, 0 }, { -0.968160239507626,
// 0.61337143270059, 0.324253423403809 }, { -0.968160239507626,
// 0.61337143270059, 0.61337143270059 }, { -0.968160239507626, 0.61337143270059,
// 0.836031107326636
//}, { -0.968160239507626, 0.61337143270059, 0.968160239507626 }, {
//-0.968160239507626, 0.836031107326636, -0.968160239507626 }, {
//-0.968160239507626, 0.836031107326636, -0.836031107326636 }, {
//-0.968160239507626, 0.836031107326636, -0.61337143270059 }, {
//-0.968160239507626, 0.836031107326636, -0.324253423403809 }, {
//-0.968160239507626, 0.836031107326636, 0 }, { -0.968160239507626,
// 0.836031107326636, 0.324253423403809 }, { -0.968160239507626,
// 0.836031107326636, 0.61337143270059 }, { -0.968160239507626,
// 0.836031107326636, 0.836031107326636 }, { -0.968160239507626,
// 0.836031107326636, 0.968160239507626 }, { -0.968160239507626,
// 0.968160239507626, -0.968160239507626 }, { -0.968160239507626,
// 0.968160239507626, -0.836031107326636 }, { -0.968160239507626,
// 0.968160239507626, -0.61337143270059 }, { -0.968160239507626,
// 0.968160239507626, -0.324253423403809 }, { -0.968160239507626,
// 0.968160239507626, 0 }, { -0.968160239507626, 0.968160239507626,
// 0.324253423403809 }, { -0.968160239507626, 0.968160239507626,
// 0.61337143270059
//}, { -0.968160239507626, 0.968160239507626, 0.836031107326636 }, {
//-0.968160239507626, 0.968160239507626, 0.968160239507626 }, {
//-0.836031107326636, -0.968160239507626, -0.968160239507626 }, {
//-0.836031107326636, -0.968160239507626, -0.836031107326636 }, {
//-0.836031107326636, -0.968160239507626, -0.61337143270059 }, {
//-0.836031107326636, -0.968160239507626, -0.324253423403809 }, {
//-0.836031107326636, -0.968160239507626, 0 }, { -0.836031107326636,
//-0.968160239507626, 0.324253423403809 }, { -0.836031107326636,
//-0.968160239507626, 0.61337143270059 }, { -0.836031107326636,
//-0.968160239507626, 0.836031107326636 }, { -0.836031107326636,
//-0.968160239507626, 0.968160239507626 }, { -0.836031107326636,
//-0.836031107326636, -0.968160239507626 }, { -0.836031107326636,
//-0.836031107326636, -0.836031107326636 }, { -0.836031107326636,
//-0.836031107326636, -0.61337143270059 }, { -0.836031107326636,
//-0.836031107326636, -0.324253423403809 }, { -0.836031107326636,
//-0.836031107326636, 0 }, { -0.836031107326636, -0.836031107326636,
// 0.324253423403809 }, { -0.836031107326636, -0.836031107326636,
// 0.61337143270059 }, { -0.836031107326636, -0.836031107326636,
// 0.836031107326636 }, { -0.836031107326636, -0.836031107326636,
// 0.968160239507626 }, { -0.836031107326636, -0.61337143270059,
//-0.968160239507626 }, { -0.836031107326636, -0.61337143270059,
//-0.836031107326636 }, { -0.836031107326636, -0.61337143270059,
//-0.61337143270059 }, { -0.836031107326636, -0.61337143270059,
//-0.324253423403809 }, { -0.836031107326636, -0.61337143270059, 0 }, {
//-0.836031107326636, -0.61337143270059, 0.324253423403809 }, {
//-0.836031107326636, -0.61337143270059, 0.61337143270059 }, {
//-0.836031107326636, -0.61337143270059, 0.836031107326636 }, {
//-0.836031107326636, -0.61337143270059, 0.968160239507626 }, {
//-0.836031107326636, -0.324253423403809, -0.968160239507626 }, {
//-0.836031107326636, -0.324253423403809, -0.836031107326636 }, {
//-0.836031107326636, -0.324253423403809, -0.61337143270059 }, {
//-0.836031107326636, -0.324253423403809, -0.324253423403809 }, {
//-0.836031107326636, -0.324253423403809, 0 }, { -0.836031107326636,
//-0.324253423403809, 0.324253423403809 }, { -0.836031107326636,
//-0.324253423403809, 0.61337143270059 }, { -0.836031107326636,
//-0.324253423403809, 0.836031107326636 }, { -0.836031107326636,
//-0.324253423403809, 0.968160239507626 }, { -0.836031107326636, 0,
//-0.968160239507626 }, { -0.836031107326636, 0, -0.836031107326636 }, {
//-0.836031107326636, 0, -0.61337143270059 }, { -0.836031107326636, 0,
//-0.324253423403809 }, { -0.836031107326636, 0, 0 }, { -0.836031107326636, 0,
// 0.324253423403809 }, { -0.836031107326636, 0, 0.61337143270059 }, {
//-0.836031107326636, 0, 0.836031107326636 }, { -0.836031107326636, 0,
// 0.968160239507626 }, { -0.836031107326636, 0.324253423403809,
//-0.968160239507626 }, { -0.836031107326636, 0.324253423403809,
//-0.836031107326636 }, { -0.836031107326636, 0.324253423403809,
//-0.61337143270059 }, { -0.836031107326636, 0.324253423403809,
//-0.324253423403809 }, { -0.836031107326636, 0.324253423403809, 0 }, {
//-0.836031107326636, 0.324253423403809, 0.324253423403809 }, {
//-0.836031107326636, 0.324253423403809, 0.61337143270059 }, {
//-0.836031107326636, 0.324253423403809, 0.836031107326636 }, {
//-0.836031107326636, 0.324253423403809, 0.968160239507626 }, {
//-0.836031107326636, 0.61337143270059, -0.968160239507626 }, {
//-0.836031107326636, 0.61337143270059, -0.836031107326636 }, {
//-0.836031107326636, 0.61337143270059, -0.61337143270059 }, {
//-0.836031107326636, 0.61337143270059, -0.324253423403809 }, {
//-0.836031107326636, 0.61337143270059, 0 }, { -0.836031107326636,
// 0.61337143270059, 0.324253423403809 }, { -0.836031107326636,
// 0.61337143270059, 0.61337143270059 }, { -0.836031107326636, 0.61337143270059,
// 0.836031107326636
//}, { -0.836031107326636, 0.61337143270059, 0.968160239507626 }, {
//-0.836031107326636, 0.836031107326636, -0.968160239507626 }, {
//-0.836031107326636, 0.836031107326636, -0.836031107326636 }, {
//-0.836031107326636, 0.836031107326636, -0.61337143270059 }, {
//-0.836031107326636, 0.836031107326636, -0.324253423403809 }, {
//-0.836031107326636, 0.836031107326636, 0 }, { -0.836031107326636,
// 0.836031107326636, 0.324253423403809 }, { -0.836031107326636,
// 0.836031107326636, 0.61337143270059 }, { -0.836031107326636,
// 0.836031107326636, 0.836031107326636 }, { -0.836031107326636,
// 0.836031107326636, 0.968160239507626 }, { -0.836031107326636,
// 0.968160239507626, -0.968160239507626 }, { -0.836031107326636,
// 0.968160239507626, -0.836031107326636 }, { -0.836031107326636,
// 0.968160239507626, -0.61337143270059 }, { -0.836031107326636,
// 0.968160239507626, -0.324253423403809 }, { -0.836031107326636,
// 0.968160239507626, 0 }, { -0.836031107326636, 0.968160239507626,
// 0.324253423403809 }, { -0.836031107326636, 0.968160239507626,
// 0.61337143270059
//}, { -0.836031107326636, 0.968160239507626, 0.836031107326636 }, {
//-0.836031107326636, 0.968160239507626, 0.968160239507626 }, {
//-0.61337143270059, -0.968160239507626, -0.968160239507626 }, {
//-0.61337143270059, -0.968160239507626, -0.836031107326636 }, {
//-0.61337143270059, -0.968160239507626, -0.61337143270059 }, {
//-0.61337143270059, -0.968160239507626, -0.324253423403809 }, {
//-0.61337143270059, -0.968160239507626, 0 }, { -0.61337143270059,
//-0.968160239507626, 0.324253423403809 }, { -0.61337143270059,
//-0.968160239507626, 0.61337143270059 }, { -0.61337143270059,
//-0.968160239507626, 0.836031107326636 }, { -0.61337143270059,
//-0.968160239507626, 0.968160239507626 }, { -0.61337143270059,
//-0.836031107326636, -0.968160239507626 }, { -0.61337143270059,
//-0.836031107326636, -0.836031107326636 }, { -0.61337143270059,
//-0.836031107326636, -0.61337143270059 }, { -0.61337143270059,
//-0.836031107326636, -0.324253423403809 }, { -0.61337143270059,
//-0.836031107326636, 0 }, { -0.61337143270059, -0.836031107326636,
// 0.324253423403809 }, { -0.61337143270059, -0.836031107326636,
// 0.61337143270059
//}, { -0.61337143270059, -0.836031107326636, 0.836031107326636 }, {
//-0.61337143270059, -0.836031107326636, 0.968160239507626 }, {
//-0.61337143270059, -0.61337143270059, -0.968160239507626 }, {
//-0.61337143270059, -0.61337143270059, -0.836031107326636 }, {
//-0.61337143270059, -0.61337143270059, -0.61337143270059 }, {
//-0.61337143270059, -0.61337143270059, -0.324253423403809 }, {
//-0.61337143270059, -0.61337143270059, 0 }, { -0.61337143270059,
//-0.61337143270059, 0.324253423403809 }, { -0.61337143270059,
//-0.61337143270059, 0.61337143270059 }, { -0.61337143270059, -0.61337143270059,
// 0.836031107326636 }, { -0.61337143270059, -0.61337143270059,
// 0.968160239507626
//}, { -0.61337143270059, -0.324253423403809, -0.968160239507626 }, {
//-0.61337143270059, -0.324253423403809, -0.836031107326636 }, {
//-0.61337143270059, -0.324253423403809, -0.61337143270059 }, {
//-0.61337143270059, -0.324253423403809, -0.324253423403809 }, {
//-0.61337143270059, -0.324253423403809, 0 }, { -0.61337143270059,
//-0.324253423403809, 0.324253423403809 }, { -0.61337143270059,
//-0.324253423403809, 0.61337143270059 }, { -0.61337143270059,
//-0.324253423403809, 0.836031107326636 }, { -0.61337143270059,
//-0.324253423403809, 0.968160239507626 }, { -0.61337143270059, 0,
//-0.968160239507626 }, { -0.61337143270059, 0, -0.836031107326636 }, {
//-0.61337143270059, 0, -0.61337143270059 }, { -0.61337143270059, 0,
//-0.324253423403809 }, { -0.61337143270059, 0, 0 }, { -0.61337143270059, 0,
// 0.324253423403809 }, { -0.61337143270059, 0, 0.61337143270059 }, {
//-0.61337143270059, 0, 0.836031107326636 }, { -0.61337143270059, 0,
// 0.968160239507626 }, { -0.61337143270059, 0.324253423403809,
//-0.968160239507626 }, { -0.61337143270059, 0.324253423403809,
//-0.836031107326636 }, { -0.61337143270059, 0.324253423403809,
//-0.61337143270059 }, { -0.61337143270059, 0.324253423403809,
//-0.324253423403809 }, { -0.61337143270059, 0.324253423403809, 0 }, {
//-0.61337143270059, 0.324253423403809, 0.324253423403809 }, {
//-0.61337143270059, 0.324253423403809, 0.61337143270059 }, { -0.61337143270059,
// 0.324253423403809, 0.836031107326636 }, { -0.61337143270059,
// 0.324253423403809, 0.968160239507626 }, { -0.61337143270059,
// 0.61337143270059, -0.968160239507626 }, { -0.61337143270059,
// 0.61337143270059, -0.836031107326636 }, { -0.61337143270059,
// 0.61337143270059, -0.61337143270059
//}, { -0.61337143270059, 0.61337143270059, -0.324253423403809 }, {
//-0.61337143270059, 0.61337143270059, 0 }, { -0.61337143270059,
// 0.61337143270059, 0.324253423403809 }, { -0.61337143270059, 0.61337143270059,
// 0.61337143270059 }, { -0.61337143270059, 0.61337143270059, 0.836031107326636
//}, { -0.61337143270059, 0.61337143270059, 0.968160239507626 }, {
//-0.61337143270059, 0.836031107326636, -0.968160239507626 }, {
//-0.61337143270059, 0.836031107326636, -0.836031107326636 }, {
//-0.61337143270059, 0.836031107326636, -0.61337143270059 }, {
//-0.61337143270059, 0.836031107326636, -0.324253423403809 }, {
//-0.61337143270059, 0.836031107326636, 0 }, { -0.61337143270059,
// 0.836031107326636, 0.324253423403809 }, { -0.61337143270059,
// 0.836031107326636, 0.61337143270059 }, { -0.61337143270059,
// 0.836031107326636, 0.836031107326636 }, { -0.61337143270059,
// 0.836031107326636, 0.968160239507626
//}, { -0.61337143270059, 0.968160239507626, -0.968160239507626 }, {
//-0.61337143270059, 0.968160239507626, -0.836031107326636 }, {
//-0.61337143270059, 0.968160239507626, -0.61337143270059 }, {
//-0.61337143270059, 0.968160239507626, -0.324253423403809 }, {
//-0.61337143270059, 0.968160239507626, 0 }, { -0.61337143270059,
// 0.968160239507626, 0.324253423403809 }, { -0.61337143270059,
// 0.968160239507626, 0.61337143270059 }, { -0.61337143270059,
// 0.968160239507626, 0.836031107326636 }, { -0.61337143270059,
// 0.968160239507626, 0.968160239507626
//}, { -0.324253423403809, -0.968160239507626, -0.968160239507626 }, {
//-0.324253423403809, -0.968160239507626, -0.836031107326636 }, {
//-0.324253423403809, -0.968160239507626, -0.61337143270059 }, {
//-0.324253423403809, -0.968160239507626, -0.324253423403809 }, {
//-0.324253423403809, -0.968160239507626, 0 }, { -0.324253423403809,
//-0.968160239507626, 0.324253423403809 }, { -0.324253423403809,
//-0.968160239507626, 0.61337143270059 }, { -0.324253423403809,
//-0.968160239507626, 0.836031107326636 }, { -0.324253423403809,
//-0.968160239507626, 0.968160239507626 }, { -0.324253423403809,
//-0.836031107326636, -0.968160239507626 }, { -0.324253423403809,
//-0.836031107326636, -0.836031107326636 }, { -0.324253423403809,
//-0.836031107326636, -0.61337143270059 }, { -0.324253423403809,
//-0.836031107326636, -0.324253423403809 }, { -0.324253423403809,
//-0.836031107326636, 0 }, { -0.324253423403809, -0.836031107326636,
// 0.324253423403809 }, { -0.324253423403809, -0.836031107326636,
// 0.61337143270059 }, { -0.324253423403809, -0.836031107326636,
// 0.836031107326636 }, { -0.324253423403809, -0.836031107326636,
// 0.968160239507626 }, { -0.324253423403809, -0.61337143270059,
//-0.968160239507626 }, { -0.324253423403809, -0.61337143270059,
//-0.836031107326636 }, { -0.324253423403809, -0.61337143270059,
//-0.61337143270059 }, { -0.324253423403809, -0.61337143270059,
//-0.324253423403809 }, { -0.324253423403809, -0.61337143270059, 0 }, {
//-0.324253423403809, -0.61337143270059, 0.324253423403809 }, {
//-0.324253423403809, -0.61337143270059, 0.61337143270059 }, {
//-0.324253423403809, -0.61337143270059, 0.836031107326636 }, {
//-0.324253423403809, -0.61337143270059, 0.968160239507626 }, {
//-0.324253423403809, -0.324253423403809, -0.968160239507626 }, {
//-0.324253423403809, -0.324253423403809, -0.836031107326636 }, {
//-0.324253423403809, -0.324253423403809, -0.61337143270059 }, {
//-0.324253423403809, -0.324253423403809, -0.324253423403809 }, {
//-0.324253423403809, -0.324253423403809, 0 }, { -0.324253423403809,
//-0.324253423403809, 0.324253423403809 }, { -0.324253423403809,
//-0.324253423403809, 0.61337143270059 }, { -0.324253423403809,
//-0.324253423403809, 0.836031107326636 }, { -0.324253423403809,
//-0.324253423403809, 0.968160239507626 }, { -0.324253423403809, 0,
//-0.968160239507626 }, { -0.324253423403809, 0, -0.836031107326636 }, {
//-0.324253423403809, 0, -0.61337143270059 }, { -0.324253423403809, 0,
//-0.324253423403809 }, { -0.324253423403809, 0, 0 }, { -0.324253423403809, 0,
// 0.324253423403809 }, { -0.324253423403809, 0, 0.61337143270059 }, {
//-0.324253423403809, 0, 0.836031107326636 }, { -0.324253423403809, 0,
// 0.968160239507626 }, { -0.324253423403809, 0.324253423403809,
//-0.968160239507626 }, { -0.324253423403809, 0.324253423403809,
//-0.836031107326636 }, { -0.324253423403809, 0.324253423403809,
//-0.61337143270059 }, { -0.324253423403809, 0.324253423403809,
//-0.324253423403809 }, { -0.324253423403809, 0.324253423403809, 0 }, {
//-0.324253423403809, 0.324253423403809, 0.324253423403809 }, {
//-0.324253423403809, 0.324253423403809, 0.61337143270059 }, {
//-0.324253423403809, 0.324253423403809, 0.836031107326636 }, {
//-0.324253423403809, 0.324253423403809, 0.968160239507626 }, {
//-0.324253423403809, 0.61337143270059, -0.968160239507626 }, {
//-0.324253423403809, 0.61337143270059, -0.836031107326636 }, {
//-0.324253423403809, 0.61337143270059, -0.61337143270059 }, {
//-0.324253423403809, 0.61337143270059, -0.324253423403809 }, {
//-0.324253423403809, 0.61337143270059, 0 }, { -0.324253423403809,
// 0.61337143270059, 0.324253423403809 }, { -0.324253423403809,
// 0.61337143270059, 0.61337143270059 }, { -0.324253423403809, 0.61337143270059,
// 0.836031107326636
//}, { -0.324253423403809, 0.61337143270059, 0.968160239507626 }, {
//-0.324253423403809, 0.836031107326636, -0.968160239507626 }, {
//-0.324253423403809, 0.836031107326636, -0.836031107326636 }, {
//-0.324253423403809, 0.836031107326636, -0.61337143270059 }, {
//-0.324253423403809, 0.836031107326636, -0.324253423403809 }, {
//-0.324253423403809, 0.836031107326636, 0 }, { -0.324253423403809,
// 0.836031107326636, 0.324253423403809 }, { -0.324253423403809,
// 0.836031107326636, 0.61337143270059 }, { -0.324253423403809,
// 0.836031107326636, 0.836031107326636 }, { -0.324253423403809,
// 0.836031107326636, 0.968160239507626 }, { -0.324253423403809,
// 0.968160239507626, -0.968160239507626 }, { -0.324253423403809,
// 0.968160239507626, -0.836031107326636 }, { -0.324253423403809,
// 0.968160239507626, -0.61337143270059 }, { -0.324253423403809,
// 0.968160239507626, -0.324253423403809 }, { -0.324253423403809,
// 0.968160239507626, 0 }, { -0.324253423403809, 0.968160239507626,
// 0.324253423403809 }, { -0.324253423403809, 0.968160239507626,
// 0.61337143270059
//}, { -0.324253423403809, 0.968160239507626, 0.836031107326636 }, {
//-0.324253423403809, 0.968160239507626, 0.968160239507626 }, { 0,
//-0.968160239507626, -0.968160239507626 }, { 0, -0.968160239507626,
//-0.836031107326636 }, { 0, -0.968160239507626, -0.61337143270059 }, { 0,
//-0.968160239507626, -0.324253423403809 }, { 0, -0.968160239507626, 0 }, { 0,
//-0.968160239507626, 0.324253423403809 }, { 0, -0.968160239507626,
// 0.61337143270059 }, { 0, -0.968160239507626, 0.836031107326636 }, { 0,
//-0.968160239507626, 0.968160239507626 }, { 0, -0.836031107326636,
//-0.968160239507626 }, { 0, -0.836031107326636, -0.836031107326636 }, { 0,
//-0.836031107326636, -0.61337143270059 }, { 0, -0.836031107326636,
//-0.324253423403809 }, { 0, -0.836031107326636, 0 }, { 0, -0.836031107326636,
// 0.324253423403809 }, { 0, -0.836031107326636, 0.61337143270059 }, { 0,
//-0.836031107326636, 0.836031107326636 }, { 0, -0.836031107326636,
// 0.968160239507626 }, { 0, -0.61337143270059, -0.968160239507626 }, { 0,
//-0.61337143270059, -0.836031107326636 }, { 0, -0.61337143270059,
//-0.61337143270059 }, { 0, -0.61337143270059, -0.324253423403809 }, { 0,
//-0.61337143270059, 0 }, { 0, -0.61337143270059, 0.324253423403809 }, { 0,
//-0.61337143270059, 0.61337143270059 }, { 0, -0.61337143270059,
// 0.836031107326636 }, { 0, -0.61337143270059, 0.968160239507626 }, { 0,
//-0.324253423403809, -0.968160239507626 }, { 0, -0.324253423403809,
//-0.836031107326636 }, { 0, -0.324253423403809, -0.61337143270059 }, { 0,
//-0.324253423403809, -0.324253423403809 }, { 0, -0.324253423403809, 0 }, { 0,
//-0.324253423403809, 0.324253423403809 }, { 0, -0.324253423403809,
// 0.61337143270059 }, { 0, -0.324253423403809, 0.836031107326636 }, { 0,
//-0.324253423403809, 0.968160239507626 }, { 0, 0, -0.968160239507626 }, { 0, 0,
//-0.836031107326636 }, { 0, 0, -0.61337143270059 }, { 0, 0, -0.324253423403809
//}, { 0, 0, 0 }, { 0, 0, 0.324253423403809 }, { 0, 0, 0.61337143270059 }, { 0,
// 0, 0.836031107326636 }, { 0, 0, 0.968160239507626 }, { 0, 0.324253423403809,
//-0.968160239507626 }, { 0, 0.324253423403809, -0.836031107326636 }, { 0,
// 0.324253423403809, -0.61337143270059 }, { 0, 0.324253423403809,
//-0.324253423403809 }, { 0, 0.324253423403809, 0 }, { 0, 0.324253423403809,
// 0.324253423403809 }, { 0, 0.324253423403809, 0.61337143270059 }, { 0,
// 0.324253423403809, 0.836031107326636 }, { 0, 0.324253423403809,
// 0.968160239507626 }, { 0, 0.61337143270059, -0.968160239507626 }, { 0,
// 0.61337143270059, -0.836031107326636 }, { 0, 0.61337143270059,
//-0.61337143270059 }, { 0, 0.61337143270059, -0.324253423403809 }, { 0,
// 0.61337143270059, 0 }, { 0, 0.61337143270059, 0.324253423403809 }, { 0,
// 0.61337143270059, 0.61337143270059 }, { 0, 0.61337143270059,
// 0.836031107326636
//}, { 0, 0.61337143270059, 0.968160239507626 }, { 0, 0.836031107326636,
//-0.968160239507626 }, { 0, 0.836031107326636, -0.836031107326636 }, { 0,
// 0.836031107326636, -0.61337143270059 }, { 0, 0.836031107326636,
//-0.324253423403809 }, { 0, 0.836031107326636, 0 }, { 0, 0.836031107326636,
// 0.324253423403809 }, { 0, 0.836031107326636, 0.61337143270059 }, { 0,
// 0.836031107326636, 0.836031107326636 }, { 0, 0.836031107326636,
// 0.968160239507626 }, { 0, 0.968160239507626, -0.968160239507626 }, { 0,
// 0.968160239507626, -0.836031107326636 }, { 0, 0.968160239507626,
//-0.61337143270059 }, { 0, 0.968160239507626, -0.324253423403809 }, { 0,
// 0.968160239507626, 0 }, { 0, 0.968160239507626, 0.324253423403809 }, { 0,
// 0.968160239507626, 0.61337143270059 }, { 0, 0.968160239507626,
// 0.836031107326636 }, { 0, 0.968160239507626, 0.968160239507626 }, {
// 0.324253423403809, -0.968160239507626, -0.968160239507626 }, {
// 0.324253423403809, -0.968160239507626, -0.836031107326636 }, {
// 0.324253423403809, -0.968160239507626, -0.61337143270059 }, {
// 0.324253423403809, -0.968160239507626, -0.324253423403809 }, {
// 0.324253423403809, -0.968160239507626, 0 }, { 0.324253423403809,
//-0.968160239507626, 0.324253423403809 }, { 0.324253423403809,
//-0.968160239507626, 0.61337143270059 }, { 0.324253423403809,
//-0.968160239507626, 0.836031107326636 }, { 0.324253423403809,
//-0.968160239507626, 0.968160239507626 }, { 0.324253423403809,
//-0.836031107326636, -0.968160239507626 }, { 0.324253423403809,
//-0.836031107326636, -0.836031107326636 }, { 0.324253423403809,
//-0.836031107326636, -0.61337143270059 }, { 0.324253423403809,
//-0.836031107326636, -0.324253423403809 }, { 0.324253423403809,
//-0.836031107326636, 0 }, { 0.324253423403809, -0.836031107326636,
// 0.324253423403809 }, { 0.324253423403809, -0.836031107326636,
// 0.61337143270059
//}, { 0.324253423403809, -0.836031107326636, 0.836031107326636 }, {
// 0.324253423403809, -0.836031107326636, 0.968160239507626 }, {
// 0.324253423403809, -0.61337143270059, -0.968160239507626 }, {
// 0.324253423403809, -0.61337143270059, -0.836031107326636 }, {
// 0.324253423403809, -0.61337143270059, -0.61337143270059 }, {
// 0.324253423403809, -0.61337143270059, -0.324253423403809 }, {
// 0.324253423403809, -0.61337143270059, 0 }, { 0.324253423403809,
//-0.61337143270059, 0.324253423403809 }, { 0.324253423403809,
//-0.61337143270059, 0.61337143270059 }, { 0.324253423403809, -0.61337143270059,
// 0.836031107326636 }, { 0.324253423403809, -0.61337143270059,
// 0.968160239507626
//}, { 0.324253423403809, -0.324253423403809, -0.968160239507626 }, {
// 0.324253423403809, -0.324253423403809, -0.836031107326636 }, {
// 0.324253423403809, -0.324253423403809, -0.61337143270059 }, {
// 0.324253423403809, -0.324253423403809, -0.324253423403809 }, {
// 0.324253423403809, -0.324253423403809, 0 }, { 0.324253423403809,
//-0.324253423403809, 0.324253423403809 }, { 0.324253423403809,
//-0.324253423403809, 0.61337143270059 }, { 0.324253423403809,
//-0.324253423403809, 0.836031107326636 }, { 0.324253423403809,
//-0.324253423403809, 0.968160239507626 }, { 0.324253423403809, 0,
//-0.968160239507626 }, { 0.324253423403809, 0, -0.836031107326636 }, {
// 0.324253423403809, 0, -0.61337143270059 }, { 0.324253423403809, 0,
//-0.324253423403809 }, { 0.324253423403809, 0, 0 }, { 0.324253423403809, 0,
// 0.324253423403809 }, { 0.324253423403809, 0, 0.61337143270059 }, {
// 0.324253423403809, 0, 0.836031107326636 }, { 0.324253423403809, 0,
// 0.968160239507626 }, { 0.324253423403809, 0.324253423403809,
//-0.968160239507626 }, { 0.324253423403809, 0.324253423403809,
//-0.836031107326636 }, { 0.324253423403809, 0.324253423403809,
//-0.61337143270059 }, { 0.324253423403809, 0.324253423403809,
//-0.324253423403809 }, { 0.324253423403809, 0.324253423403809, 0 }, {
// 0.324253423403809, 0.324253423403809, 0.324253423403809 }, {
// 0.324253423403809, 0.324253423403809, 0.61337143270059 }, {
// 0.324253423403809, 0.324253423403809, 0.836031107326636 }, {
// 0.324253423403809, 0.324253423403809, 0.968160239507626 }, {
// 0.324253423403809, 0.61337143270059, -0.968160239507626 }, {
// 0.324253423403809, 0.61337143270059, -0.836031107326636 }, {
// 0.324253423403809, 0.61337143270059, -0.61337143270059
//}, { 0.324253423403809, 0.61337143270059, -0.324253423403809 }, {
// 0.324253423403809, 0.61337143270059, 0 }, { 0.324253423403809,
// 0.61337143270059, 0.324253423403809 }, { 0.324253423403809, 0.61337143270059,
// 0.61337143270059 }, { 0.324253423403809, 0.61337143270059, 0.836031107326636
//}, { 0.324253423403809, 0.61337143270059, 0.968160239507626 }, {
// 0.324253423403809, 0.836031107326636, -0.968160239507626 }, {
// 0.324253423403809, 0.836031107326636, -0.836031107326636 }, {
// 0.324253423403809, 0.836031107326636, -0.61337143270059 }, {
// 0.324253423403809, 0.836031107326636, -0.324253423403809 }, {
// 0.324253423403809, 0.836031107326636, 0 }, { 0.324253423403809,
// 0.836031107326636, 0.324253423403809 }, { 0.324253423403809,
// 0.836031107326636, 0.61337143270059 }, { 0.324253423403809,
// 0.836031107326636, 0.836031107326636 }, { 0.324253423403809,
// 0.836031107326636, 0.968160239507626
//}, { 0.324253423403809, 0.968160239507626, -0.968160239507626 }, {
// 0.324253423403809, 0.968160239507626, -0.836031107326636 }, {
// 0.324253423403809, 0.968160239507626, -0.61337143270059 }, {
// 0.324253423403809, 0.968160239507626, -0.324253423403809 }, {
// 0.324253423403809, 0.968160239507626, 0 }, { 0.324253423403809,
// 0.968160239507626, 0.324253423403809 }, { 0.324253423403809,
// 0.968160239507626, 0.61337143270059 }, { 0.324253423403809,
// 0.968160239507626, 0.836031107326636 }, { 0.324253423403809,
// 0.968160239507626, 0.968160239507626
//}, { 0.61337143270059, -0.968160239507626, -0.968160239507626 }, {
// 0.61337143270059, -0.968160239507626, -0.836031107326636 }, {
// 0.61337143270059, -0.968160239507626, -0.61337143270059 }, {
// 0.61337143270059, -0.968160239507626, -0.324253423403809 }, {
// 0.61337143270059, -0.968160239507626, 0 }, { 0.61337143270059,
//-0.968160239507626, 0.324253423403809 }, { 0.61337143270059,
// -0.968160239507626, 0.61337143270059
//}, { 0.61337143270059, -0.968160239507626, 0.836031107326636 }, {
// 0.61337143270059, -0.968160239507626, 0.968160239507626 }, {
// 0.61337143270059, -0.836031107326636, -0.968160239507626 }, {
// 0.61337143270059, -0.836031107326636, -0.836031107326636 }, {
// 0.61337143270059, -0.836031107326636, -0.61337143270059 }, {
// 0.61337143270059, -0.836031107326636, -0.324253423403809 }, {
// 0.61337143270059, -0.836031107326636, 0 }, { 0.61337143270059,
//-0.836031107326636,
// 0.324253423403809 }, { 0.61337143270059, -0.836031107326636, 0.61337143270059
//}, { 0.61337143270059, -0.836031107326636, 0.836031107326636 }, {
// 0.61337143270059, -0.836031107326636, 0.968160239507626 }, {
// 0.61337143270059, -0.61337143270059, -0.968160239507626 }, {
// 0.61337143270059, -0.61337143270059, -0.836031107326636 }, {
// 0.61337143270059, -0.61337143270059, -0.61337143270059 }, { 0.61337143270059,
//-0.61337143270059, -0.324253423403809 }, { 0.61337143270059,
//-0.61337143270059, 0 }, {
// 0.61337143270059, -0.61337143270059, 0.324253423403809 }, { 0.61337143270059,
//-0.61337143270059, 0.61337143270059 }, { 0.61337143270059, -0.61337143270059,
// 0.836031107326636 }, { 0.61337143270059, -0.61337143270059, 0.968160239507626
//}, { 0.61337143270059, -0.324253423403809, -0.968160239507626 }, {
// 0.61337143270059, -0.324253423403809, -0.836031107326636 }, {
// 0.61337143270059, -0.324253423403809, -0.61337143270059 }, {
// 0.61337143270059, -0.324253423403809, -0.324253423403809 }, {
// 0.61337143270059, -0.324253423403809, 0 }, { 0.61337143270059,
//-0.324253423403809, 0.324253423403809 }, { 0.61337143270059,
// -0.324253423403809, 0.61337143270059
//}, { 0.61337143270059, -0.324253423403809, 0.836031107326636 }, {
// 0.61337143270059, -0.324253423403809, 0.968160239507626 }, {
// 0.61337143270059, 0, -0.968160239507626 }, { 0.61337143270059, 0,
// -0.836031107326636 }, { 0.61337143270059, 0, -0.61337143270059 }, {
// 0.61337143270059, 0, -0.324253423403809 }, { 0.61337143270059, 0, 0 }, {
// 0.61337143270059, 0, 0.324253423403809 }, { 0.61337143270059, 0,
// 0.61337143270059 }, { 0.61337143270059, 0, 0.836031107326636 }, {
// 0.61337143270059, 0, 0.968160239507626 }, { 0.61337143270059,
// 0.324253423403809, -0.968160239507626
//}, { 0.61337143270059, 0.324253423403809, -0.836031107326636 }, {
// 0.61337143270059, 0.324253423403809, -0.61337143270059 }, { 0.61337143270059,
// 0.324253423403809, -0.324253423403809 }, { 0.61337143270059,
// 0.324253423403809, 0 }, { 0.61337143270059, 0.324253423403809,
// 0.324253423403809 }, { 0.61337143270059, 0.324253423403809, 0.61337143270059
//}, { 0.61337143270059, 0.324253423403809, 0.836031107326636 }, {
// 0.61337143270059, 0.324253423403809, 0.968160239507626 }, { 0.61337143270059,
// 0.61337143270059, -0.968160239507626 }, { 0.61337143270059, 0.61337143270059,
//-0.836031107326636 }, { 0.61337143270059, 0.61337143270059, -0.61337143270059
//}, { 0.61337143270059, 0.61337143270059, -0.324253423403809 }, {
// 0.61337143270059, 0.61337143270059, 0 }, { 0.61337143270059,
// 0.61337143270059, 0.324253423403809 }, { 0.61337143270059, 0.61337143270059,
// 0.61337143270059 }, { 0.61337143270059, 0.61337143270059, 0.836031107326636
//}, { 0.61337143270059, 0.61337143270059, 0.968160239507626 }, {
// 0.61337143270059, 0.836031107326636, -0.968160239507626 }, {
// 0.61337143270059, 0.836031107326636, -0.836031107326636 }, {
// 0.61337143270059, 0.836031107326636, -0.61337143270059
//}, { 0.61337143270059, 0.836031107326636, -0.324253423403809 }, {
// 0.61337143270059, 0.836031107326636, 0 }, { 0.61337143270059,
// 0.836031107326636, 0.324253423403809 }, { 0.61337143270059,
// 0.836031107326636, 0.61337143270059 }, { 0.61337143270059, 0.836031107326636,
// 0.836031107326636
//}, { 0.61337143270059, 0.836031107326636, 0.968160239507626 }, {
// 0.61337143270059, 0.968160239507626, -0.968160239507626 }, {
// 0.61337143270059, 0.968160239507626, -0.836031107326636 }, {
// 0.61337143270059, 0.968160239507626, -0.61337143270059 }, { 0.61337143270059,
// 0.968160239507626, -0.324253423403809 }, { 0.61337143270059,
// 0.968160239507626, 0 }, { 0.61337143270059, 0.968160239507626,
// 0.324253423403809 }, { 0.61337143270059, 0.968160239507626, 0.61337143270059
// }, { 0.61337143270059, 0.968160239507626, 0.836031107326636 }, {
// 0.61337143270059, 0.968160239507626, 0.968160239507626
//}, { 0.836031107326636, -0.968160239507626, -0.968160239507626 }, {
// 0.836031107326636, -0.968160239507626, -0.836031107326636 }, {
// 0.836031107326636, -0.968160239507626, -0.61337143270059 }, {
// 0.836031107326636, -0.968160239507626, -0.324253423403809 }, {
// 0.836031107326636, -0.968160239507626, 0 }, { 0.836031107326636,
//-0.968160239507626, 0.324253423403809 }, { 0.836031107326636,
//-0.968160239507626, 0.61337143270059 }, { 0.836031107326636,
//-0.968160239507626, 0.836031107326636 }, { 0.836031107326636,
//-0.968160239507626, 0.968160239507626 }, { 0.836031107326636,
//-0.836031107326636, -0.968160239507626 }, { 0.836031107326636,
//-0.836031107326636, -0.836031107326636 }, { 0.836031107326636,
//-0.836031107326636, -0.61337143270059 }, { 0.836031107326636,
//-0.836031107326636, -0.324253423403809 }, { 0.836031107326636,
//-0.836031107326636, 0 }, { 0.836031107326636, -0.836031107326636,
// 0.324253423403809 }, { 0.836031107326636, -0.836031107326636,
// 0.61337143270059
//}, { 0.836031107326636, -0.836031107326636, 0.836031107326636 }, {
// 0.836031107326636, -0.836031107326636, 0.968160239507626 }, {
// 0.836031107326636, -0.61337143270059, -0.968160239507626 }, {
// 0.836031107326636, -0.61337143270059, -0.836031107326636 }, {
// 0.836031107326636, -0.61337143270059, -0.61337143270059 }, {
// 0.836031107326636, -0.61337143270059, -0.324253423403809 }, {
// 0.836031107326636, -0.61337143270059, 0 }, { 0.836031107326636,
//-0.61337143270059, 0.324253423403809 }, { 0.836031107326636,
//-0.61337143270059, 0.61337143270059 }, { 0.836031107326636, -0.61337143270059,
// 0.836031107326636 }, { 0.836031107326636, -0.61337143270059,
// 0.968160239507626
//}, { 0.836031107326636, -0.324253423403809, -0.968160239507626 }, {
// 0.836031107326636, -0.324253423403809, -0.836031107326636 }, {
// 0.836031107326636, -0.324253423403809, -0.61337143270059 }, {
// 0.836031107326636, -0.324253423403809, -0.324253423403809 }, {
// 0.836031107326636, -0.324253423403809, 0 }, { 0.836031107326636,
//-0.324253423403809, 0.324253423403809 }, { 0.836031107326636,
//-0.324253423403809, 0.61337143270059 }, { 0.836031107326636,
//-0.324253423403809, 0.836031107326636 }, { 0.836031107326636,
//-0.324253423403809, 0.968160239507626 }, { 0.836031107326636, 0,
//-0.968160239507626 }, { 0.836031107326636, 0, -0.836031107326636 }, {
// 0.836031107326636, 0, -0.61337143270059 }, { 0.836031107326636, 0,
//-0.324253423403809 }, { 0.836031107326636, 0, 0 }, { 0.836031107326636, 0,
// 0.324253423403809 }, { 0.836031107326636, 0, 0.61337143270059 }, {
// 0.836031107326636, 0, 0.836031107326636 }, { 0.836031107326636, 0,
// 0.968160239507626 }, { 0.836031107326636, 0.324253423403809,
//-0.968160239507626 }, { 0.836031107326636, 0.324253423403809,
//-0.836031107326636 }, { 0.836031107326636, 0.324253423403809,
//-0.61337143270059 }, { 0.836031107326636, 0.324253423403809,
//-0.324253423403809 }, { 0.836031107326636, 0.324253423403809, 0 }, {
// 0.836031107326636, 0.324253423403809, 0.324253423403809 }, {
// 0.836031107326636, 0.324253423403809, 0.61337143270059 }, {
// 0.836031107326636, 0.324253423403809, 0.836031107326636 }, {
// 0.836031107326636, 0.324253423403809, 0.968160239507626 }, {
// 0.836031107326636, 0.61337143270059, -0.968160239507626 }, {
// 0.836031107326636, 0.61337143270059, -0.836031107326636 }, {
// 0.836031107326636, 0.61337143270059, -0.61337143270059
//}, { 0.836031107326636, 0.61337143270059, -0.324253423403809 }, {
// 0.836031107326636, 0.61337143270059, 0 }, { 0.836031107326636,
// 0.61337143270059, 0.324253423403809 }, { 0.836031107326636, 0.61337143270059,
// 0.61337143270059 }, { 0.836031107326636, 0.61337143270059, 0.836031107326636
//}, { 0.836031107326636, 0.61337143270059, 0.968160239507626 }, {
// 0.836031107326636, 0.836031107326636, -0.968160239507626 }, {
// 0.836031107326636, 0.836031107326636, -0.836031107326636 }, {
// 0.836031107326636, 0.836031107326636, -0.61337143270059 }, {
// 0.836031107326636, 0.836031107326636, -0.324253423403809 }, {
// 0.836031107326636, 0.836031107326636, 0 }, { 0.836031107326636,
// 0.836031107326636, 0.324253423403809 }, { 0.836031107326636,
// 0.836031107326636, 0.61337143270059 }, { 0.836031107326636,
// 0.836031107326636, 0.836031107326636 }, { 0.836031107326636,
// 0.836031107326636, 0.968160239507626
//}, { 0.836031107326636, 0.968160239507626, -0.968160239507626 }, {
// 0.836031107326636, 0.968160239507626, -0.836031107326636 }, {
// 0.836031107326636, 0.968160239507626, -0.61337143270059 }, {
// 0.836031107326636, 0.968160239507626, -0.324253423403809 }, {
// 0.836031107326636, 0.968160239507626, 0 }, { 0.836031107326636,
// 0.968160239507626, 0.324253423403809 }, { 0.836031107326636,
// 0.968160239507626, 0.61337143270059 }, { 0.836031107326636,
// 0.968160239507626, 0.836031107326636 }, { 0.836031107326636,
// 0.968160239507626, 0.968160239507626
//}, { 0.968160239507626, -0.968160239507626, -0.968160239507626 }, {
// 0.968160239507626, -0.968160239507626, -0.836031107326636 }, {
// 0.968160239507626, -0.968160239507626, -0.61337143270059 }, {
// 0.968160239507626, -0.968160239507626, -0.324253423403809 }, {
// 0.968160239507626, -0.968160239507626, 0 }, { 0.968160239507626,
//-0.968160239507626, 0.324253423403809 }, { 0.968160239507626,
//-0.968160239507626, 0.61337143270059 }, { 0.968160239507626,
//-0.968160239507626, 0.836031107326636 }, { 0.968160239507626,
//-0.968160239507626, 0.968160239507626 }, { 0.968160239507626,
//-0.836031107326636, -0.968160239507626 }, { 0.968160239507626,
//-0.836031107326636, -0.836031107326636 }, { 0.968160239507626,
//-0.836031107326636, -0.61337143270059 }, { 0.968160239507626,
//-0.836031107326636, -0.324253423403809 }, { 0.968160239507626,
//-0.836031107326636, 0 }, { 0.968160239507626, -0.836031107326636,
// 0.324253423403809 }, { 0.968160239507626, -0.836031107326636,
// 0.61337143270059
//}, { 0.968160239507626, -0.836031107326636, 0.836031107326636 }, {
// 0.968160239507626, -0.836031107326636, 0.968160239507626 }, {
// 0.968160239507626, -0.61337143270059, -0.968160239507626 }, {
// 0.968160239507626, -0.61337143270059, -0.836031107326636 }, {
// 0.968160239507626, -0.61337143270059, -0.61337143270059 }, {
// 0.968160239507626, -0.61337143270059, -0.324253423403809 }, {
// 0.968160239507626, -0.61337143270059, 0 }, { 0.968160239507626,
//-0.61337143270059, 0.324253423403809 }, { 0.968160239507626,
//-0.61337143270059, 0.61337143270059 }, { 0.968160239507626, -0.61337143270059,
// 0.836031107326636 }, { 0.968160239507626, -0.61337143270059,
// 0.968160239507626
//}, { 0.968160239507626, -0.324253423403809, -0.968160239507626 }, {
// 0.968160239507626, -0.324253423403809, -0.836031107326636 }, {
// 0.968160239507626, -0.324253423403809, -0.61337143270059 }, {
// 0.968160239507626, -0.324253423403809, -0.324253423403809 }, {
// 0.968160239507626, -0.324253423403809, 0 }, { 0.968160239507626,
//-0.324253423403809, 0.324253423403809 }, { 0.968160239507626,
//-0.324253423403809, 0.61337143270059 }, { 0.968160239507626,
//-0.324253423403809, 0.836031107326636 }, { 0.968160239507626,
//-0.324253423403809, 0.968160239507626 }, { 0.968160239507626, 0,
//-0.968160239507626 }, { 0.968160239507626, 0, -0.836031107326636 }, {
// 0.968160239507626, 0, -0.61337143270059 }, { 0.968160239507626, 0,
//-0.324253423403809 }, { 0.968160239507626, 0, 0 }, { 0.968160239507626, 0,
// 0.324253423403809 }, { 0.968160239507626, 0, 0.61337143270059 }, {
// 0.968160239507626, 0, 0.836031107326636 }, { 0.968160239507626, 0,
// 0.968160239507626 }, { 0.968160239507626, 0.324253423403809,
//-0.968160239507626 }, { 0.968160239507626, 0.324253423403809,
//-0.836031107326636 }, { 0.968160239507626, 0.324253423403809,
//-0.61337143270059 }, { 0.968160239507626, 0.324253423403809,
//-0.324253423403809 }, { 0.968160239507626, 0.324253423403809, 0 }, {
// 0.968160239507626, 0.324253423403809, 0.324253423403809 }, {
// 0.968160239507626, 0.324253423403809, 0.61337143270059 }, {
// 0.968160239507626, 0.324253423403809, 0.836031107326636 }, {
// 0.968160239507626, 0.324253423403809, 0.968160239507626 }, {
// 0.968160239507626, 0.61337143270059, -0.968160239507626 }, {
// 0.968160239507626, 0.61337143270059, -0.836031107326636 }, {
// 0.968160239507626, 0.61337143270059, -0.61337143270059
//}, { 0.968160239507626, 0.61337143270059, -0.324253423403809 }, {
// 0.968160239507626, 0.61337143270059, 0 }, { 0.968160239507626,
// 0.61337143270059, 0.324253423403809 }, { 0.968160239507626, 0.61337143270059,
// 0.61337143270059 }, { 0.968160239507626, 0.61337143270059, 0.836031107326636
//}, { 0.968160239507626, 0.61337143270059, 0.968160239507626 }, {
// 0.968160239507626, 0.836031107326636, -0.968160239507626 }, {
// 0.968160239507626, 0.836031107326636, -0.836031107326636 }, {
// 0.968160239507626, 0.836031107326636, -0.61337143270059 }, {
// 0.968160239507626, 0.836031107326636, -0.324253423403809 }, {
// 0.968160239507626, 0.836031107326636, 0 }, { 0.968160239507626,
// 0.836031107326636, 0.324253423403809 }, { 0.968160239507626,
// 0.836031107326636, 0.61337143270059 }, { 0.968160239507626,
// 0.836031107326636, 0.836031107326636 }, { 0.968160239507626,
// 0.836031107326636, 0.968160239507626
//}, { 0.968160239507626, 0.968160239507626, -0.968160239507626 }, {
// 0.968160239507626, 0.968160239507626, -0.836031107326636 }, {
// 0.968160239507626, 0.968160239507626, -0.61337143270059 }, {
// 0.968160239507626, 0.968160239507626, -0.324253423403809 }, {
// 0.968160239507626, 0.968160239507626, 0 }, { 0.968160239507626,
// 0.968160239507626, 0.324253423403809 }, { 0.968160239507626,
// 0.968160239507626, 0.61337143270059 }, { 0.968160239507626,
// 0.968160239507626, 0.836031107326636 }, { 0.968160239507626,
// 0.968160239507626, 0.968160239507626 } } , { 0.0005368601019997297,
// 0.001193276159092649, 0.001721470784014517, 0.002063216801989402,
// 0.002181404712903652, 0.002063216801989402, 0.001721470784014517,
// 0.001193276159092649, 0.0005368601019997297, 0.001193276159092649,
// 0.002652288718336573, 0.003826304166555643, 0.004585901264934326,
// 0.004848596920397677, 0.004585901264934326, 0.003826304166555643,
// 0.002652288718336573, 0.001193276159092649, 0.001721470784014517,
// 0.003826304166555643, 0.005519988632377539, 0.006615815614687656,
// 0.006994791505994799, 0.006615815614687656, 0.005519988632377539,
// 0.003826304166555643, 0.001721470784014517, 0.002063216801989402,
// 0.004585901264934326, 0.006615815614687656, 0.007929185939046591,
// 0.008383396044588071, 0.007929185939046591, 0.006615815614687656,
// 0.004585901264934326, 0.002063216801989402, 0.002181404712903652,
// 0.004848596920397677, 0.006994791505994799, 0.008383396044588071,
// 0.008863624813528533, 0.008383396044588071, 0.006994791505994799,
// 0.004848596920397677, 0.002181404712903652, 0.002063216801989402,
// 0.004585901264934326, 0.006615815614687656, 0.007929185939046591,
// 0.008383396044588071, 0.007929185939046591, 0.006615815614687656,
// 0.004585901264934326, 0.002063216801989402, 0.001721470784014517,
// 0.003826304166555643, 0.005519988632377539, 0.006615815614687656,
// 0.006994791505994799, 0.006615815614687656, 0.005519988632377539,
// 0.003826304166555643, 0.001721470784014517, 0.001193276159092649,
// 0.002652288718336573, 0.003826304166555643, 0.004585901264934326,
// 0.004848596920397677, 0.004585901264934326, 0.003826304166555643,
// 0.002652288718336573, 0.001193276159092649, 0.0005368601019997297,
// 0.001193276159092649, 0.001721470784014517, 0.002063216801989402,
// 0.002181404712903652, 0.002063216801989402, 0.001721470784014517,
// 0.001193276159092649, 0.0005368601019997297, 0.001193276159092649,
// 0.002652288718336573, 0.003826304166555643, 0.004585901264934326,
// 0.004848596920397677, 0.004585901264934326, 0.003826304166555643,
// 0.002652288718336573, 0.001193276159092649, 0.002652288718336573,
// 0.005895228352475006, 0.008504706388835012, 0.01019305891239746,
// 0.01077695118078173, 0.01019305891239746, 0.008504706388835012,
// 0.005895228352475006, 0.002652288718336573, 0.003826304166555643,
// 0.008504706388835012, 0.01226925005032665, 0.01470493898979223,
// 0.01554728673417337, 0.01470493898979223, 0.01226925005032665,
// 0.008504706388835012, 0.003826304166555643, 0.004585901264934327,
// 0.01019305891239746, 0.01470493898979223, 0.01762416038523518,
// 0.0186337308266644, 0.01762416038523518, 0.01470493898979223,
// 0.01019305891239746, 0.004585901264934327, 0.004848596920397677,
// 0.01077695118078173, 0.01554728673417337, 0.0186337308266644,
// 0.01970113281603286, 0.0186337308266644, 0.01554728673417337,
// 0.01077695118078173, 0.004848596920397677, 0.004585901264934327,
// 0.01019305891239746, 0.01470493898979223, 0.01762416038523518,
// 0.0186337308266644, 0.01762416038523518, 0.01470493898979223,
// 0.01019305891239746, 0.004585901264934327, 0.003826304166555643,
// 0.008504706388835012, 0.01226925005032665, 0.01470493898979223,
// 0.01554728673417337, 0.01470493898979223, 0.01226925005032665,
// 0.008504706388835012, 0.003826304166555643, 0.002652288718336573,
// 0.005895228352475006, 0.008504706388835012, 0.01019305891239746,
// 0.01077695118078173, 0.01019305891239746, 0.008504706388835012,
// 0.005895228352475006, 0.002652288718336573, 0.001193276159092649,
// 0.002652288718336573, 0.003826304166555643, 0.004585901264934326,
// 0.004848596920397677, 0.004585901264934326, 0.003826304166555643,
// 0.002652288718336573, 0.001193276159092649, 0.001721470784014517,
// 0.003826304166555643, 0.005519988632377539, 0.006615815614687656,
// 0.006994791505994799, 0.006615815614687656, 0.005519988632377539,
// 0.003826304166555643, 0.001721470784014517, 0.003826304166555643,
// 0.008504706388835012, 0.01226925005032665, 0.01470493898979223,
// 0.01554728673417337, 0.01470493898979223, 0.01226925005032665,
// 0.008504706388835012, 0.003826304166555643, 0.00551998863237754,
// 0.01226925005032665, 0.01770014035935002, 0.02121396850071301,
// 0.02242917507371224, 0.02121396850071301, 0.01770014035935002,
// 0.01226925005032665, 0.00551998863237754, 0.006615815614687656,
// 0.01470493898979223, 0.02121396850071301, 0.02542536106565484,
// 0.02688181019193803, 0.02542536106565484, 0.02121396850071301,
// 0.01470493898979223, 0.006615815614687656, 0.0069947915059948,
// 0.01554728673417337, 0.02242917507371224, 0.02688181019193803,
// 0.02842168956143284, 0.02688181019193803, 0.02242917507371224,
// 0.01554728673417337, 0.0069947915059948, 0.006615815614687656,
// 0.01470493898979223, 0.02121396850071301, 0.02542536106565484,
// 0.02688181019193803, 0.02542536106565484, 0.02121396850071301,
// 0.01470493898979223, 0.006615815614687656, 0.00551998863237754,
// 0.01226925005032665, 0.01770014035935002, 0.02121396850071301,
// 0.02242917507371224, 0.02121396850071301, 0.01770014035935002,
// 0.01226925005032665, 0.00551998863237754, 0.003826304166555643,
// 0.008504706388835012, 0.01226925005032665, 0.01470493898979223,
// 0.01554728673417337, 0.01470493898979223, 0.01226925005032665,
// 0.008504706388835012, 0.003826304166555643, 0.001721470784014517,
// 0.003826304166555643, 0.005519988632377539, 0.006615815614687656,
// 0.006994791505994799, 0.006615815614687656, 0.005519988632377539,
// 0.003826304166555643, 0.001721470784014517, 0.002063216801989402,
// 0.004585901264934326, 0.006615815614687656, 0.007929185939046591,
// 0.008383396044588071, 0.007929185939046591, 0.006615815614687656,
// 0.004585901264934326, 0.002063216801989402, 0.004585901264934327,
// 0.01019305891239746, 0.01470493898979223, 0.01762416038523518,
// 0.0186337308266644, 0.01762416038523518, 0.01470493898979223,
// 0.01019305891239746, 0.004585901264934327, 0.006615815614687656,
// 0.01470493898979223, 0.02121396850071301, 0.02542536106565484,
// 0.02688181019193803, 0.02542536106565484, 0.02121396850071301,
// 0.01470493898979223, 0.006615815614687656, 0.007929185939046591,
// 0.01762416038523518, 0.02542536106565484, 0.03047279839682941,
// 0.03221838149733523, 0.03047279839682941, 0.02542536106565484,
// 0.01762416038523518, 0.007929185939046591, 0.008383396044588071,
// 0.0186337308266644, 0.02688181019193803, 0.03221838149733524,
// 0.03406395739538762, 0.03221838149733524, 0.02688181019193803,
// 0.0186337308266644, 0.008383396044588071, 0.007929185939046591,
// 0.01762416038523518, 0.02542536106565484, 0.03047279839682941,
// 0.03221838149733523, 0.03047279839682941, 0.02542536106565484,
// 0.01762416038523518, 0.007929185939046591, 0.006615815614687656,
// 0.01470493898979223, 0.02121396850071301, 0.02542536106565484,
// 0.02688181019193803, 0.02542536106565484, 0.02121396850071301,
// 0.01470493898979223, 0.006615815614687656, 0.004585901264934327,
// 0.01019305891239746, 0.01470493898979223, 0.01762416038523518,
// 0.0186337308266644, 0.01762416038523518, 0.01470493898979223,
// 0.01019305891239746, 0.004585901264934327, 0.002063216801989402,
// 0.004585901264934326, 0.006615815614687656, 0.007929185939046591,
// 0.008383396044588071, 0.007929185939046591, 0.006615815614687656,
// 0.004585901264934326, 0.002063216801989402, 0.002181404712903652,
// 0.004848596920397677, 0.006994791505994799, 0.008383396044588071,
// 0.008863624813528533, 0.008383396044588071, 0.006994791505994799,
// 0.004848596920397677, 0.002181404712903652, 0.004848596920397677,
// 0.01077695118078173, 0.01554728673417337, 0.0186337308266644,
// 0.01970113281603286, 0.0186337308266644, 0.01554728673417337,
// 0.01077695118078173, 0.004848596920397677, 0.0069947915059948,
// 0.01554728673417337, 0.02242917507371224, 0.02688181019193803,
// 0.02842168956143284, 0.02688181019193803, 0.02242917507371224,
// 0.01554728673417337, 0.0069947915059948, 0.008383396044588071,
// 0.0186337308266644, 0.02688181019193803, 0.03221838149733524,
// 0.03406395739538762, 0.03221838149733524, 0.02688181019193803,
// 0.0186337308266644, 0.008383396044588071, 0.008863624813528533,
// 0.01970113281603286, 0.02842168956143284, 0.03406395739538763,
// 0.03601525401053295, 0.03406395739538763, 0.02842168956143284,
// 0.01970113281603286, 0.008863624813528533, 0.008383396044588071,
// 0.0186337308266644, 0.02688181019193803, 0.03221838149733524,
// 0.03406395739538762, 0.03221838149733524, 0.02688181019193803,
// 0.0186337308266644, 0.008383396044588071, 0.0069947915059948,
// 0.01554728673417337, 0.02242917507371224, 0.02688181019193803,
// 0.02842168956143284, 0.02688181019193803, 0.02242917507371224,
// 0.01554728673417337, 0.0069947915059948, 0.004848596920397677,
// 0.01077695118078173, 0.01554728673417337, 0.0186337308266644,
// 0.01970113281603286, 0.0186337308266644, 0.01554728673417337,
// 0.01077695118078173, 0.004848596920397677, 0.002181404712903652,
// 0.004848596920397677, 0.006994791505994799, 0.008383396044588071,
// 0.008863624813528533, 0.008383396044588071, 0.006994791505994799,
// 0.004848596920397677, 0.002181404712903652, 0.002063216801989402,
// 0.004585901264934326, 0.006615815614687656, 0.007929185939046591,
// 0.008383396044588071, 0.007929185939046591, 0.006615815614687656,
// 0.004585901264934326, 0.002063216801989402, 0.004585901264934327,
// 0.01019305891239746, 0.01470493898979223, 0.01762416038523518,
// 0.0186337308266644, 0.01762416038523518, 0.01470493898979223,
// 0.01019305891239746, 0.004585901264934327, 0.006615815614687656,
// 0.01470493898979223, 0.02121396850071301, 0.02542536106565484,
// 0.02688181019193803, 0.02542536106565484, 0.02121396850071301,
// 0.01470493898979223, 0.006615815614687656, 0.007929185939046591,
// 0.01762416038523518, 0.02542536106565484, 0.03047279839682941,
// 0.03221838149733523, 0.03047279839682941, 0.02542536106565484,
// 0.01762416038523518, 0.007929185939046591, 0.008383396044588071,
// 0.0186337308266644, 0.02688181019193803, 0.03221838149733524,
// 0.03406395739538762, 0.03221838149733524, 0.02688181019193803,
// 0.0186337308266644, 0.008383396044588071, 0.007929185939046591,
// 0.01762416038523518, 0.02542536106565484, 0.03047279839682941,
// 0.03221838149733523, 0.03047279839682941, 0.02542536106565484,
// 0.01762416038523518, 0.007929185939046591, 0.006615815614687656,
// 0.01470493898979223, 0.02121396850071301, 0.02542536106565484,
// 0.02688181019193803, 0.02542536106565484, 0.02121396850071301,
// 0.01470493898979223, 0.006615815614687656, 0.004585901264934327,
// 0.01019305891239746, 0.01470493898979223, 0.01762416038523518,
// 0.0186337308266644, 0.01762416038523518, 0.01470493898979223,
// 0.01019305891239746, 0.004585901264934327, 0.002063216801989402,
// 0.004585901264934326, 0.006615815614687656, 0.007929185939046591,
// 0.008383396044588071, 0.007929185939046591, 0.006615815614687656,
// 0.004585901264934326, 0.002063216801989402, 0.001721470784014517,
// 0.003826304166555643, 0.005519988632377539, 0.006615815614687656,
// 0.006994791505994799, 0.006615815614687656, 0.005519988632377539,
// 0.003826304166555643, 0.001721470784014517, 0.003826304166555643,
// 0.008504706388835012, 0.01226925005032665, 0.01470493898979223,
// 0.01554728673417337, 0.01470493898979223, 0.01226925005032665,
// 0.008504706388835012, 0.003826304166555643, 0.00551998863237754,
// 0.01226925005032665, 0.01770014035935002, 0.02121396850071301,
// 0.02242917507371224, 0.02121396850071301, 0.01770014035935002,
// 0.01226925005032665, 0.00551998863237754, 0.006615815614687656,
// 0.01470493898979223, 0.02121396850071301, 0.02542536106565484,
// 0.02688181019193803, 0.02542536106565484, 0.02121396850071301,
// 0.01470493898979223, 0.006615815614687656, 0.0069947915059948,
// 0.01554728673417337, 0.02242917507371224, 0.02688181019193803,
// 0.02842168956143284, 0.02688181019193803, 0.02242917507371224,
// 0.01554728673417337, 0.0069947915059948, 0.006615815614687656,
// 0.01470493898979223, 0.02121396850071301, 0.02542536106565484,
// 0.02688181019193803, 0.02542536106565484, 0.02121396850071301,
// 0.01470493898979223, 0.006615815614687656, 0.00551998863237754,
// 0.01226925005032665, 0.01770014035935002, 0.02121396850071301,
// 0.02242917507371224, 0.02121396850071301, 0.01770014035935002,
// 0.01226925005032665, 0.00551998863237754, 0.003826304166555643,
// 0.008504706388835012, 0.01226925005032665, 0.01470493898979223,
// 0.01554728673417337, 0.01470493898979223, 0.01226925005032665,
// 0.008504706388835012, 0.003826304166555643, 0.001721470784014517,
// 0.003826304166555643, 0.005519988632377539, 0.006615815614687656,
// 0.006994791505994799, 0.006615815614687656, 0.005519988632377539,
// 0.003826304166555643, 0.001721470784014517, 0.001193276159092649,
// 0.002652288718336573, 0.003826304166555643, 0.004585901264934326,
// 0.004848596920397677, 0.004585901264934326, 0.003826304166555643,
// 0.002652288718336573, 0.001193276159092649, 0.002652288718336573,
// 0.005895228352475006, 0.008504706388835012, 0.01019305891239746,
// 0.01077695118078173, 0.01019305891239746, 0.008504706388835012,
// 0.005895228352475006, 0.002652288718336573, 0.003826304166555643,
// 0.008504706388835012, 0.01226925005032665, 0.01470493898979223,
// 0.01554728673417337, 0.01470493898979223, 0.01226925005032665,
// 0.008504706388835012, 0.003826304166555643, 0.004585901264934327,
// 0.01019305891239746, 0.01470493898979223, 0.01762416038523518,
// 0.0186337308266644, 0.01762416038523518, 0.01470493898979223,
// 0.01019305891239746, 0.004585901264934327, 0.004848596920397677,
// 0.01077695118078173, 0.01554728673417337, 0.0186337308266644,
// 0.01970113281603286, 0.0186337308266644, 0.01554728673417337,
// 0.01077695118078173, 0.004848596920397677, 0.004585901264934327,
// 0.01019305891239746, 0.01470493898979223, 0.01762416038523518,
// 0.0186337308266644, 0.01762416038523518, 0.01470493898979223,
// 0.01019305891239746, 0.004585901264934327, 0.003826304166555643,
// 0.008504706388835012, 0.01226925005032665, 0.01470493898979223,
// 0.01554728673417337, 0.01470493898979223, 0.01226925005032665,
// 0.008504706388835012, 0.003826304166555643, 0.002652288718336573,
// 0.005895228352475006, 0.008504706388835012, 0.01019305891239746,
// 0.01077695118078173, 0.01019305891239746, 0.008504706388835012,
// 0.005895228352475006, 0.002652288718336573, 0.001193276159092649,
// 0.002652288718336573, 0.003826304166555643, 0.004585901264934326,
// 0.004848596920397677, 0.004585901264934326, 0.003826304166555643,
// 0.002652288718336573, 0.001193276159092649, 0.0005368601019997297,
// 0.001193276159092649, 0.001721470784014517, 0.002063216801989402,
// 0.002181404712903652, 0.002063216801989402, 0.001721470784014517,
// 0.001193276159092649, 0.0005368601019997297, 0.001193276159092649,
// 0.002652288718336573, 0.003826304166555643, 0.004585901264934326,
// 0.004848596920397677, 0.004585901264934326, 0.003826304166555643,
// 0.002652288718336573, 0.001193276159092649, 0.001721470784014517,
// 0.003826304166555643, 0.005519988632377539, 0.006615815614687656,
// 0.006994791505994799, 0.006615815614687656, 0.005519988632377539,
// 0.003826304166555643, 0.001721470784014517, 0.002063216801989402,
// 0.004585901264934326, 0.006615815614687656, 0.007929185939046591,
// 0.008383396044588071, 0.007929185939046591, 0.006615815614687656,
// 0.004585901264934326, 0.002063216801989402, 0.002181404712903652,
// 0.004848596920397677, 0.006994791505994799, 0.008383396044588071,
// 0.008863624813528533, 0.008383396044588071, 0.006994791505994799,
// 0.004848596920397677, 0.002181404712903652, 0.002063216801989402,
// 0.004585901264934326, 0.006615815614687656, 0.007929185939046591,
// 0.008383396044588071, 0.007929185939046591, 0.006615815614687656,
// 0.004585901264934326, 0.002063216801989402, 0.001721470784014517,
// 0.003826304166555643, 0.005519988632377539, 0.006615815614687656,
// 0.006994791505994799, 0.006615815614687656, 0.005519988632377539,
// 0.003826304166555643, 0.001721470784014517, 0.001193276159092649,
// 0.002652288718336573, 0.003826304166555643, 0.004585901264934326,
// 0.004848596920397677, 0.004585901264934326, 0.003826304166555643,
// 0.002652288718336573, 0.001193276159092649, 0.0005368601019997297,
// 0.001193276159092649, 0.001721470784014517, 0.002063216801989402,
// 0.002181404712903652, 0.002063216801989402, 0.001721470784014517,
// 0.001193276159092649, 0.0005368601019997297 } };
// case 18: 					case 19:	return { { {
// -0.973906528517172, -0.973906528517172, -0.973906528517172 }, {
//-0.973906528517172, -0.973906528517172, -0.865063366688985 }, {
//-0.973906528517172, -0.973906528517172, -0.679409568299024 }, {
//-0.973906528517172, -0.973906528517172, -0.433395394129247 }, {
//-0.973906528517172, -0.973906528517172, -0.148874338981631 }, {
//-0.973906528517172, -0.973906528517172, 0.148874338981631 }, {
// -0.973906528517172, -0.973906528517172, 0.433395394129247 }, {
// -0.973906528517172, -0.973906528517172, 0.679409568299024 }, {
// -0.973906528517172, -0.973906528517172, 0.865063366688985 }, {
// -0.973906528517172, -0.973906528517172, 0.973906528517172 }, {
// -0.973906528517172, -0.865063366688985, -0.973906528517172 }, {
//-0.973906528517172, -0.865063366688985, -0.865063366688985 }, {
//-0.973906528517172, -0.865063366688985, -0.679409568299024 }, {
//-0.973906528517172, -0.865063366688985, -0.433395394129247 }, {
//-0.973906528517172, -0.865063366688985, -0.148874338981631 }, {
//-0.973906528517172, -0.865063366688985, 0.148874338981631 }, {
// -0.973906528517172, -0.865063366688985, 0.433395394129247 }, {
// -0.973906528517172, -0.865063366688985, 0.679409568299024 }, {
// -0.973906528517172, -0.865063366688985, 0.865063366688985 }, {
// -0.973906528517172, -0.865063366688985, 0.973906528517172 }, {
// -0.973906528517172, -0.679409568299024, -0.973906528517172 }, {
//-0.973906528517172, -0.679409568299024, -0.865063366688985 }, {
//-0.973906528517172, -0.679409568299024, -0.679409568299024 }, {
//-0.973906528517172, -0.679409568299024, -0.433395394129247 }, {
//-0.973906528517172, -0.679409568299024, -0.148874338981631 }, {
//-0.973906528517172, -0.679409568299024, 0.148874338981631 }, {
// -0.973906528517172, -0.679409568299024, 0.433395394129247 }, {
// -0.973906528517172, -0.679409568299024, 0.679409568299024 }, {
// -0.973906528517172, -0.679409568299024, 0.865063366688985 }, {
// -0.973906528517172, -0.679409568299024, 0.973906528517172 }, {
// -0.973906528517172, -0.433395394129247, -0.973906528517172 }, {
//-0.973906528517172, -0.433395394129247, -0.865063366688985 }, {
//-0.973906528517172, -0.433395394129247, -0.679409568299024 }, {
//-0.973906528517172, -0.433395394129247, -0.433395394129247 }, {
//-0.973906528517172, -0.433395394129247, -0.148874338981631 }, {
//-0.973906528517172, -0.433395394129247, 0.148874338981631 }, {
// -0.973906528517172, -0.433395394129247, 0.433395394129247 }, {
// -0.973906528517172, -0.433395394129247, 0.679409568299024 }, {
// -0.973906528517172, -0.433395394129247, 0.865063366688985 }, {
// -0.973906528517172, -0.433395394129247, 0.973906528517172 }, {
// -0.973906528517172, -0.148874338981631, -0.973906528517172 }, {
//-0.973906528517172, -0.148874338981631, -0.865063366688985 }, {
//-0.973906528517172, -0.148874338981631, -0.679409568299024 }, {
//-0.973906528517172, -0.148874338981631, -0.433395394129247 }, {
//-0.973906528517172, -0.148874338981631, -0.148874338981631 }, {
//-0.973906528517172, -0.148874338981631, 0.148874338981631 }, {
// -0.973906528517172, -0.148874338981631, 0.433395394129247 }, {
// -0.973906528517172, -0.148874338981631, 0.679409568299024 }, {
// -0.973906528517172, -0.148874338981631, 0.865063366688985 }, {
// -0.973906528517172, -0.148874338981631, 0.973906528517172 }, {
// -0.973906528517172, 0.148874338981631, -0.973906528517172 }, {
//-0.973906528517172, 0.148874338981631, -0.865063366688985 }, {
//-0.973906528517172, 0.148874338981631, -0.679409568299024 }, {
//-0.973906528517172, 0.148874338981631, -0.433395394129247 }, {
//-0.973906528517172, 0.148874338981631, -0.148874338981631 }, {
//-0.973906528517172, 0.148874338981631, 0.148874338981631 }, {
// -0.973906528517172, 0.148874338981631, 0.433395394129247 }, {
// -0.973906528517172, 0.148874338981631, 0.679409568299024 }, {
// -0.973906528517172, 0.148874338981631, 0.865063366688985 }, {
// -0.973906528517172, 0.148874338981631, 0.973906528517172 }, {
// -0.973906528517172, 0.433395394129247, -0.973906528517172 }, {
//-0.973906528517172, 0.433395394129247, -0.865063366688985 }, {
//-0.973906528517172, 0.433395394129247, -0.679409568299024 }, {
//-0.973906528517172, 0.433395394129247, -0.433395394129247 }, {
//-0.973906528517172, 0.433395394129247, -0.148874338981631 }, {
//-0.973906528517172, 0.433395394129247, 0.148874338981631 }, {
// -0.973906528517172, 0.433395394129247, 0.433395394129247 }, {
// -0.973906528517172, 0.433395394129247, 0.679409568299024 }, {
// -0.973906528517172, 0.433395394129247, 0.865063366688985 }, {
// -0.973906528517172, 0.433395394129247, 0.973906528517172 }, {
// -0.973906528517172, 0.679409568299024, -0.973906528517172 }, {
//-0.973906528517172, 0.679409568299024, -0.865063366688985 }, {
//-0.973906528517172, 0.679409568299024, -0.679409568299024 }, {
//-0.973906528517172, 0.679409568299024, -0.433395394129247 }, {
//-0.973906528517172, 0.679409568299024, -0.148874338981631 }, {
//-0.973906528517172, 0.679409568299024, 0.148874338981631 }, {
// -0.973906528517172, 0.679409568299024, 0.433395394129247 }, {
// -0.973906528517172, 0.679409568299024, 0.679409568299024 }, {
// -0.973906528517172, 0.679409568299024, 0.865063366688985 }, {
// -0.973906528517172, 0.679409568299024, 0.973906528517172 }, {
// -0.973906528517172, 0.865063366688985, -0.973906528517172 }, {
//-0.973906528517172, 0.865063366688985, -0.865063366688985 }, {
//-0.973906528517172, 0.865063366688985, -0.679409568299024 }, {
//-0.973906528517172, 0.865063366688985, -0.433395394129247 }, {
//-0.973906528517172, 0.865063366688985, -0.148874338981631 }, {
//-0.973906528517172, 0.865063366688985, 0.148874338981631 }, {
// -0.973906528517172, 0.865063366688985, 0.433395394129247 }, {
// -0.973906528517172, 0.865063366688985, 0.679409568299024 }, {
// -0.973906528517172, 0.865063366688985, 0.865063366688985 }, {
// -0.973906528517172, 0.865063366688985, 0.973906528517172 }, {
// -0.973906528517172, 0.973906528517172, -0.973906528517172 }, {
//-0.973906528517172, 0.973906528517172, -0.865063366688985 }, {
//-0.973906528517172, 0.973906528517172, -0.679409568299024 }, {
//-0.973906528517172, 0.973906528517172, -0.433395394129247 }, {
//-0.973906528517172, 0.973906528517172, -0.148874338981631 }, {
//-0.973906528517172, 0.973906528517172, 0.148874338981631 }, {
// -0.973906528517172, 0.973906528517172, 0.433395394129247 }, {
// -0.973906528517172, 0.973906528517172, 0.679409568299024 }, {
// -0.973906528517172, 0.973906528517172, 0.865063366688985 }, {
// -0.973906528517172, 0.973906528517172, 0.973906528517172 }, {
// -0.865063366688985, -0.973906528517172, -0.973906528517172 }, {
//-0.865063366688985, -0.973906528517172, -0.865063366688985 }, {
//-0.865063366688985, -0.973906528517172, -0.679409568299024 }, {
//-0.865063366688985, -0.973906528517172, -0.433395394129247 }, {
//-0.865063366688985, -0.973906528517172, -0.148874338981631 }, {
//-0.865063366688985, -0.973906528517172, 0.148874338981631 }, {
// -0.865063366688985, -0.973906528517172, 0.433395394129247 }, {
// -0.865063366688985, -0.973906528517172, 0.679409568299024 }, {
// -0.865063366688985, -0.973906528517172, 0.865063366688985 }, {
// -0.865063366688985, -0.973906528517172, 0.973906528517172 }, {
// -0.865063366688985, -0.865063366688985, -0.973906528517172 }, {
//-0.865063366688985, -0.865063366688985, -0.865063366688985 }, {
//-0.865063366688985, -0.865063366688985, -0.679409568299024 }, {
//-0.865063366688985, -0.865063366688985, -0.433395394129247 }, {
//-0.865063366688985, -0.865063366688985, -0.148874338981631 }, {
//-0.865063366688985, -0.865063366688985, 0.148874338981631 }, {
// -0.865063366688985, -0.865063366688985, 0.433395394129247 }, {
// -0.865063366688985, -0.865063366688985, 0.679409568299024 }, {
// -0.865063366688985, -0.865063366688985, 0.865063366688985 }, {
// -0.865063366688985, -0.865063366688985, 0.973906528517172 }, {
// -0.865063366688985, -0.679409568299024, -0.973906528517172 }, {
//-0.865063366688985, -0.679409568299024, -0.865063366688985 }, {
//-0.865063366688985, -0.679409568299024, -0.679409568299024 }, {
//-0.865063366688985, -0.679409568299024, -0.433395394129247 }, {
//-0.865063366688985, -0.679409568299024, -0.148874338981631 }, {
//-0.865063366688985, -0.679409568299024, 0.148874338981631 }, {
// -0.865063366688985, -0.679409568299024, 0.433395394129247 }, {
// -0.865063366688985, -0.679409568299024, 0.679409568299024 }, {
// -0.865063366688985, -0.679409568299024, 0.865063366688985 }, {
// -0.865063366688985, -0.679409568299024, 0.973906528517172 }, {
// -0.865063366688985, -0.433395394129247, -0.973906528517172 }, {
//-0.865063366688985, -0.433395394129247, -0.865063366688985 }, {
//-0.865063366688985, -0.433395394129247, -0.679409568299024 }, {
//-0.865063366688985, -0.433395394129247, -0.433395394129247 }, {
//-0.865063366688985, -0.433395394129247, -0.148874338981631 }, {
//-0.865063366688985, -0.433395394129247, 0.148874338981631 }, {
// -0.865063366688985, -0.433395394129247, 0.433395394129247 }, {
// -0.865063366688985, -0.433395394129247, 0.679409568299024 }, {
// -0.865063366688985, -0.433395394129247, 0.865063366688985 }, {
// -0.865063366688985, -0.433395394129247, 0.973906528517172 }, {
// -0.865063366688985, -0.148874338981631, -0.973906528517172 }, {
//-0.865063366688985, -0.148874338981631, -0.865063366688985 }, {
//-0.865063366688985, -0.148874338981631, -0.679409568299024 }, {
//-0.865063366688985, -0.148874338981631, -0.433395394129247 }, {
//-0.865063366688985, -0.148874338981631, -0.148874338981631 }, {
//-0.865063366688985, -0.148874338981631, 0.148874338981631 }, {
// -0.865063366688985, -0.148874338981631, 0.433395394129247 }, {
// -0.865063366688985, -0.148874338981631, 0.679409568299024 }, {
// -0.865063366688985, -0.148874338981631, 0.865063366688985 }, {
// -0.865063366688985, -0.148874338981631, 0.973906528517172 }, {
// -0.865063366688985, 0.148874338981631, -0.973906528517172 }, {
//-0.865063366688985, 0.148874338981631, -0.865063366688985 }, {
//-0.865063366688985, 0.148874338981631, -0.679409568299024 }, {
//-0.865063366688985, 0.148874338981631, -0.433395394129247 }, {
//-0.865063366688985, 0.148874338981631, -0.148874338981631 }, {
//-0.865063366688985, 0.148874338981631, 0.148874338981631 }, {
// -0.865063366688985, 0.148874338981631, 0.433395394129247 }, {
// -0.865063366688985, 0.148874338981631, 0.679409568299024 }, {
// -0.865063366688985, 0.148874338981631, 0.865063366688985 }, {
// -0.865063366688985, 0.148874338981631, 0.973906528517172 }, {
// -0.865063366688985, 0.433395394129247, -0.973906528517172 }, {
//-0.865063366688985, 0.433395394129247, -0.865063366688985 }, {
//-0.865063366688985, 0.433395394129247, -0.679409568299024 }, {
//-0.865063366688985, 0.433395394129247, -0.433395394129247 }, {
//-0.865063366688985, 0.433395394129247, -0.148874338981631 }, {
//-0.865063366688985, 0.433395394129247, 0.148874338981631 }, {
// -0.865063366688985, 0.433395394129247, 0.433395394129247 }, {
// -0.865063366688985, 0.433395394129247, 0.679409568299024 }, {
// -0.865063366688985, 0.433395394129247, 0.865063366688985 }, {
// -0.865063366688985, 0.433395394129247, 0.973906528517172 }, {
// -0.865063366688985, 0.679409568299024, -0.973906528517172 }, {
//-0.865063366688985, 0.679409568299024, -0.865063366688985 }, {
//-0.865063366688985, 0.679409568299024, -0.679409568299024 }, {
//-0.865063366688985, 0.679409568299024, -0.433395394129247 }, {
//-0.865063366688985, 0.679409568299024, -0.148874338981631 }, {
//-0.865063366688985, 0.679409568299024, 0.148874338981631 }, {
// -0.865063366688985, 0.679409568299024, 0.433395394129247 }, {
// -0.865063366688985, 0.679409568299024, 0.679409568299024 }, {
// -0.865063366688985, 0.679409568299024, 0.865063366688985 }, {
// -0.865063366688985, 0.679409568299024, 0.973906528517172 }, {
// -0.865063366688985, 0.865063366688985, -0.973906528517172 }, {
//-0.865063366688985, 0.865063366688985, -0.865063366688985 }, {
//-0.865063366688985, 0.865063366688985, -0.679409568299024 }, {
//-0.865063366688985, 0.865063366688985, -0.433395394129247 }, {
//-0.865063366688985, 0.865063366688985, -0.148874338981631 }, {
//-0.865063366688985, 0.865063366688985, 0.148874338981631 }, {
// -0.865063366688985, 0.865063366688985, 0.433395394129247 }, {
// -0.865063366688985, 0.865063366688985, 0.679409568299024 }, {
// -0.865063366688985, 0.865063366688985, 0.865063366688985 }, {
// -0.865063366688985, 0.865063366688985, 0.973906528517172 }, {
// -0.865063366688985, 0.973906528517172, -0.973906528517172 }, {
//-0.865063366688985, 0.973906528517172, -0.865063366688985 }, {
//-0.865063366688985, 0.973906528517172, -0.679409568299024 }, {
//-0.865063366688985, 0.973906528517172, -0.433395394129247 }, {
//-0.865063366688985, 0.973906528517172, -0.148874338981631 }, {
//-0.865063366688985, 0.973906528517172, 0.148874338981631 }, {
// -0.865063366688985, 0.973906528517172, 0.433395394129247 }, {
// -0.865063366688985, 0.973906528517172, 0.679409568299024 }, {
// -0.865063366688985, 0.973906528517172, 0.865063366688985 }, {
// -0.865063366688985, 0.973906528517172, 0.973906528517172 }, {
// -0.679409568299024, -0.973906528517172, -0.973906528517172 }, {
//-0.679409568299024, -0.973906528517172, -0.865063366688985 }, {
//-0.679409568299024, -0.973906528517172, -0.679409568299024 }, {
//-0.679409568299024, -0.973906528517172, -0.433395394129247 }, {
//-0.679409568299024, -0.973906528517172, -0.148874338981631 }, {
//-0.679409568299024, -0.973906528517172, 0.148874338981631 }, {
// -0.679409568299024, -0.973906528517172, 0.433395394129247 }, {
// -0.679409568299024, -0.973906528517172, 0.679409568299024 }, {
// -0.679409568299024, -0.973906528517172, 0.865063366688985 }, {
// -0.679409568299024, -0.973906528517172, 0.973906528517172 }, {
// -0.679409568299024, -0.865063366688985, -0.973906528517172 }, {
//-0.679409568299024, -0.865063366688985, -0.865063366688985 }, {
//-0.679409568299024, -0.865063366688985, -0.679409568299024 }, {
//-0.679409568299024, -0.865063366688985, -0.433395394129247 }, {
//-0.679409568299024, -0.865063366688985, -0.148874338981631 }, {
//-0.679409568299024, -0.865063366688985, 0.148874338981631 }, {
// -0.679409568299024, -0.865063366688985, 0.433395394129247 }, {
// -0.679409568299024, -0.865063366688985, 0.679409568299024 }, {
// -0.679409568299024, -0.865063366688985, 0.865063366688985 }, {
// -0.679409568299024, -0.865063366688985, 0.973906528517172 }, {
// -0.679409568299024, -0.679409568299024, -0.973906528517172 }, {
//-0.679409568299024, -0.679409568299024, -0.865063366688985 }, {
//-0.679409568299024, -0.679409568299024, -0.679409568299024 }, {
//-0.679409568299024, -0.679409568299024, -0.433395394129247 }, {
//-0.679409568299024, -0.679409568299024, -0.148874338981631 }, {
//-0.679409568299024, -0.679409568299024, 0.148874338981631 }, {
// -0.679409568299024, -0.679409568299024, 0.433395394129247 }, {
// -0.679409568299024, -0.679409568299024, 0.679409568299024 }, {
// -0.679409568299024, -0.679409568299024, 0.865063366688985 }, {
// -0.679409568299024, -0.679409568299024, 0.973906528517172 }, {
// -0.679409568299024, -0.433395394129247, -0.973906528517172 }, {
//-0.679409568299024, -0.433395394129247, -0.865063366688985 }, {
//-0.679409568299024, -0.433395394129247, -0.679409568299024 }, {
//-0.679409568299024, -0.433395394129247, -0.433395394129247 }, {
//-0.679409568299024, -0.433395394129247, -0.148874338981631 }, {
//-0.679409568299024, -0.433395394129247, 0.148874338981631 }, {
// -0.679409568299024, -0.433395394129247, 0.433395394129247 }, {
// -0.679409568299024, -0.433395394129247, 0.679409568299024 }, {
// -0.679409568299024, -0.433395394129247, 0.865063366688985 }, {
// -0.679409568299024, -0.433395394129247, 0.973906528517172 }, {
// -0.679409568299024, -0.148874338981631, -0.973906528517172 }, {
//-0.679409568299024, -0.148874338981631, -0.865063366688985 }, {
//-0.679409568299024, -0.148874338981631, -0.679409568299024 }, {
//-0.679409568299024, -0.148874338981631, -0.433395394129247 }, {
//-0.679409568299024, -0.148874338981631, -0.148874338981631 }, {
//-0.679409568299024, -0.148874338981631, 0.148874338981631 }, {
// -0.679409568299024, -0.148874338981631, 0.433395394129247 }, {
// -0.679409568299024, -0.148874338981631, 0.679409568299024 }, {
// -0.679409568299024, -0.148874338981631, 0.865063366688985 }, {
// -0.679409568299024, -0.148874338981631, 0.973906528517172 }, {
// -0.679409568299024, 0.148874338981631, -0.973906528517172 }, {
//-0.679409568299024, 0.148874338981631, -0.865063366688985 }, {
//-0.679409568299024, 0.148874338981631, -0.679409568299024 }, {
//-0.679409568299024, 0.148874338981631, -0.433395394129247 }, {
//-0.679409568299024, 0.148874338981631, -0.148874338981631 }, {
//-0.679409568299024, 0.148874338981631, 0.148874338981631 }, {
// -0.679409568299024, 0.148874338981631, 0.433395394129247 }, {
// -0.679409568299024, 0.148874338981631, 0.679409568299024 }, {
// -0.679409568299024, 0.148874338981631, 0.865063366688985 }, {
// -0.679409568299024, 0.148874338981631, 0.973906528517172 }, {
// -0.679409568299024, 0.433395394129247, -0.973906528517172 }, {
//-0.679409568299024, 0.433395394129247, -0.865063366688985 }, {
//-0.679409568299024, 0.433395394129247, -0.679409568299024 }, {
//-0.679409568299024, 0.433395394129247, -0.433395394129247 }, {
//-0.679409568299024, 0.433395394129247, -0.148874338981631 }, {
//-0.679409568299024, 0.433395394129247, 0.148874338981631 }, {
// -0.679409568299024, 0.433395394129247, 0.433395394129247 }, {
// -0.679409568299024, 0.433395394129247, 0.679409568299024 }, {
// -0.679409568299024, 0.433395394129247, 0.865063366688985 }, {
// -0.679409568299024, 0.433395394129247, 0.973906528517172 }, {
// -0.679409568299024, 0.679409568299024, -0.973906528517172 }, {
//-0.679409568299024, 0.679409568299024, -0.865063366688985 }, {
//-0.679409568299024, 0.679409568299024, -0.679409568299024 }, {
//-0.679409568299024, 0.679409568299024, -0.433395394129247 }, {
//-0.679409568299024, 0.679409568299024, -0.148874338981631 }, {
//-0.679409568299024, 0.679409568299024, 0.148874338981631 }, {
// -0.679409568299024, 0.679409568299024, 0.433395394129247 }, {
// -0.679409568299024, 0.679409568299024, 0.679409568299024 }, {
// -0.679409568299024, 0.679409568299024, 0.865063366688985 }, {
// -0.679409568299024, 0.679409568299024, 0.973906528517172 }, {
// -0.679409568299024, 0.865063366688985, -0.973906528517172 }, {
//-0.679409568299024, 0.865063366688985, -0.865063366688985 }, {
//-0.679409568299024, 0.865063366688985, -0.679409568299024 }, {
//-0.679409568299024, 0.865063366688985, -0.433395394129247 }, {
//-0.679409568299024, 0.865063366688985, -0.148874338981631 }, {
//-0.679409568299024, 0.865063366688985, 0.148874338981631 }, {
// -0.679409568299024, 0.865063366688985, 0.433395394129247 }, {
// -0.679409568299024, 0.865063366688985, 0.679409568299024 }, {
// -0.679409568299024, 0.865063366688985, 0.865063366688985 }, {
// -0.679409568299024, 0.865063366688985, 0.973906528517172 }, {
// -0.679409568299024, 0.973906528517172, -0.973906528517172 }, {
//-0.679409568299024, 0.973906528517172, -0.865063366688985 }, {
//-0.679409568299024, 0.973906528517172, -0.679409568299024 }, {
//-0.679409568299024, 0.973906528517172, -0.433395394129247 }, {
//-0.679409568299024, 0.973906528517172, -0.148874338981631 }, {
//-0.679409568299024, 0.973906528517172, 0.148874338981631 }, {
// -0.679409568299024, 0.973906528517172, 0.433395394129247 }, {
// -0.679409568299024, 0.973906528517172, 0.679409568299024 }, {
// -0.679409568299024, 0.973906528517172, 0.865063366688985 }, {
// -0.679409568299024, 0.973906528517172, 0.973906528517172 }, {
// -0.433395394129247, -0.973906528517172, -0.973906528517172 }, {
//-0.433395394129247, -0.973906528517172, -0.865063366688985 }, {
//-0.433395394129247, -0.973906528517172, -0.679409568299024 }, {
//-0.433395394129247, -0.973906528517172, -0.433395394129247 }, {
//-0.433395394129247, -0.973906528517172, -0.148874338981631 }, {
//-0.433395394129247, -0.973906528517172, 0.148874338981631 }, {
// -0.433395394129247, -0.973906528517172, 0.433395394129247 }, {
// -0.433395394129247, -0.973906528517172, 0.679409568299024 }, {
// -0.433395394129247, -0.973906528517172, 0.865063366688985 }, {
// -0.433395394129247, -0.973906528517172, 0.973906528517172 }, {
// -0.433395394129247, -0.865063366688985, -0.973906528517172 }, {
//-0.433395394129247, -0.865063366688985, -0.865063366688985 }, {
//-0.433395394129247, -0.865063366688985, -0.679409568299024 }, {
//-0.433395394129247, -0.865063366688985, -0.433395394129247 }, {
//-0.433395394129247, -0.865063366688985, -0.148874338981631 }, {
//-0.433395394129247, -0.865063366688985, 0.148874338981631 }, {
// -0.433395394129247, -0.865063366688985, 0.433395394129247 }, {
// -0.433395394129247, -0.865063366688985, 0.679409568299024 }, {
// -0.433395394129247, -0.865063366688985, 0.865063366688985 }, {
// -0.433395394129247, -0.865063366688985, 0.973906528517172 }, {
// -0.433395394129247, -0.679409568299024, -0.973906528517172 }, {
//-0.433395394129247, -0.679409568299024, -0.865063366688985 }, {
//-0.433395394129247, -0.679409568299024, -0.679409568299024 }, {
//-0.433395394129247, -0.679409568299024, -0.433395394129247 }, {
//-0.433395394129247, -0.679409568299024, -0.148874338981631 }, {
//-0.433395394129247, -0.679409568299024, 0.148874338981631 }, {
// -0.433395394129247, -0.679409568299024, 0.433395394129247 }, {
// -0.433395394129247, -0.679409568299024, 0.679409568299024 }, {
// -0.433395394129247, -0.679409568299024, 0.865063366688985 }, {
// -0.433395394129247, -0.679409568299024, 0.973906528517172 }, {
// -0.433395394129247, -0.433395394129247, -0.973906528517172 }, {
//-0.433395394129247, -0.433395394129247, -0.865063366688985 }, {
//-0.433395394129247, -0.433395394129247, -0.679409568299024 }, {
//-0.433395394129247, -0.433395394129247, -0.433395394129247 }, {
//-0.433395394129247, -0.433395394129247, -0.148874338981631 }, {
//-0.433395394129247, -0.433395394129247, 0.148874338981631 }, {
// -0.433395394129247, -0.433395394129247, 0.433395394129247 }, {
// -0.433395394129247, -0.433395394129247, 0.679409568299024 }, {
// -0.433395394129247, -0.433395394129247, 0.865063366688985 }, {
// -0.433395394129247, -0.433395394129247, 0.973906528517172 }, {
// -0.433395394129247, -0.148874338981631, -0.973906528517172 }, {
//-0.433395394129247, -0.148874338981631, -0.865063366688985 }, {
//-0.433395394129247, -0.148874338981631, -0.679409568299024 }, {
//-0.433395394129247, -0.148874338981631, -0.433395394129247 }, {
//-0.433395394129247, -0.148874338981631, -0.148874338981631 }, {
//-0.433395394129247, -0.148874338981631, 0.148874338981631 }, {
// -0.433395394129247, -0.148874338981631, 0.433395394129247 }, {
// -0.433395394129247, -0.148874338981631, 0.679409568299024 }, {
// -0.433395394129247, -0.148874338981631, 0.865063366688985 }, {
// -0.433395394129247, -0.148874338981631, 0.973906528517172 }, {
// -0.433395394129247, 0.148874338981631, -0.973906528517172 }, {
//-0.433395394129247, 0.148874338981631, -0.865063366688985 }, {
//-0.433395394129247, 0.148874338981631, -0.679409568299024 }, {
//-0.433395394129247, 0.148874338981631, -0.433395394129247 }, {
//-0.433395394129247, 0.148874338981631, -0.148874338981631 }, {
//-0.433395394129247, 0.148874338981631, 0.148874338981631 }, {
// -0.433395394129247, 0.148874338981631, 0.433395394129247 }, {
// -0.433395394129247, 0.148874338981631, 0.679409568299024 }, {
// -0.433395394129247, 0.148874338981631, 0.865063366688985 }, {
// -0.433395394129247, 0.148874338981631, 0.973906528517172 }, {
// -0.433395394129247, 0.433395394129247, -0.973906528517172 }, {
//-0.433395394129247, 0.433395394129247, -0.865063366688985 }, {
//-0.433395394129247, 0.433395394129247, -0.679409568299024 }, {
//-0.433395394129247, 0.433395394129247, -0.433395394129247 }, {
//-0.433395394129247, 0.433395394129247, -0.148874338981631 }, {
//-0.433395394129247, 0.433395394129247, 0.148874338981631 }, {
// -0.433395394129247, 0.433395394129247, 0.433395394129247 }, {
// -0.433395394129247, 0.433395394129247, 0.679409568299024 }, {
// -0.433395394129247, 0.433395394129247, 0.865063366688985 }, {
// -0.433395394129247, 0.433395394129247, 0.973906528517172 }, {
// -0.433395394129247, 0.679409568299024, -0.973906528517172 }, {
//-0.433395394129247, 0.679409568299024, -0.865063366688985 }, {
//-0.433395394129247, 0.679409568299024, -0.679409568299024 }, {
//-0.433395394129247, 0.679409568299024, -0.433395394129247 }, {
//-0.433395394129247, 0.679409568299024, -0.148874338981631 }, {
//-0.433395394129247, 0.679409568299024, 0.148874338981631 }, {
// -0.433395394129247, 0.679409568299024, 0.433395394129247 }, {
// -0.433395394129247, 0.679409568299024, 0.679409568299024 }, {
// -0.433395394129247, 0.679409568299024, 0.865063366688985 }, {
// -0.433395394129247, 0.679409568299024, 0.973906528517172 }, {
// -0.433395394129247, 0.865063366688985, -0.973906528517172 }, {
//-0.433395394129247, 0.865063366688985, -0.865063366688985 }, {
//-0.433395394129247, 0.865063366688985, -0.679409568299024 }, {
//-0.433395394129247, 0.865063366688985, -0.433395394129247 }, {
//-0.433395394129247, 0.865063366688985, -0.148874338981631 }, {
//-0.433395394129247, 0.865063366688985, 0.148874338981631 }, {
// -0.433395394129247, 0.865063366688985, 0.433395394129247 }, {
// -0.433395394129247, 0.865063366688985, 0.679409568299024 }, {
// -0.433395394129247, 0.865063366688985, 0.865063366688985 }, {
// -0.433395394129247, 0.865063366688985, 0.973906528517172 }, {
// -0.433395394129247, 0.973906528517172, -0.973906528517172 }, {
//-0.433395394129247, 0.973906528517172, -0.865063366688985 }, {
//-0.433395394129247, 0.973906528517172, -0.679409568299024 }, {
//-0.433395394129247, 0.973906528517172, -0.433395394129247 }, {
//-0.433395394129247, 0.973906528517172, -0.148874338981631 }, {
//-0.433395394129247, 0.973906528517172, 0.148874338981631 }, {
// -0.433395394129247, 0.973906528517172, 0.433395394129247 }, {
// -0.433395394129247, 0.973906528517172, 0.679409568299024 }, {
// -0.433395394129247, 0.973906528517172, 0.865063366688985 }, {
// -0.433395394129247, 0.973906528517172, 0.973906528517172 }, {
// -0.148874338981631, -0.973906528517172, -0.973906528517172 }, {
//-0.148874338981631, -0.973906528517172, -0.865063366688985 }, {
//-0.148874338981631, -0.973906528517172, -0.679409568299024 }, {
//-0.148874338981631, -0.973906528517172, -0.433395394129247 }, {
//-0.148874338981631, -0.973906528517172, -0.148874338981631 }, {
//-0.148874338981631, -0.973906528517172, 0.148874338981631 }, {
// -0.148874338981631, -0.973906528517172, 0.433395394129247 }, {
// -0.148874338981631, -0.973906528517172, 0.679409568299024 }, {
// -0.148874338981631, -0.973906528517172, 0.865063366688985 }, {
// -0.148874338981631, -0.973906528517172, 0.973906528517172 }, {
// -0.148874338981631, -0.865063366688985, -0.973906528517172 }, {
//-0.148874338981631, -0.865063366688985, -0.865063366688985 }, {
//-0.148874338981631, -0.865063366688985, -0.679409568299024 }, {
//-0.148874338981631, -0.865063366688985, -0.433395394129247 }, {
//-0.148874338981631, -0.865063366688985, -0.148874338981631 }, {
//-0.148874338981631, -0.865063366688985, 0.148874338981631 }, {
// -0.148874338981631, -0.865063366688985, 0.433395394129247 }, {
// -0.148874338981631, -0.865063366688985, 0.679409568299024 }, {
// -0.148874338981631, -0.865063366688985, 0.865063366688985 }, {
// -0.148874338981631, -0.865063366688985, 0.973906528517172 }, {
// -0.148874338981631, -0.679409568299024, -0.973906528517172 }, {
//-0.148874338981631, -0.679409568299024, -0.865063366688985 }, {
//-0.148874338981631, -0.679409568299024, -0.679409568299024 }, {
//-0.148874338981631, -0.679409568299024, -0.433395394129247 }, {
//-0.148874338981631, -0.679409568299024, -0.148874338981631 }, {
//-0.148874338981631, -0.679409568299024, 0.148874338981631 }, {
// -0.148874338981631, -0.679409568299024, 0.433395394129247 }, {
// -0.148874338981631, -0.679409568299024, 0.679409568299024 }, {
// -0.148874338981631, -0.679409568299024, 0.865063366688985 }, {
// -0.148874338981631, -0.679409568299024, 0.973906528517172 }, {
// -0.148874338981631, -0.433395394129247, -0.973906528517172 }, {
//-0.148874338981631, -0.433395394129247, -0.865063366688985 }, {
//-0.148874338981631, -0.433395394129247, -0.679409568299024 }, {
//-0.148874338981631, -0.433395394129247, -0.433395394129247 }, {
//-0.148874338981631, -0.433395394129247, -0.148874338981631 }, {
//-0.148874338981631, -0.433395394129247, 0.148874338981631 }, {
// -0.148874338981631, -0.433395394129247, 0.433395394129247 }, {
// -0.148874338981631, -0.433395394129247, 0.679409568299024 }, {
// -0.148874338981631, -0.433395394129247, 0.865063366688985 }, {
// -0.148874338981631, -0.433395394129247, 0.973906528517172 }, {
// -0.148874338981631, -0.148874338981631, -0.973906528517172 }, {
//-0.148874338981631, -0.148874338981631, -0.865063366688985 }, {
//-0.148874338981631, -0.148874338981631, -0.679409568299024 }, {
//-0.148874338981631, -0.148874338981631, -0.433395394129247 }, {
//-0.148874338981631, -0.148874338981631, -0.148874338981631 }, {
//-0.148874338981631, -0.148874338981631, 0.148874338981631 }, {
// -0.148874338981631, -0.148874338981631, 0.433395394129247 }, {
// -0.148874338981631, -0.148874338981631, 0.679409568299024 }, {
// -0.148874338981631, -0.148874338981631, 0.865063366688985 }, {
// -0.148874338981631, -0.148874338981631, 0.973906528517172 }, {
// -0.148874338981631, 0.148874338981631, -0.973906528517172 }, {
//-0.148874338981631, 0.148874338981631, -0.865063366688985 }, {
//-0.148874338981631, 0.148874338981631, -0.679409568299024 }, {
//-0.148874338981631, 0.148874338981631, -0.433395394129247 }, {
//-0.148874338981631, 0.148874338981631, -0.148874338981631 }, {
//-0.148874338981631, 0.148874338981631, 0.148874338981631 }, {
// -0.148874338981631, 0.148874338981631, 0.433395394129247 }, {
// -0.148874338981631, 0.148874338981631, 0.679409568299024 }, {
// -0.148874338981631, 0.148874338981631, 0.865063366688985 }, {
// -0.148874338981631, 0.148874338981631, 0.973906528517172 }, {
// -0.148874338981631, 0.433395394129247, -0.973906528517172 }, {
//-0.148874338981631, 0.433395394129247, -0.865063366688985 }, {
//-0.148874338981631, 0.433395394129247, -0.679409568299024 }, {
//-0.148874338981631, 0.433395394129247, -0.433395394129247 }, {
//-0.148874338981631, 0.433395394129247, -0.148874338981631 }, {
//-0.148874338981631, 0.433395394129247, 0.148874338981631 }, {
// -0.148874338981631, 0.433395394129247, 0.433395394129247 }, {
// -0.148874338981631, 0.433395394129247, 0.679409568299024 }, {
// -0.148874338981631, 0.433395394129247, 0.865063366688985 }, {
// -0.148874338981631, 0.433395394129247, 0.973906528517172 }, {
// -0.148874338981631, 0.679409568299024, -0.973906528517172 }, {
//-0.148874338981631, 0.679409568299024, -0.865063366688985 }, {
//-0.148874338981631, 0.679409568299024, -0.679409568299024 }, {
//-0.148874338981631, 0.679409568299024, -0.433395394129247 }, {
//-0.148874338981631, 0.679409568299024, -0.148874338981631 }, {
//-0.148874338981631, 0.679409568299024, 0.148874338981631 }, {
// -0.148874338981631, 0.679409568299024, 0.433395394129247 }, {
// -0.148874338981631, 0.679409568299024, 0.679409568299024 }, {
// -0.148874338981631, 0.679409568299024, 0.865063366688985 }, {
// -0.148874338981631, 0.679409568299024, 0.973906528517172 }, {
// -0.148874338981631, 0.865063366688985, -0.973906528517172 }, {
//-0.148874338981631, 0.865063366688985, -0.865063366688985 }, {
//-0.148874338981631, 0.865063366688985, -0.679409568299024 }, {
//-0.148874338981631, 0.865063366688985, -0.433395394129247 }, {
//-0.148874338981631, 0.865063366688985, -0.148874338981631 }, {
//-0.148874338981631, 0.865063366688985, 0.148874338981631 }, {
// -0.148874338981631, 0.865063366688985, 0.433395394129247 }, {
// -0.148874338981631, 0.865063366688985, 0.679409568299024 }, {
// -0.148874338981631, 0.865063366688985, 0.865063366688985 }, {
// -0.148874338981631, 0.865063366688985, 0.973906528517172 }, {
// -0.148874338981631, 0.973906528517172, -0.973906528517172 }, {
//-0.148874338981631, 0.973906528517172, -0.865063366688985 }, {
//-0.148874338981631, 0.973906528517172, -0.679409568299024 }, {
//-0.148874338981631, 0.973906528517172, -0.433395394129247 }, {
//-0.148874338981631, 0.973906528517172, -0.148874338981631 }, {
//-0.148874338981631, 0.973906528517172, 0.148874338981631 }, {
// -0.148874338981631, 0.973906528517172, 0.433395394129247 }, {
// -0.148874338981631, 0.973906528517172, 0.679409568299024 }, {
// -0.148874338981631, 0.973906528517172, 0.865063366688985 }, {
// -0.148874338981631, 0.973906528517172, 0.973906528517172 }, {
// 0.148874338981631, -0.973906528517172, -0.973906528517172 }, {
// 0.148874338981631, -0.973906528517172, -0.865063366688985 }, {
// 0.148874338981631, -0.973906528517172, -0.679409568299024 }, {
// 0.148874338981631, -0.973906528517172, -0.433395394129247 }, {
// 0.148874338981631, -0.973906528517172, -0.148874338981631 }, {
// 0.148874338981631, -0.973906528517172, 0.148874338981631 }, {
// 0.148874338981631, -0.973906528517172, 0.433395394129247 }, {
// 0.148874338981631, -0.973906528517172, 0.679409568299024 }, {
// 0.148874338981631, -0.973906528517172, 0.865063366688985 }, {
// 0.148874338981631, -0.973906528517172, 0.973906528517172 }, {
// 0.148874338981631, -0.865063366688985, -0.973906528517172 }, {
// 0.148874338981631, -0.865063366688985, -0.865063366688985 }, {
// 0.148874338981631, -0.865063366688985, -0.679409568299024 }, {
// 0.148874338981631, -0.865063366688985, -0.433395394129247 }, {
// 0.148874338981631, -0.865063366688985, -0.148874338981631 }, {
// 0.148874338981631, -0.865063366688985, 0.148874338981631 }, {
// 0.148874338981631, -0.865063366688985, 0.433395394129247 }, {
// 0.148874338981631, -0.865063366688985, 0.679409568299024 }, {
// 0.148874338981631, -0.865063366688985, 0.865063366688985 }, {
// 0.148874338981631, -0.865063366688985, 0.973906528517172 }, {
// 0.148874338981631, -0.679409568299024, -0.973906528517172 }, {
// 0.148874338981631, -0.679409568299024, -0.865063366688985 }, {
// 0.148874338981631, -0.679409568299024, -0.679409568299024 }, {
// 0.148874338981631, -0.679409568299024, -0.433395394129247 }, {
// 0.148874338981631, -0.679409568299024, -0.148874338981631 }, {
// 0.148874338981631, -0.679409568299024, 0.148874338981631 }, {
// 0.148874338981631, -0.679409568299024, 0.433395394129247 }, {
// 0.148874338981631, -0.679409568299024, 0.679409568299024 }, {
// 0.148874338981631, -0.679409568299024, 0.865063366688985 }, {
// 0.148874338981631, -0.679409568299024, 0.973906528517172 }, {
// 0.148874338981631, -0.433395394129247, -0.973906528517172 }, {
// 0.148874338981631, -0.433395394129247, -0.865063366688985 }, {
// 0.148874338981631, -0.433395394129247, -0.679409568299024 }, {
// 0.148874338981631, -0.433395394129247, -0.433395394129247 }, {
// 0.148874338981631, -0.433395394129247, -0.148874338981631 }, {
// 0.148874338981631, -0.433395394129247, 0.148874338981631 }, {
// 0.148874338981631, -0.433395394129247, 0.433395394129247 }, {
// 0.148874338981631, -0.433395394129247, 0.679409568299024 }, {
// 0.148874338981631, -0.433395394129247, 0.865063366688985 }, {
// 0.148874338981631, -0.433395394129247, 0.973906528517172 }, {
// 0.148874338981631, -0.148874338981631, -0.973906528517172 }, {
// 0.148874338981631, -0.148874338981631, -0.865063366688985 }, {
// 0.148874338981631, -0.148874338981631, -0.679409568299024 }, {
// 0.148874338981631, -0.148874338981631, -0.433395394129247 }, {
// 0.148874338981631, -0.148874338981631, -0.148874338981631 }, {
// 0.148874338981631, -0.148874338981631, 0.148874338981631 }, {
// 0.148874338981631, -0.148874338981631, 0.433395394129247 }, {
// 0.148874338981631, -0.148874338981631, 0.679409568299024 }, {
// 0.148874338981631, -0.148874338981631, 0.865063366688985 }, {
// 0.148874338981631, -0.148874338981631, 0.973906528517172 }, {
// 0.148874338981631, 0.148874338981631, -0.973906528517172 }, {
// 0.148874338981631, 0.148874338981631, -0.865063366688985 }, {
// 0.148874338981631, 0.148874338981631, -0.679409568299024 }, {
// 0.148874338981631, 0.148874338981631, -0.433395394129247 }, {
// 0.148874338981631, 0.148874338981631, -0.148874338981631 }, {
// 0.148874338981631, 0.148874338981631, 0.148874338981631 }, {
// 0.148874338981631, 0.148874338981631, 0.433395394129247
//}, { 0.148874338981631, 0.148874338981631, 0.679409568299024 }, {
// 0.148874338981631, 0.148874338981631, 0.865063366688985 }, {
// 0.148874338981631, 0.148874338981631, 0.973906528517172 }, {
// 0.148874338981631, 0.433395394129247, -0.973906528517172 }, {
// 0.148874338981631, 0.433395394129247, -0.865063366688985 }, {
// 0.148874338981631, 0.433395394129247, -0.679409568299024 }, {
// 0.148874338981631, 0.433395394129247, -0.433395394129247 }, {
// 0.148874338981631, 0.433395394129247, -0.148874338981631 }, {
// 0.148874338981631, 0.433395394129247, 0.148874338981631 }, {
// 0.148874338981631, 0.433395394129247, 0.433395394129247 }, {
// 0.148874338981631, 0.433395394129247, 0.679409568299024 }, {
// 0.148874338981631, 0.433395394129247, 0.865063366688985 }, {
// 0.148874338981631, 0.433395394129247, 0.973906528517172 }, {
// 0.148874338981631, 0.679409568299024, -0.973906528517172 }, {
// 0.148874338981631, 0.679409568299024, -0.865063366688985 }, {
// 0.148874338981631, 0.679409568299024, -0.679409568299024 }, {
// 0.148874338981631, 0.679409568299024, -0.433395394129247 }, {
// 0.148874338981631, 0.679409568299024, -0.148874338981631 }, {
// 0.148874338981631, 0.679409568299024, 0.148874338981631 }, {
// 0.148874338981631, 0.679409568299024, 0.433395394129247 }, {
// 0.148874338981631, 0.679409568299024, 0.679409568299024 }, {
// 0.148874338981631, 0.679409568299024, 0.865063366688985 }, {
// 0.148874338981631, 0.679409568299024, 0.973906528517172 }, {
// 0.148874338981631, 0.865063366688985, -0.973906528517172 }, {
// 0.148874338981631, 0.865063366688985, -0.865063366688985 }, {
// 0.148874338981631, 0.865063366688985, -0.679409568299024 }, {
// 0.148874338981631, 0.865063366688985, -0.433395394129247 }, {
// 0.148874338981631, 0.865063366688985, -0.148874338981631 }, {
// 0.148874338981631, 0.865063366688985, 0.148874338981631 }, {
// 0.148874338981631, 0.865063366688985, 0.433395394129247 }, {
// 0.148874338981631, 0.865063366688985, 0.679409568299024 }, {
// 0.148874338981631, 0.865063366688985, 0.865063366688985 }, {
// 0.148874338981631, 0.865063366688985, 0.973906528517172 }, {
// 0.148874338981631, 0.973906528517172, -0.973906528517172 }, {
// 0.148874338981631, 0.973906528517172, -0.865063366688985 }, {
// 0.148874338981631, 0.973906528517172, -0.679409568299024 }, {
// 0.148874338981631, 0.973906528517172, -0.433395394129247 }, {
// 0.148874338981631, 0.973906528517172, -0.148874338981631 }, {
// 0.148874338981631, 0.973906528517172, 0.148874338981631 }, {
// 0.148874338981631, 0.973906528517172, 0.433395394129247 }, {
// 0.148874338981631, 0.973906528517172, 0.679409568299024 }, {
// 0.148874338981631, 0.973906528517172, 0.865063366688985 }, {
// 0.148874338981631, 0.973906528517172, 0.973906528517172 }, {
// 0.433395394129247, -0.973906528517172, -0.973906528517172 }, {
// 0.433395394129247, -0.973906528517172, -0.865063366688985 }, {
// 0.433395394129247, -0.973906528517172, -0.679409568299024 }, {
// 0.433395394129247, -0.973906528517172, -0.433395394129247 }, {
// 0.433395394129247, -0.973906528517172, -0.148874338981631 }, {
// 0.433395394129247, -0.973906528517172, 0.148874338981631 }, {
// 0.433395394129247, -0.973906528517172, 0.433395394129247 }, {
// 0.433395394129247, -0.973906528517172, 0.679409568299024 }, {
// 0.433395394129247, -0.973906528517172, 0.865063366688985 }, {
// 0.433395394129247, -0.973906528517172, 0.973906528517172 }, {
// 0.433395394129247, -0.865063366688985, -0.973906528517172 }, {
// 0.433395394129247, -0.865063366688985, -0.865063366688985 }, {
// 0.433395394129247, -0.865063366688985, -0.679409568299024 }, {
// 0.433395394129247, -0.865063366688985, -0.433395394129247 }, {
// 0.433395394129247, -0.865063366688985, -0.148874338981631 }, {
// 0.433395394129247, -0.865063366688985, 0.148874338981631 }, {
// 0.433395394129247, -0.865063366688985, 0.433395394129247 }, {
// 0.433395394129247, -0.865063366688985, 0.679409568299024 }, {
// 0.433395394129247, -0.865063366688985, 0.865063366688985 }, {
// 0.433395394129247, -0.865063366688985, 0.973906528517172 }, {
// 0.433395394129247, -0.679409568299024, -0.973906528517172 }, {
// 0.433395394129247, -0.679409568299024, -0.865063366688985 }, {
// 0.433395394129247, -0.679409568299024, -0.679409568299024 }, {
// 0.433395394129247, -0.679409568299024, -0.433395394129247 }, {
// 0.433395394129247, -0.679409568299024, -0.148874338981631 }, {
// 0.433395394129247, -0.679409568299024, 0.148874338981631 }, {
// 0.433395394129247, -0.679409568299024, 0.433395394129247 }, {
// 0.433395394129247, -0.679409568299024, 0.679409568299024 }, {
// 0.433395394129247, -0.679409568299024, 0.865063366688985 }, {
// 0.433395394129247, -0.679409568299024, 0.973906528517172 }, {
// 0.433395394129247, -0.433395394129247, -0.973906528517172 }, {
// 0.433395394129247, -0.433395394129247, -0.865063366688985 }, {
// 0.433395394129247, -0.433395394129247, -0.679409568299024 }, {
// 0.433395394129247, -0.433395394129247, -0.433395394129247 }, {
// 0.433395394129247, -0.433395394129247, -0.148874338981631 }, {
// 0.433395394129247, -0.433395394129247, 0.148874338981631 }, {
// 0.433395394129247, -0.433395394129247, 0.433395394129247 }, {
// 0.433395394129247, -0.433395394129247, 0.679409568299024 }, {
// 0.433395394129247, -0.433395394129247, 0.865063366688985 }, {
// 0.433395394129247, -0.433395394129247, 0.973906528517172 }, {
// 0.433395394129247, -0.148874338981631, -0.973906528517172 }, {
// 0.433395394129247, -0.148874338981631, -0.865063366688985 }, {
// 0.433395394129247, -0.148874338981631, -0.679409568299024 }, {
// 0.433395394129247, -0.148874338981631, -0.433395394129247 }, {
// 0.433395394129247, -0.148874338981631, -0.148874338981631 }, {
// 0.433395394129247, -0.148874338981631, 0.148874338981631 }, {
// 0.433395394129247, -0.148874338981631, 0.433395394129247 }, {
// 0.433395394129247, -0.148874338981631, 0.679409568299024 }, {
// 0.433395394129247, -0.148874338981631, 0.865063366688985 }, {
// 0.433395394129247, -0.148874338981631, 0.973906528517172 }, {
// 0.433395394129247, 0.148874338981631, -0.973906528517172 }, {
// 0.433395394129247, 0.148874338981631, -0.865063366688985 }, {
// 0.433395394129247, 0.148874338981631, -0.679409568299024 }, {
// 0.433395394129247, 0.148874338981631, -0.433395394129247 }, {
// 0.433395394129247, 0.148874338981631, -0.148874338981631 }, {
// 0.433395394129247, 0.148874338981631, 0.148874338981631 }, {
// 0.433395394129247, 0.148874338981631, 0.433395394129247 }, {
// 0.433395394129247, 0.148874338981631, 0.679409568299024 }, {
// 0.433395394129247, 0.148874338981631, 0.865063366688985 }, {
// 0.433395394129247, 0.148874338981631, 0.973906528517172 }, {
// 0.433395394129247, 0.433395394129247, -0.973906528517172 }, {
// 0.433395394129247, 0.433395394129247, -0.865063366688985 }, {
// 0.433395394129247, 0.433395394129247, -0.679409568299024 }, {
// 0.433395394129247, 0.433395394129247, -0.433395394129247 }, {
// 0.433395394129247, 0.433395394129247, -0.148874338981631 }, {
// 0.433395394129247, 0.433395394129247, 0.148874338981631 }, {
// 0.433395394129247, 0.433395394129247, 0.433395394129247 }, {
// 0.433395394129247, 0.433395394129247, 0.679409568299024 }, {
// 0.433395394129247, 0.433395394129247, 0.865063366688985 }, {
// 0.433395394129247, 0.433395394129247, 0.973906528517172 }, {
// 0.433395394129247, 0.679409568299024, -0.973906528517172 }, {
// 0.433395394129247, 0.679409568299024, -0.865063366688985 }, {
// 0.433395394129247, 0.679409568299024, -0.679409568299024 }, {
// 0.433395394129247, 0.679409568299024, -0.433395394129247 }, {
// 0.433395394129247, 0.679409568299024, -0.148874338981631 }, {
// 0.433395394129247, 0.679409568299024, 0.148874338981631 }, {
// 0.433395394129247, 0.679409568299024, 0.433395394129247 }, {
// 0.433395394129247, 0.679409568299024, 0.679409568299024 }, {
// 0.433395394129247, 0.679409568299024, 0.865063366688985 }, {
// 0.433395394129247, 0.679409568299024, 0.973906528517172 }, {
// 0.433395394129247, 0.865063366688985, -0.973906528517172 }, {
// 0.433395394129247, 0.865063366688985, -0.865063366688985 }, {
// 0.433395394129247, 0.865063366688985, -0.679409568299024 }, {
// 0.433395394129247, 0.865063366688985, -0.433395394129247 }, {
// 0.433395394129247, 0.865063366688985, -0.148874338981631 }, {
// 0.433395394129247, 0.865063366688985, 0.148874338981631 }, {
// 0.433395394129247, 0.865063366688985, 0.433395394129247 }, {
// 0.433395394129247, 0.865063366688985, 0.679409568299024 }, {
// 0.433395394129247, 0.865063366688985, 0.865063366688985 }, {
// 0.433395394129247, 0.865063366688985, 0.973906528517172 }, {
// 0.433395394129247, 0.973906528517172, -0.973906528517172 }, {
// 0.433395394129247, 0.973906528517172, -0.865063366688985 }, {
// 0.433395394129247, 0.973906528517172, -0.679409568299024 }, {
// 0.433395394129247, 0.973906528517172, -0.433395394129247 }, {
// 0.433395394129247, 0.973906528517172, -0.148874338981631 }, {
// 0.433395394129247, 0.973906528517172, 0.148874338981631 }, {
// 0.433395394129247, 0.973906528517172, 0.433395394129247 }, {
// 0.433395394129247, 0.973906528517172, 0.679409568299024 }, {
// 0.433395394129247, 0.973906528517172, 0.865063366688985 }, {
// 0.433395394129247, 0.973906528517172, 0.973906528517172 }, {
// 0.679409568299024, -0.973906528517172, -0.973906528517172 }, {
// 0.679409568299024, -0.973906528517172, -0.865063366688985 }, {
// 0.679409568299024, -0.973906528517172, -0.679409568299024 }, {
// 0.679409568299024, -0.973906528517172, -0.433395394129247 }, {
// 0.679409568299024, -0.973906528517172, -0.148874338981631 }, {
// 0.679409568299024, -0.973906528517172, 0.148874338981631 }, {
// 0.679409568299024, -0.973906528517172, 0.433395394129247 }, {
// 0.679409568299024, -0.973906528517172, 0.679409568299024 }, {
// 0.679409568299024, -0.973906528517172, 0.865063366688985 }, {
// 0.679409568299024, -0.973906528517172, 0.973906528517172 }, {
// 0.679409568299024, -0.865063366688985, -0.973906528517172 }, {
// 0.679409568299024, -0.865063366688985, -0.865063366688985 }, {
// 0.679409568299024, -0.865063366688985, -0.679409568299024 }, {
// 0.679409568299024, -0.865063366688985, -0.433395394129247 }, {
// 0.679409568299024, -0.865063366688985, -0.148874338981631 }, {
// 0.679409568299024, -0.865063366688985, 0.148874338981631 }, {
// 0.679409568299024, -0.865063366688985, 0.433395394129247 }, {
// 0.679409568299024, -0.865063366688985, 0.679409568299024 }, {
// 0.679409568299024, -0.865063366688985, 0.865063366688985 }, {
// 0.679409568299024, -0.865063366688985, 0.973906528517172 }, {
// 0.679409568299024, -0.679409568299024, -0.973906528517172 }, {
// 0.679409568299024, -0.679409568299024, -0.865063366688985 }, {
// 0.679409568299024, -0.679409568299024, -0.679409568299024 }, {
// 0.679409568299024, -0.679409568299024, -0.433395394129247 }, {
// 0.679409568299024, -0.679409568299024, -0.148874338981631 }, {
// 0.679409568299024, -0.679409568299024, 0.148874338981631 }, {
// 0.679409568299024, -0.679409568299024, 0.433395394129247 }, {
// 0.679409568299024, -0.679409568299024, 0.679409568299024 }, {
// 0.679409568299024, -0.679409568299024, 0.865063366688985 }, {
// 0.679409568299024, -0.679409568299024, 0.973906528517172 }, {
// 0.679409568299024, -0.433395394129247, -0.973906528517172 }, {
// 0.679409568299024, -0.433395394129247, -0.865063366688985 }, {
// 0.679409568299024, -0.433395394129247, -0.679409568299024 }, {
// 0.679409568299024, -0.433395394129247, -0.433395394129247 }, {
// 0.679409568299024, -0.433395394129247, -0.148874338981631 }, {
// 0.679409568299024, -0.433395394129247, 0.148874338981631 }, {
// 0.679409568299024, -0.433395394129247, 0.433395394129247 }, {
// 0.679409568299024, -0.433395394129247, 0.679409568299024 }, {
// 0.679409568299024, -0.433395394129247, 0.865063366688985 }, {
// 0.679409568299024, -0.433395394129247, 0.973906528517172 }, {
// 0.679409568299024, -0.148874338981631, -0.973906528517172 }, {
// 0.679409568299024, -0.148874338981631, -0.865063366688985 }, {
// 0.679409568299024, -0.148874338981631, -0.679409568299024 }, {
// 0.679409568299024, -0.148874338981631, -0.433395394129247 }, {
// 0.679409568299024, -0.148874338981631, -0.148874338981631 }, {
// 0.679409568299024, -0.148874338981631, 0.148874338981631 }, {
// 0.679409568299024, -0.148874338981631, 0.433395394129247 }, {
// 0.679409568299024, -0.148874338981631, 0.679409568299024 }, {
// 0.679409568299024, -0.148874338981631, 0.865063366688985 }, {
// 0.679409568299024, -0.148874338981631, 0.973906528517172 }, {
// 0.679409568299024, 0.148874338981631, -0.973906528517172 }, {
// 0.679409568299024, 0.148874338981631, -0.865063366688985 }, {
// 0.679409568299024, 0.148874338981631, -0.679409568299024 }, {
// 0.679409568299024, 0.148874338981631, -0.433395394129247 }, {
// 0.679409568299024, 0.148874338981631, -0.148874338981631 }, {
// 0.679409568299024, 0.148874338981631, 0.148874338981631 }, {
// 0.679409568299024, 0.148874338981631, 0.433395394129247 }, {
// 0.679409568299024, 0.148874338981631, 0.679409568299024 }, {
// 0.679409568299024, 0.148874338981631, 0.865063366688985 }, {
// 0.679409568299024, 0.148874338981631, 0.973906528517172 }, {
// 0.679409568299024, 0.433395394129247, -0.973906528517172 }, {
// 0.679409568299024, 0.433395394129247, -0.865063366688985 }, {
// 0.679409568299024, 0.433395394129247, -0.679409568299024 }, {
// 0.679409568299024, 0.433395394129247, -0.433395394129247 }, {
// 0.679409568299024, 0.433395394129247, -0.148874338981631 }, {
// 0.679409568299024, 0.433395394129247, 0.148874338981631 }, {
// 0.679409568299024, 0.433395394129247, 0.433395394129247 }, {
// 0.679409568299024, 0.433395394129247, 0.679409568299024 }, {
// 0.679409568299024, 0.433395394129247, 0.865063366688985 }, {
// 0.679409568299024, 0.433395394129247, 0.973906528517172 }, {
// 0.679409568299024, 0.679409568299024, -0.973906528517172 }, {
// 0.679409568299024, 0.679409568299024, -0.865063366688985 }, {
// 0.679409568299024, 0.679409568299024, -0.679409568299024 }, {
// 0.679409568299024, 0.679409568299024, -0.433395394129247 }, {
// 0.679409568299024, 0.679409568299024, -0.148874338981631 }, {
// 0.679409568299024, 0.679409568299024, 0.148874338981631 }, {
// 0.679409568299024, 0.679409568299024, 0.433395394129247 }, {
// 0.679409568299024, 0.679409568299024, 0.679409568299024 }, {
// 0.679409568299024, 0.679409568299024, 0.865063366688985 }, {
// 0.679409568299024, 0.679409568299024, 0.973906528517172 }, {
// 0.679409568299024, 0.865063366688985, -0.973906528517172 }, {
// 0.679409568299024, 0.865063366688985, -0.865063366688985 }, {
// 0.679409568299024, 0.865063366688985, -0.679409568299024 }, {
// 0.679409568299024, 0.865063366688985, -0.433395394129247 }, {
// 0.679409568299024, 0.865063366688985, -0.148874338981631 }, {
// 0.679409568299024, 0.865063366688985, 0.148874338981631 }, {
// 0.679409568299024, 0.865063366688985, 0.433395394129247 }, {
// 0.679409568299024, 0.865063366688985, 0.679409568299024 }, {
// 0.679409568299024, 0.865063366688985, 0.865063366688985 }, {
// 0.679409568299024, 0.865063366688985, 0.973906528517172 }, {
// 0.679409568299024, 0.973906528517172, -0.973906528517172 }, {
// 0.679409568299024, 0.973906528517172, -0.865063366688985 }, {
// 0.679409568299024, 0.973906528517172, -0.679409568299024 }, {
// 0.679409568299024, 0.973906528517172, -0.433395394129247 }, {
// 0.679409568299024, 0.973906528517172, -0.148874338981631 }, {
// 0.679409568299024, 0.973906528517172, 0.148874338981631 }, {
// 0.679409568299024, 0.973906528517172, 0.433395394129247 }, {
// 0.679409568299024, 0.973906528517172, 0.679409568299024 }, {
// 0.679409568299024, 0.973906528517172, 0.865063366688985 }, {
// 0.679409568299024, 0.973906528517172, 0.973906528517172 }, {
// 0.865063366688985, -0.973906528517172, -0.973906528517172 }, {
// 0.865063366688985, -0.973906528517172, -0.865063366688985 }, {
// 0.865063366688985, -0.973906528517172, -0.679409568299024 }, {
// 0.865063366688985, -0.973906528517172, -0.433395394129247 }, {
// 0.865063366688985, -0.973906528517172, -0.148874338981631 }, {
// 0.865063366688985, -0.973906528517172, 0.148874338981631 }, {
// 0.865063366688985, -0.973906528517172, 0.433395394129247 }, {
// 0.865063366688985, -0.973906528517172, 0.679409568299024 }, {
// 0.865063366688985, -0.973906528517172, 0.865063366688985 }, {
// 0.865063366688985, -0.973906528517172, 0.973906528517172 }, {
// 0.865063366688985, -0.865063366688985, -0.973906528517172 }, {
// 0.865063366688985, -0.865063366688985, -0.865063366688985 }, {
// 0.865063366688985, -0.865063366688985, -0.679409568299024 }, {
// 0.865063366688985, -0.865063366688985, -0.433395394129247 }, {
// 0.865063366688985, -0.865063366688985, -0.148874338981631 }, {
// 0.865063366688985, -0.865063366688985, 0.148874338981631 }, {
// 0.865063366688985, -0.865063366688985, 0.433395394129247 }, {
// 0.865063366688985, -0.865063366688985, 0.679409568299024 }, {
// 0.865063366688985, -0.865063366688985, 0.865063366688985 }, {
// 0.865063366688985, -0.865063366688985, 0.973906528517172 }, {
// 0.865063366688985, -0.679409568299024, -0.973906528517172 }, {
// 0.865063366688985, -0.679409568299024, -0.865063366688985 }, {
// 0.865063366688985, -0.679409568299024, -0.679409568299024 }, {
// 0.865063366688985, -0.679409568299024, -0.433395394129247 }, {
// 0.865063366688985, -0.679409568299024, -0.148874338981631 }, {
// 0.865063366688985, -0.679409568299024, 0.148874338981631 }, {
// 0.865063366688985, -0.679409568299024, 0.433395394129247 }, {
// 0.865063366688985, -0.679409568299024, 0.679409568299024 }, {
// 0.865063366688985, -0.679409568299024, 0.865063366688985 }, {
// 0.865063366688985, -0.679409568299024, 0.973906528517172 }, {
// 0.865063366688985, -0.433395394129247, -0.973906528517172 }, {
// 0.865063366688985, -0.433395394129247, -0.865063366688985 }, {
// 0.865063366688985, -0.433395394129247, -0.679409568299024 }, {
// 0.865063366688985, -0.433395394129247, -0.433395394129247 }, {
// 0.865063366688985, -0.433395394129247, -0.148874338981631 }, {
// 0.865063366688985, -0.433395394129247, 0.148874338981631 }, {
// 0.865063366688985, -0.433395394129247, 0.433395394129247 }, {
// 0.865063366688985, -0.433395394129247, 0.679409568299024 }, {
// 0.865063366688985, -0.433395394129247, 0.865063366688985 }, {
// 0.865063366688985, -0.433395394129247, 0.973906528517172 }, {
// 0.865063366688985, -0.148874338981631, -0.973906528517172 }, {
// 0.865063366688985, -0.148874338981631, -0.865063366688985 }, {
// 0.865063366688985, -0.148874338981631, -0.679409568299024 }, {
// 0.865063366688985, -0.148874338981631, -0.433395394129247 }, {
// 0.865063366688985, -0.148874338981631, -0.148874338981631 }, {
// 0.865063366688985, -0.148874338981631, 0.148874338981631 }, {
// 0.865063366688985, -0.148874338981631, 0.433395394129247 }, {
// 0.865063366688985, -0.148874338981631, 0.679409568299024 }, {
// 0.865063366688985, -0.148874338981631, 0.865063366688985 }, {
// 0.865063366688985, -0.148874338981631, 0.973906528517172 }, {
// 0.865063366688985, 0.148874338981631, -0.973906528517172 }, {
// 0.865063366688985, 0.148874338981631, -0.865063366688985 }, {
// 0.865063366688985, 0.148874338981631, -0.679409568299024 }, {
// 0.865063366688985, 0.148874338981631, -0.433395394129247 }, {
// 0.865063366688985, 0.148874338981631, -0.148874338981631 }, {
// 0.865063366688985, 0.148874338981631, 0.148874338981631 }, {
// 0.865063366688985, 0.148874338981631, 0.433395394129247 }, {
// 0.865063366688985, 0.148874338981631, 0.679409568299024 }, {
// 0.865063366688985, 0.148874338981631, 0.865063366688985 }, {
// 0.865063366688985, 0.148874338981631, 0.973906528517172 }, {
// 0.865063366688985, 0.433395394129247, -0.973906528517172 }, {
// 0.865063366688985, 0.433395394129247, -0.865063366688985 }, {
// 0.865063366688985, 0.433395394129247, -0.679409568299024 }, {
// 0.865063366688985, 0.433395394129247, -0.433395394129247 }, {
// 0.865063366688985, 0.433395394129247, -0.148874338981631 }, {
// 0.865063366688985, 0.433395394129247, 0.148874338981631 }, {
// 0.865063366688985, 0.433395394129247, 0.433395394129247 }, {
// 0.865063366688985, 0.433395394129247, 0.679409568299024 }, {
// 0.865063366688985, 0.433395394129247, 0.865063366688985 }, {
// 0.865063366688985, 0.433395394129247, 0.973906528517172 }, {
// 0.865063366688985, 0.679409568299024, -0.973906528517172 }, {
// 0.865063366688985, 0.679409568299024, -0.865063366688985 }, {
// 0.865063366688985, 0.679409568299024, -0.679409568299024 }, {
// 0.865063366688985, 0.679409568299024, -0.433395394129247 }, {
// 0.865063366688985, 0.679409568299024, -0.148874338981631 }, {
// 0.865063366688985, 0.679409568299024, 0.148874338981631 }, {
// 0.865063366688985, 0.679409568299024, 0.433395394129247 }, {
// 0.865063366688985, 0.679409568299024, 0.679409568299024 }, {
// 0.865063366688985, 0.679409568299024, 0.865063366688985 }, {
// 0.865063366688985, 0.679409568299024, 0.973906528517172 }, {
// 0.865063366688985, 0.865063366688985, -0.973906528517172 }, {
// 0.865063366688985, 0.865063366688985, -0.865063366688985 }, {
// 0.865063366688985, 0.865063366688985, -0.679409568299024 }, {
// 0.865063366688985, 0.865063366688985, -0.433395394129247 }, {
// 0.865063366688985, 0.865063366688985, -0.148874338981631 }, {
// 0.865063366688985, 0.865063366688985, 0.148874338981631 }, {
// 0.865063366688985, 0.865063366688985, 0.433395394129247 }, {
// 0.865063366688985, 0.865063366688985, 0.679409568299024 }, {
// 0.865063366688985, 0.865063366688985, 0.865063366688985 }, {
// 0.865063366688985, 0.865063366688985, 0.973906528517172 }, {
// 0.865063366688985, 0.973906528517172, -0.973906528517172 }, {
// 0.865063366688985, 0.973906528517172, -0.865063366688985 }, {
// 0.865063366688985, 0.973906528517172, -0.679409568299024 }, {
// 0.865063366688985, 0.973906528517172, -0.433395394129247 }, {
// 0.865063366688985, 0.973906528517172, -0.148874338981631 }, {
// 0.865063366688985, 0.973906528517172, 0.148874338981631 }, {
// 0.865063366688985, 0.973906528517172, 0.433395394129247 }, {
// 0.865063366688985, 0.973906528517172, 0.679409568299024 }, {
// 0.865063366688985, 0.973906528517172, 0.865063366688985 }, {
// 0.865063366688985, 0.973906528517172, 0.973906528517172 }, {
// 0.973906528517172, -0.973906528517172, -0.973906528517172 }, {
// 0.973906528517172, -0.973906528517172, -0.865063366688985 }, {
// 0.973906528517172, -0.973906528517172, -0.679409568299024 }, {
// 0.973906528517172, -0.973906528517172, -0.433395394129247 }, {
// 0.973906528517172, -0.973906528517172, -0.148874338981631 }, {
// 0.973906528517172, -0.973906528517172, 0.148874338981631 }, {
// 0.973906528517172, -0.973906528517172, 0.433395394129247 }, {
// 0.973906528517172, -0.973906528517172, 0.679409568299024 }, {
// 0.973906528517172, -0.973906528517172, 0.865063366688985 }, {
// 0.973906528517172, -0.973906528517172, 0.973906528517172 }, {
// 0.973906528517172, -0.865063366688985, -0.973906528517172 }, {
// 0.973906528517172, -0.865063366688985, -0.865063366688985 }, {
// 0.973906528517172, -0.865063366688985, -0.679409568299024 }, {
// 0.973906528517172, -0.865063366688985, -0.433395394129247 }, {
// 0.973906528517172, -0.865063366688985, -0.148874338981631 }, {
// 0.973906528517172, -0.865063366688985, 0.148874338981631 }, {
// 0.973906528517172, -0.865063366688985, 0.433395394129247 }, {
// 0.973906528517172, -0.865063366688985, 0.679409568299024 }, {
// 0.973906528517172, -0.865063366688985, 0.865063366688985 }, {
// 0.973906528517172, -0.865063366688985, 0.973906528517172 }, {
// 0.973906528517172, -0.679409568299024, -0.973906528517172 }, {
// 0.973906528517172, -0.679409568299024, -0.865063366688985 }, {
// 0.973906528517172, -0.679409568299024, -0.679409568299024 }, {
// 0.973906528517172, -0.679409568299024, -0.433395394129247 }, {
// 0.973906528517172, -0.679409568299024, -0.148874338981631 }, {
// 0.973906528517172, -0.679409568299024, 0.148874338981631 }, {
// 0.973906528517172, -0.679409568299024, 0.433395394129247 }, {
// 0.973906528517172, -0.679409568299024, 0.679409568299024 }, {
// 0.973906528517172, -0.679409568299024, 0.865063366688985 }, {
// 0.973906528517172, -0.679409568299024, 0.973906528517172 }, {
// 0.973906528517172, -0.433395394129247, -0.973906528517172 }, {
// 0.973906528517172, -0.433395394129247, -0.865063366688985 }, {
// 0.973906528517172, -0.433395394129247, -0.679409568299024 }, {
// 0.973906528517172, -0.433395394129247, -0.433395394129247 }, {
// 0.973906528517172, -0.433395394129247, -0.148874338981631 }, {
// 0.973906528517172, -0.433395394129247, 0.148874338981631 }, {
// 0.973906528517172, -0.433395394129247, 0.433395394129247 }, {
// 0.973906528517172, -0.433395394129247, 0.679409568299024 }, {
// 0.973906528517172, -0.433395394129247, 0.865063366688985 }, {
// 0.973906528517172, -0.433395394129247, 0.973906528517172 }, {
// 0.973906528517172, -0.148874338981631, -0.973906528517172 }, {
// 0.973906528517172, -0.148874338981631, -0.865063366688985 }, {
// 0.973906528517172, -0.148874338981631, -0.679409568299024 }, {
// 0.973906528517172, -0.148874338981631, -0.433395394129247 }, {
// 0.973906528517172, -0.148874338981631, -0.148874338981631 }, {
// 0.973906528517172, -0.148874338981631, 0.148874338981631 }, {
// 0.973906528517172, -0.148874338981631, 0.433395394129247 }, {
// 0.973906528517172, -0.148874338981631, 0.679409568299024 }, {
// 0.973906528517172, -0.148874338981631, 0.865063366688985 }, {
// 0.973906528517172, -0.148874338981631, 0.973906528517172 }, {
// 0.973906528517172, 0.148874338981631, -0.973906528517172 }, {
// 0.973906528517172, 0.148874338981631, -0.865063366688985 }, {
// 0.973906528517172, 0.148874338981631, -0.679409568299024 }, {
// 0.973906528517172, 0.148874338981631, -0.433395394129247 }, {
// 0.973906528517172, 0.148874338981631, -0.148874338981631 }, {
// 0.973906528517172, 0.148874338981631, 0.148874338981631 }, {
// 0.973906528517172, 0.148874338981631, 0.433395394129247 }, {
// 0.973906528517172, 0.148874338981631, 0.679409568299024 }, {
// 0.973906528517172, 0.148874338981631, 0.865063366688985 }, {
// 0.973906528517172, 0.148874338981631, 0.973906528517172 }, {
// 0.973906528517172, 0.433395394129247, -0.973906528517172 }, {
// 0.973906528517172, 0.433395394129247, -0.865063366688985 }, {
// 0.973906528517172, 0.433395394129247, -0.679409568299024 }, {
// 0.973906528517172, 0.433395394129247, -0.433395394129247 }, {
// 0.973906528517172, 0.433395394129247, -0.148874338981631 }, {
// 0.973906528517172, 0.433395394129247, 0.148874338981631 }, {
// 0.973906528517172, 0.433395394129247, 0.433395394129247 }, {
// 0.973906528517172, 0.433395394129247, 0.679409568299024 }, {
// 0.973906528517172, 0.433395394129247, 0.865063366688985 }, {
// 0.973906528517172, 0.433395394129247, 0.973906528517172 }, {
// 0.973906528517172, 0.679409568299024, -0.973906528517172 }, {
// 0.973906528517172, 0.679409568299024, -0.865063366688985 }, {
// 0.973906528517172, 0.679409568299024, -0.679409568299024 }, {
// 0.973906528517172, 0.679409568299024, -0.433395394129247 }, {
// 0.973906528517172, 0.679409568299024, -0.148874338981631 }, {
// 0.973906528517172, 0.679409568299024, 0.148874338981631 }, {
// 0.973906528517172, 0.679409568299024, 0.433395394129247 }, {
// 0.973906528517172, 0.679409568299024, 0.679409568299024 }, {
// 0.973906528517172, 0.679409568299024, 0.865063366688985 }, {
// 0.973906528517172, 0.679409568299024, 0.973906528517172 }, {
// 0.973906528517172, 0.865063366688985, -0.973906528517172 }, {
// 0.973906528517172, 0.865063366688985, -0.865063366688985 }, {
// 0.973906528517172, 0.865063366688985, -0.679409568299024 }, {
// 0.973906528517172, 0.865063366688985, -0.433395394129247 }, {
// 0.973906528517172, 0.865063366688985, -0.148874338981631 }, {
// 0.973906528517172, 0.865063366688985, 0.148874338981631 }, {
// 0.973906528517172, 0.865063366688985, 0.433395394129247 }, {
// 0.973906528517172, 0.865063366688985, 0.679409568299024 }, {
// 0.973906528517172, 0.865063366688985, 0.865063366688985 }, {
// 0.973906528517172, 0.865063366688985, 0.973906528517172 }, {
// 0.973906528517172, 0.973906528517172, -0.973906528517172 }, {
// 0.973906528517172, 0.973906528517172, -0.865063366688985 }, {
// 0.973906528517172, 0.973906528517172, -0.679409568299024 }, {
// 0.973906528517172, 0.973906528517172, -0.433395394129247 }, {
// 0.973906528517172, 0.973906528517172, -0.148874338981631 }, {
// 0.973906528517172, 0.973906528517172, 0.148874338981631 }, {
// 0.973906528517172, 0.973906528517172, 0.433395394129247 }, {
// 0.973906528517172, 0.973906528517172, 0.679409568299024 }, {
// 0.973906528517172, 0.973906528517172, 0.865063366688985 }, {
// 0.973906528517172, 0.973906528517172, 0.973906528517172 } } , {
// 0.0002963586692327501, 0.0006643214323718632, 0.0009738538125414616,
// 0.001196908918378898, 0.001313625319402651, 0.001313625319402651,
// 0.001196908918378898, 0.0009738538125414616, 0.0006643214323718632,
// 0.0002963586692327501, 0.0006643214323718632, 0.001489151529297777,
// 0.002183003322775247, 0.002683006537769466, 0.002944639534401979,
// 0.002944639534401979, 0.002683006537769466, 0.002183003322775247,
// 0.001489151529297777, 0.0006643214323718632, 0.0009738538125414615,
// 0.002183003322775247, 0.003200146804062971, 0.003933120351923071,
// 0.004316658017338323, 0.004316658017338323, 0.003933120351923071,
// 0.003200146804062971, 0.002183003322775247, 0.0009738538125414615,
// 0.001196908918378898, 0.002683006537769466, 0.003933120351923071,
// 0.004833976892269804, 0.005305361453646353, 0.005305361453646353,
// 0.004833976892269804, 0.003933120351923071, 0.002683006537769466,
// 0.001196908918378898, 0.001313625319402651, 0.002944639534401979,
// 0.004316658017338323, 0.005305361453646353, 0.005822713012726075,
// 0.005822713012726075, 0.005305361453646353, 0.004316658017338323,
// 0.002944639534401979, 0.001313625319402651, 0.001313625319402651,
// 0.002944639534401979, 0.004316658017338323, 0.005305361453646353,
// 0.005822713012726075, 0.005822713012726075, 0.005305361453646353,
// 0.004316658017338323, 0.002944639534401979, 0.001313625319402651,
// 0.001196908918378898, 0.002683006537769466, 0.003933120351923071,
// 0.004833976892269804, 0.005305361453646353, 0.005305361453646353,
// 0.004833976892269804, 0.003933120351923071, 0.002683006537769466,
// 0.001196908918378898, 0.0009738538125414615, 0.002183003322775247,
// 0.003200146804062971, 0.003933120351923071, 0.004316658017338323,
// 0.004316658017338323, 0.003933120351923071, 0.003200146804062971,
// 0.002183003322775247, 0.0009738538125414615, 0.0006643214323718632,
// 0.001489151529297777, 0.002183003322775247, 0.002683006537769466,
// 0.002944639534401979, 0.002944639534401979, 0.002683006537769466,
// 0.002183003322775247, 0.001489151529297777, 0.0006643214323718632,
// 0.0002963586692327501, 0.0006643214323718632, 0.0009738538125414616,
// 0.001196908918378898, 0.001313625319402651, 0.001313625319402651,
// 0.001196908918378898, 0.0009738538125414616, 0.0006643214323718632,
// 0.0002963586692327501, 0.0006643214323718632, 0.001489151529297777,
// 0.002183003322775247, 0.002683006537769466, 0.002944639534401979,
// 0.002944639534401979, 0.002683006537769466, 0.002183003322775247,
// 0.001489151529297777, 0.0006643214323718632, 0.001489151529297777,
// 0.003338101360500127, 0.004893448529827343, 0.006014262214257229,
// 0.006600742129046395, 0.006600742129046395, 0.006014262214257229,
// 0.004893448529827343, 0.003338101360500127, 0.001489151529297777,
// 0.002183003322775246, 0.004893448529827342, 0.00717349053489548,
// 0.008816533535681429, 0.00967627653532729, 0.00967627653532729,
// 0.008816533535681429, 0.00717349053489548, 0.004893448529827342,
// 0.002183003322775246, 0.002683006537769466, 0.006014262214257229,
// 0.008816533535681431, 0.01083590522740236, 0.01189256696711922,
// 0.01189256696711922, 0.01083590522740236, 0.008816533535681431,
// 0.006014262214257229, 0.002683006537769466, 0.002944639534401979,
// 0.006600742129046395, 0.009676276535327292, 0.01189256696711922,
// 0.01305226892440442, 0.01305226892440442, 0.01189256696711922,
// 0.009676276535327292, 0.006600742129046395, 0.002944639534401979,
// 0.002944639534401979, 0.006600742129046395, 0.009676276535327292,
// 0.01189256696711922, 0.01305226892440442, 0.01305226892440442,
// 0.01189256696711922, 0.009676276535327292, 0.006600742129046395,
// 0.002944639534401979, 0.002683006537769466, 0.006014262214257229,
// 0.008816533535681431, 0.01083590522740236, 0.01189256696711922,
// 0.01189256696711922, 0.01083590522740236, 0.008816533535681431,
// 0.006014262214257229, 0.002683006537769466, 0.002183003322775246,
// 0.004893448529827342, 0.00717349053489548, 0.008816533535681429,
// 0.00967627653532729, 0.00967627653532729, 0.008816533535681429,
// 0.00717349053489548, 0.004893448529827342, 0.002183003322775246,
// 0.001489151529297777, 0.003338101360500127, 0.004893448529827343,
// 0.006014262214257229, 0.006600742129046395, 0.006600742129046395,
// 0.006014262214257229, 0.004893448529827343, 0.003338101360500127,
// 0.001489151529297777, 0.0006643214323718632, 0.001489151529297777,
// 0.002183003322775247, 0.002683006537769466, 0.002944639534401979,
// 0.002944639534401979, 0.002683006537769466, 0.002183003322775247,
// 0.001489151529297777, 0.0006643214323718632, 0.0009738538125414615,
// 0.002183003322775247, 0.003200146804062971, 0.003933120351923071,
// 0.004316658017338323, 0.004316658017338323, 0.003933120351923071,
// 0.003200146804062971, 0.002183003322775247, 0.0009738538125414615,
// 0.002183003322775246, 0.004893448529827342, 0.00717349053489548,
// 0.008816533535681429, 0.00967627653532729, 0.00967627653532729,
// 0.008816533535681429, 0.00717349053489548, 0.004893448529827342,
// 0.002183003322775246, 0.003200146804062971, 0.00717349053489548,
// 0.01051588999875527, 0.01292448862663951, 0.01418481827613106,
// 0.01418481827613106, 0.01292448862663951, 0.01051588999875527,
// 0.00717349053489548, 0.003200146804062971, 0.003933120351923071,
// 0.008816533535681429, 0.01292448862663951, 0.01588476165877603,
// 0.01743376190721892, 0.01743376190721892, 0.01588476165877603,
// 0.01292448862663951, 0.008816533535681429, 0.003933120351923071,
// 0.004316658017338323, 0.00967627653532729, 0.01418481827613106,
// 0.01743376190721892, 0.0191338126920953, 0.0191338126920953,
// 0.01743376190721892, 0.01418481827613106, 0.00967627653532729,
// 0.004316658017338323, 0.004316658017338323, 0.00967627653532729,
// 0.01418481827613106, 0.01743376190721892, 0.0191338126920953,
// 0.0191338126920953, 0.01743376190721892, 0.01418481827613106,
// 0.00967627653532729, 0.004316658017338323, 0.003933120351923071,
// 0.008816533535681429, 0.01292448862663951, 0.01588476165877603,
// 0.01743376190721892, 0.01743376190721892, 0.01588476165877603,
// 0.01292448862663951, 0.008816533535681429, 0.003933120351923071,
// 0.003200146804062971, 0.00717349053489548, 0.01051588999875527,
// 0.01292448862663951, 0.01418481827613106, 0.01418481827613106,
// 0.01292448862663951, 0.01051588999875527, 0.00717349053489548,
// 0.003200146804062971, 0.002183003322775246, 0.004893448529827342,
// 0.00717349053489548, 0.008816533535681429, 0.00967627653532729,
// 0.00967627653532729, 0.008816533535681429, 0.00717349053489548,
// 0.004893448529827342, 0.002183003322775246, 0.0009738538125414615,
// 0.002183003322775247, 0.003200146804062971, 0.003933120351923071,
// 0.004316658017338323, 0.004316658017338323, 0.003933120351923071,
// 0.003200146804062971, 0.002183003322775247, 0.0009738538125414615,
// 0.001196908918378898, 0.002683006537769466, 0.003933120351923071,
// 0.004833976892269804, 0.005305361453646353, 0.005305361453646353,
// 0.004833976892269804, 0.003933120351923071, 0.002683006537769466,
// 0.001196908918378898, 0.002683006537769466, 0.006014262214257229,
// 0.008816533535681431, 0.01083590522740236, 0.01189256696711922,
// 0.01189256696711922, 0.01083590522740236, 0.008816533535681431,
// 0.006014262214257229, 0.002683006537769466, 0.003933120351923071,
// 0.008816533535681429, 0.01292448862663951, 0.01588476165877603,
// 0.01743376190721892, 0.01743376190721892, 0.01588476165877603,
// 0.01292448862663951, 0.008816533535681429, 0.003933120351923071,
// 0.004833976892269804, 0.01083590522740236, 0.01588476165877603,
// 0.01952306665627265, 0.02142685569324734, 0.02142685569324734,
// 0.01952306665627265, 0.01588476165877603, 0.01083590522740236,
// 0.004833976892269804, 0.005305361453646352, 0.01189256696711922,
// 0.01743376190721892, 0.02142685569324734, 0.02351629244433974,
// 0.02351629244433974, 0.02142685569324734, 0.01743376190721892,
// 0.01189256696711922, 0.005305361453646352, 0.005305361453646352,
// 0.01189256696711922, 0.01743376190721892, 0.02142685569324734,
// 0.02351629244433974, 0.02351629244433974, 0.02142685569324734,
// 0.01743376190721892, 0.01189256696711922, 0.005305361453646352,
// 0.004833976892269804, 0.01083590522740236, 0.01588476165877603,
// 0.01952306665627265, 0.02142685569324734, 0.02142685569324734,
// 0.01952306665627265, 0.01588476165877603, 0.01083590522740236,
// 0.004833976892269804, 0.003933120351923071, 0.008816533535681429,
// 0.01292448862663951, 0.01588476165877603, 0.01743376190721892,
// 0.01743376190721892, 0.01588476165877603, 0.01292448862663951,
// 0.008816533535681429, 0.003933120351923071, 0.002683006537769466,
// 0.006014262214257229, 0.008816533535681431, 0.01083590522740236,
// 0.01189256696711922, 0.01189256696711922, 0.01083590522740236,
// 0.008816533535681431, 0.006014262214257229, 0.002683006537769466,
// 0.001196908918378898, 0.002683006537769466, 0.003933120351923071,
// 0.004833976892269804, 0.005305361453646353, 0.005305361453646353,
// 0.004833976892269804, 0.003933120351923071, 0.002683006537769466,
// 0.001196908918378898, 0.001313625319402651, 0.002944639534401979,
// 0.004316658017338323, 0.005305361453646353, 0.005822713012726075,
// 0.005822713012726075, 0.005305361453646353, 0.004316658017338323,
// 0.002944639534401979, 0.001313625319402651, 0.002944639534401979,
// 0.006600742129046395, 0.009676276535327292, 0.01189256696711922,
// 0.01305226892440442, 0.01305226892440442, 0.01189256696711922,
// 0.009676276535327292, 0.006600742129046395, 0.002944639534401979,
// 0.004316658017338323, 0.00967627653532729, 0.01418481827613106,
// 0.01743376190721892, 0.0191338126920953, 0.0191338126920953,
// 0.01743376190721892, 0.01418481827613106, 0.00967627653532729,
// 0.004316658017338323, 0.005305361453646352, 0.01189256696711922,
// 0.01743376190721892, 0.02142685569324734, 0.02351629244433974,
// 0.02351629244433974, 0.02142685569324734, 0.01743376190721892,
// 0.01189256696711922, 0.005305361453646352, 0.005822713012726075,
// 0.01305226892440442, 0.01913381269209531, 0.02351629244433975,
// 0.02580948031969027, 0.02580948031969027, 0.02351629244433975,
// 0.01913381269209531, 0.01305226892440442, 0.005822713012726075,
// 0.005822713012726075, 0.01305226892440442, 0.01913381269209531,
// 0.02351629244433975, 0.02580948031969027, 0.02580948031969027,
// 0.02351629244433975, 0.01913381269209531, 0.01305226892440442,
// 0.005822713012726075, 0.005305361453646352, 0.01189256696711922,
// 0.01743376190721892, 0.02142685569324734, 0.02351629244433974,
// 0.02351629244433974, 0.02142685569324734, 0.01743376190721892,
// 0.01189256696711922, 0.005305361453646352, 0.004316658017338323,
// 0.00967627653532729, 0.01418481827613106, 0.01743376190721892,
// 0.0191338126920953, 0.0191338126920953, 0.01743376190721892,
// 0.01418481827613106, 0.00967627653532729, 0.004316658017338323,
// 0.002944639534401979, 0.006600742129046395, 0.009676276535327292,
// 0.01189256696711922, 0.01305226892440442, 0.01305226892440442,
// 0.01189256696711922, 0.009676276535327292, 0.006600742129046395,
// 0.002944639534401979, 0.001313625319402651, 0.002944639534401979,
// 0.004316658017338323, 0.005305361453646353, 0.005822713012726075,
// 0.005822713012726075, 0.005305361453646353, 0.004316658017338323,
// 0.002944639534401979, 0.001313625319402651, 0.001313625319402651,
// 0.002944639534401979, 0.004316658017338323, 0.005305361453646353,
// 0.005822713012726075, 0.005822713012726075, 0.005305361453646353,
// 0.004316658017338323, 0.002944639534401979, 0.001313625319402651,
// 0.002944639534401979, 0.006600742129046395, 0.009676276535327292,
// 0.01189256696711922, 0.01305226892440442, 0.01305226892440442,
// 0.01189256696711922, 0.009676276535327292, 0.006600742129046395,
// 0.002944639534401979, 0.004316658017338323, 0.00967627653532729,
// 0.01418481827613106, 0.01743376190721892, 0.0191338126920953,
// 0.0191338126920953, 0.01743376190721892, 0.01418481827613106,
// 0.00967627653532729, 0.004316658017338323, 0.005305361453646352,
// 0.01189256696711922, 0.01743376190721892, 0.02142685569324734,
// 0.02351629244433974, 0.02351629244433974, 0.02142685569324734,
// 0.01743376190721892, 0.01189256696711922, 0.005305361453646352,
// 0.005822713012726075, 0.01305226892440442, 0.01913381269209531,
// 0.02351629244433975, 0.02580948031969027, 0.02580948031969027,
// 0.02351629244433975, 0.01913381269209531, 0.01305226892440442,
// 0.005822713012726075, 0.005822713012726075, 0.01305226892440442,
// 0.01913381269209531, 0.02351629244433975, 0.02580948031969027,
// 0.02580948031969027, 0.02351629244433975, 0.01913381269209531,
// 0.01305226892440442, 0.005822713012726075, 0.005305361453646352,
// 0.01189256696711922, 0.01743376190721892, 0.02142685569324734,
// 0.02351629244433974, 0.02351629244433974, 0.02142685569324734,
// 0.01743376190721892, 0.01189256696711922, 0.005305361453646352,
// 0.004316658017338323, 0.00967627653532729, 0.01418481827613106,
// 0.01743376190721892, 0.0191338126920953, 0.0191338126920953,
// 0.01743376190721892, 0.01418481827613106, 0.00967627653532729,
// 0.004316658017338323, 0.002944639534401979, 0.006600742129046395,
// 0.009676276535327292, 0.01189256696711922, 0.01305226892440442,
// 0.01305226892440442, 0.01189256696711922, 0.009676276535327292,
// 0.006600742129046395, 0.002944639534401979, 0.001313625319402651,
// 0.002944639534401979, 0.004316658017338323, 0.005305361453646353,
// 0.005822713012726075, 0.005822713012726075, 0.005305361453646353,
// 0.004316658017338323, 0.002944639534401979, 0.001313625319402651,
// 0.001196908918378898, 0.002683006537769466, 0.003933120351923071,
// 0.004833976892269804, 0.005305361453646353, 0.005305361453646353,
// 0.004833976892269804, 0.003933120351923071, 0.002683006537769466,
// 0.001196908918378898, 0.002683006537769466, 0.006014262214257229,
// 0.008816533535681431, 0.01083590522740236, 0.01189256696711922,
// 0.01189256696711922, 0.01083590522740236, 0.008816533535681431,
// 0.006014262214257229, 0.002683006537769466, 0.003933120351923071,
// 0.008816533535681429, 0.01292448862663951, 0.01588476165877603,
// 0.01743376190721892, 0.01743376190721892, 0.01588476165877603,
// 0.01292448862663951, 0.008816533535681429, 0.003933120351923071,
// 0.004833976892269804, 0.01083590522740236, 0.01588476165877603,
// 0.01952306665627265, 0.02142685569324734, 0.02142685569324734,
// 0.01952306665627265, 0.01588476165877603, 0.01083590522740236,
// 0.004833976892269804, 0.005305361453646352, 0.01189256696711922,
// 0.01743376190721892, 0.02142685569324734, 0.02351629244433974,
// 0.02351629244433974, 0.02142685569324734, 0.01743376190721892,
// 0.01189256696711922, 0.005305361453646352, 0.005305361453646352,
// 0.01189256696711922, 0.01743376190721892, 0.02142685569324734,
// 0.02351629244433974, 0.02351629244433974, 0.02142685569324734,
// 0.01743376190721892, 0.01189256696711922, 0.005305361453646352,
// 0.004833976892269804, 0.01083590522740236, 0.01588476165877603,
// 0.01952306665627265, 0.02142685569324734, 0.02142685569324734,
// 0.01952306665627265, 0.01588476165877603, 0.01083590522740236,
// 0.004833976892269804, 0.003933120351923071, 0.008816533535681429,
// 0.01292448862663951, 0.01588476165877603, 0.01743376190721892,
// 0.01743376190721892, 0.01588476165877603, 0.01292448862663951,
// 0.008816533535681429, 0.003933120351923071, 0.002683006537769466,
// 0.006014262214257229, 0.008816533535681431, 0.01083590522740236,
// 0.01189256696711922, 0.01189256696711922, 0.01083590522740236,
// 0.008816533535681431, 0.006014262214257229, 0.002683006537769466,
// 0.001196908918378898, 0.002683006537769466, 0.003933120351923071,
// 0.004833976892269804, 0.005305361453646353, 0.005305361453646353,
// 0.004833976892269804, 0.003933120351923071, 0.002683006537769466,
// 0.001196908918378898, 0.0009738538125414615, 0.002183003322775247,
// 0.003200146804062971, 0.003933120351923071, 0.004316658017338323,
// 0.004316658017338323, 0.003933120351923071, 0.003200146804062971,
// 0.002183003322775247, 0.0009738538125414615, 0.002183003322775246,
// 0.004893448529827342, 0.00717349053489548, 0.008816533535681429,
// 0.00967627653532729, 0.00967627653532729, 0.008816533535681429,
// 0.00717349053489548, 0.004893448529827342, 0.002183003322775246,
// 0.003200146804062971, 0.00717349053489548, 0.01051588999875527,
// 0.01292448862663951, 0.01418481827613106, 0.01418481827613106,
// 0.01292448862663951, 0.01051588999875527, 0.00717349053489548,
// 0.003200146804062971, 0.003933120351923071, 0.008816533535681429,
// 0.01292448862663951, 0.01588476165877603, 0.01743376190721892,
// 0.01743376190721892, 0.01588476165877603, 0.01292448862663951,
// 0.008816533535681429, 0.003933120351923071, 0.004316658017338323,
// 0.00967627653532729, 0.01418481827613106, 0.01743376190721892,
// 0.0191338126920953, 0.0191338126920953, 0.01743376190721892,
// 0.01418481827613106, 0.00967627653532729, 0.004316658017338323,
// 0.004316658017338323, 0.00967627653532729, 0.01418481827613106,
// 0.01743376190721892, 0.0191338126920953, 0.0191338126920953,
// 0.01743376190721892, 0.01418481827613106, 0.00967627653532729,
// 0.004316658017338323, 0.003933120351923071, 0.008816533535681429,
// 0.01292448862663951, 0.01588476165877603, 0.01743376190721892,
// 0.01743376190721892, 0.01588476165877603, 0.01292448862663951,
// 0.008816533535681429, 0.003933120351923071, 0.003200146804062971,
// 0.00717349053489548, 0.01051588999875527, 0.01292448862663951,
// 0.01418481827613106, 0.01418481827613106, 0.01292448862663951,
// 0.01051588999875527, 0.00717349053489548, 0.003200146804062971,
// 0.002183003322775246, 0.004893448529827342, 0.00717349053489548,
// 0.008816533535681429, 0.00967627653532729, 0.00967627653532729,
// 0.008816533535681429, 0.00717349053489548, 0.004893448529827342,
// 0.002183003322775246, 0.0009738538125414615, 0.002183003322775247,
// 0.003200146804062971, 0.003933120351923071, 0.004316658017338323,
// 0.004316658017338323, 0.003933120351923071, 0.003200146804062971,
// 0.002183003322775247, 0.0009738538125414615, 0.0006643214323718632,
// 0.001489151529297777, 0.002183003322775247, 0.002683006537769466,
// 0.002944639534401979, 0.002944639534401979, 0.002683006537769466,
// 0.002183003322775247, 0.001489151529297777, 0.0006643214323718632,
// 0.001489151529297777, 0.003338101360500127, 0.004893448529827343,
// 0.006014262214257229, 0.006600742129046395, 0.006600742129046395,
// 0.006014262214257229, 0.004893448529827343, 0.003338101360500127,
// 0.001489151529297777, 0.002183003322775246, 0.004893448529827342,
// 0.00717349053489548, 0.008816533535681429, 0.00967627653532729,
// 0.00967627653532729, 0.008816533535681429, 0.00717349053489548,
// 0.004893448529827342, 0.002183003322775246, 0.002683006537769466,
// 0.006014262214257229, 0.008816533535681431, 0.01083590522740236,
// 0.01189256696711922, 0.01189256696711922, 0.01083590522740236,
// 0.008816533535681431, 0.006014262214257229, 0.002683006537769466,
// 0.002944639534401979, 0.006600742129046395, 0.009676276535327292,
// 0.01189256696711922, 0.01305226892440442, 0.01305226892440442,
// 0.01189256696711922, 0.009676276535327292, 0.006600742129046395,
// 0.002944639534401979, 0.002944639534401979, 0.006600742129046395,
// 0.009676276535327292, 0.01189256696711922, 0.01305226892440442,
// 0.01305226892440442, 0.01189256696711922, 0.009676276535327292,
// 0.006600742129046395, 0.002944639534401979, 0.002683006537769466,
// 0.006014262214257229, 0.008816533535681431, 0.01083590522740236,
// 0.01189256696711922, 0.01189256696711922, 0.01083590522740236,
// 0.008816533535681431, 0.006014262214257229, 0.002683006537769466,
// 0.002183003322775246, 0.004893448529827342, 0.00717349053489548,
// 0.008816533535681429, 0.00967627653532729, 0.00967627653532729,
// 0.008816533535681429, 0.00717349053489548, 0.004893448529827342,
// 0.002183003322775246, 0.001489151529297777, 0.003338101360500127,
// 0.004893448529827343, 0.006014262214257229, 0.006600742129046395,
// 0.006600742129046395, 0.006014262214257229, 0.004893448529827343,
// 0.003338101360500127, 0.001489151529297777, 0.0006643214323718632,
// 0.001489151529297777, 0.002183003322775247, 0.002683006537769466,
// 0.002944639534401979, 0.002944639534401979, 0.002683006537769466,
// 0.002183003322775247, 0.001489151529297777, 0.0006643214323718632,
// 0.0002963586692327501, 0.0006643214323718632, 0.0009738538125414616,
// 0.001196908918378898, 0.001313625319402651, 0.001313625319402651,
// 0.001196908918378898, 0.0009738538125414616, 0.0006643214323718632,
// 0.0002963586692327501, 0.0006643214323718632, 0.001489151529297777,
// 0.002183003322775247, 0.002683006537769466, 0.002944639534401979,
// 0.002944639534401979, 0.002683006537769466, 0.002183003322775247,
// 0.001489151529297777, 0.0006643214323718632, 0.0009738538125414615,
// 0.002183003322775247, 0.003200146804062971, 0.003933120351923071,
// 0.004316658017338323, 0.004316658017338323, 0.003933120351923071,
// 0.003200146804062971, 0.002183003322775247, 0.0009738538125414615,
// 0.001196908918378898, 0.002683006537769466, 0.003933120351923071,
// 0.004833976892269804, 0.005305361453646353, 0.005305361453646353,
// 0.004833976892269804, 0.003933120351923071, 0.002683006537769466,
// 0.001196908918378898, 0.001313625319402651, 0.002944639534401979,
// 0.004316658017338323, 0.005305361453646353, 0.005822713012726075,
// 0.005822713012726075, 0.005305361453646353, 0.004316658017338323,
// 0.002944639534401979, 0.001313625319402651, 0.001313625319402651,
// 0.002944639534401979, 0.004316658017338323, 0.005305361453646353,
// 0.005822713012726075, 0.005822713012726075, 0.005305361453646353,
// 0.004316658017338323, 0.002944639534401979, 0.001313625319402651,
// 0.001196908918378898, 0.002683006537769466, 0.003933120351923071,
// 0.004833976892269804, 0.005305361453646353, 0.005305361453646353,
// 0.004833976892269804, 0.003933120351923071, 0.002683006537769466,
// 0.001196908918378898, 0.0009738538125414615, 0.002183003322775247,
// 0.003200146804062971, 0.003933120351923071, 0.004316658017338323,
// 0.004316658017338323, 0.003933120351923071, 0.003200146804062971,
// 0.002183003322775247, 0.0009738538125414615, 0.0006643214323718632,
// 0.001489151529297777, 0.002183003322775247, 0.002683006537769466,
// 0.002944639534401979, 0.002944639534401979, 0.002683006537769466,
// 0.002183003322775247, 0.001489151529297777, 0.0006643214323718632,
// 0.0002963586692327501, 0.0006643214323718632, 0.0009738538125414616,
// 0.001196908918378898, 0.001313625319402651, 0.001313625319402651,
// 0.001196908918378898, 0.0009738538125414616, 0.0006643214323718632,
// 0.0002963586692327501 } }; 					case 20:
// case 21:	return { { { -0.978228658146057, -0.978228658146057,
//-0.978228658146057 }, { -0.978228658146057, -0.978228658146057,
//-0.887062599768095 }, { -0.978228658146057, -0.978228658146057,
//-0.730152005574049 }, { -0.978228658146057, -0.978228658146057,
//-0.519096129206812 }, { -0.978228658146057, -0.978228658146057,
//-0.269543155952345 }, { -0.978228658146057, -0.978228658146057, 0 }, {
//-0.978228658146057, -0.978228658146057, 0.269543155952345 }, {
// -0.978228658146057, -0.978228658146057, 0.519096129206812 }, {
// -0.978228658146057, -0.978228658146057, 0.730152005574049 }, {
// -0.978228658146057, -0.978228658146057, 0.887062599768095 }, {
// -0.978228658146057, -0.978228658146057, 0.978228658146057 }, {
// -0.978228658146057, -0.887062599768095, -0.978228658146057 }, {
//-0.978228658146057, -0.887062599768095, -0.887062599768095 }, {
//-0.978228658146057, -0.887062599768095, -0.730152005574049 }, {
//-0.978228658146057, -0.887062599768095, -0.519096129206812 }, {
//-0.978228658146057, -0.887062599768095, -0.269543155952345 }, {
//-0.978228658146057, -0.887062599768095, 0 }, { -0.978228658146057,
//-0.887062599768095, 0.269543155952345 }, { -0.978228658146057,
//-0.887062599768095, 0.519096129206812 }, { -0.978228658146057,
//-0.887062599768095, 0.730152005574049 }, { -0.978228658146057,
//-0.887062599768095, 0.887062599768095 }, { -0.978228658146057,
//-0.887062599768095, 0.978228658146057 }, { -0.978228658146057,
//-0.730152005574049, -0.978228658146057 }, { -0.978228658146057,
//-0.730152005574049, -0.887062599768095 }, { -0.978228658146057,
//-0.730152005574049, -0.730152005574049 }, { -0.978228658146057,
//-0.730152005574049, -0.519096129206812 }, { -0.978228658146057,
//-0.730152005574049, -0.269543155952345 }, { -0.978228658146057,
//-0.730152005574049, 0 }, { -0.978228658146057, -0.730152005574049,
// 0.269543155952345 }, { -0.978228658146057, -0.730152005574049,
// 0.519096129206812 }, { -0.978228658146057, -0.730152005574049,
// 0.730152005574049 }, { -0.978228658146057, -0.730152005574049,
// 0.887062599768095 }, { -0.978228658146057, -0.730152005574049,
// 0.978228658146057 }, { -0.978228658146057, -0.519096129206812,
//-0.978228658146057 }, { -0.978228658146057, -0.519096129206812,
//-0.887062599768095 }, { -0.978228658146057, -0.519096129206812,
//-0.730152005574049 }, { -0.978228658146057, -0.519096129206812,
//-0.519096129206812 }, { -0.978228658146057, -0.519096129206812,
//-0.269543155952345 }, { -0.978228658146057, -0.519096129206812, 0 }, {
//-0.978228658146057, -0.519096129206812, 0.269543155952345 }, {
// -0.978228658146057, -0.519096129206812, 0.519096129206812 }, {
// -0.978228658146057, -0.519096129206812, 0.730152005574049 }, {
// -0.978228658146057, -0.519096129206812, 0.887062599768095 }, {
// -0.978228658146057, -0.519096129206812, 0.978228658146057 }, {
// -0.978228658146057, -0.269543155952345, -0.978228658146057 }, {
//-0.978228658146057, -0.269543155952345, -0.887062599768095 }, {
//-0.978228658146057, -0.269543155952345, -0.730152005574049 }, {
//-0.978228658146057, -0.269543155952345, -0.519096129206812 }, {
//-0.978228658146057, -0.269543155952345, -0.269543155952345 }, {
//-0.978228658146057, -0.269543155952345, 0 }, { -0.978228658146057,
//-0.269543155952345, 0.269543155952345 }, { -0.978228658146057,
//-0.269543155952345, 0.519096129206812 }, { -0.978228658146057,
//-0.269543155952345, 0.730152005574049 }, { -0.978228658146057,
//-0.269543155952345, 0.887062599768095 }, { -0.978228658146057,
//-0.269543155952345, 0.978228658146057 }, { -0.978228658146057, 0,
//-0.978228658146057 }, { -0.978228658146057, 0, -0.887062599768095 }, {
//-0.978228658146057, 0, -0.730152005574049 }, { -0.978228658146057, 0,
//-0.519096129206812 }, { -0.978228658146057, 0, -0.269543155952345 }, {
//-0.978228658146057, 0, 0 }, { -0.978228658146057, 0, 0.269543155952345 }, {
// -0.978228658146057, 0, 0.519096129206812 }, { -0.978228658146057, 0,
// 0.730152005574049 }, { -0.978228658146057, 0, 0.887062599768095 }, {
// -0.978228658146057, 0, 0.978228658146057 }, { -0.978228658146057,
// 0.269543155952345, -0.978228658146057 }, { -0.978228658146057,
// 0.269543155952345, -0.887062599768095 }, { -0.978228658146057,
// 0.269543155952345, -0.730152005574049 }, { -0.978228658146057,
// 0.269543155952345, -0.519096129206812 }, { -0.978228658146057,
// 0.269543155952345, -0.269543155952345 }, { -0.978228658146057,
// 0.269543155952345, 0 }, { -0.978228658146057, 0.269543155952345,
// 0.269543155952345 }, { -0.978228658146057, 0.269543155952345,
// 0.519096129206812 }, { -0.978228658146057, 0.269543155952345,
// 0.730152005574049 }, { -0.978228658146057, 0.269543155952345,
// 0.887062599768095 }, { -0.978228658146057, 0.269543155952345,
// 0.978228658146057 }, { -0.978228658146057, 0.519096129206812,
// -0.978228658146057 }, { -0.978228658146057, 0.519096129206812,
// -0.887062599768095 }, { -0.978228658146057, 0.519096129206812,
// -0.730152005574049 }, { -0.978228658146057, 0.519096129206812,
// -0.519096129206812 }, { -0.978228658146057, 0.519096129206812,
// -0.269543155952345 }, { -0.978228658146057, 0.519096129206812, 0 }, {
// -0.978228658146057, 0.519096129206812, 0.269543155952345 }, {
// -0.978228658146057, 0.519096129206812, 0.519096129206812 }, {
// -0.978228658146057, 0.519096129206812, 0.730152005574049 }, {
// -0.978228658146057, 0.519096129206812, 0.887062599768095 }, {
// -0.978228658146057, 0.519096129206812, 0.978228658146057 }, {
// -0.978228658146057, 0.730152005574049, -0.978228658146057 }, {
//-0.978228658146057, 0.730152005574049, -0.887062599768095 }, {
//-0.978228658146057, 0.730152005574049, -0.730152005574049 }, {
//-0.978228658146057, 0.730152005574049, -0.519096129206812 }, {
//-0.978228658146057, 0.730152005574049, -0.269543155952345 }, {
//-0.978228658146057, 0.730152005574049, 0 }, { -0.978228658146057,
// 0.730152005574049, 0.269543155952345 }, { -0.978228658146057,
// 0.730152005574049, 0.519096129206812 }, { -0.978228658146057,
// 0.730152005574049, 0.730152005574049 }, { -0.978228658146057,
// 0.730152005574049, 0.887062599768095 }, { -0.978228658146057,
// 0.730152005574049, 0.978228658146057 }, { -0.978228658146057,
// 0.887062599768095, -0.978228658146057 }, { -0.978228658146057,
// 0.887062599768095, -0.887062599768095 }, { -0.978228658146057,
// 0.887062599768095, -0.730152005574049 }, { -0.978228658146057,
// 0.887062599768095, -0.519096129206812 }, { -0.978228658146057,
// 0.887062599768095, -0.269543155952345 }, { -0.978228658146057,
// 0.887062599768095, 0 }, { -0.978228658146057, 0.887062599768095,
// 0.269543155952345 }, { -0.978228658146057, 0.887062599768095,
// 0.519096129206812 }, { -0.978228658146057, 0.887062599768095,
// 0.730152005574049 }, { -0.978228658146057, 0.887062599768095,
// 0.887062599768095 }, { -0.978228658146057, 0.887062599768095,
// 0.978228658146057 }, { -0.978228658146057, 0.978228658146057,
// -0.978228658146057 }, { -0.978228658146057, 0.978228658146057,
// -0.887062599768095 }, { -0.978228658146057, 0.978228658146057,
// -0.730152005574049 }, { -0.978228658146057, 0.978228658146057,
// -0.519096129206812 }, { -0.978228658146057, 0.978228658146057,
// -0.269543155952345 }, { -0.978228658146057, 0.978228658146057, 0 }, {
// -0.978228658146057, 0.978228658146057, 0.269543155952345 }, {
// -0.978228658146057, 0.978228658146057, 0.519096129206812 }, {
// -0.978228658146057, 0.978228658146057, 0.730152005574049 }, {
// -0.978228658146057, 0.978228658146057, 0.887062599768095 }, {
// -0.978228658146057, 0.978228658146057, 0.978228658146057 }, {
// -0.887062599768095, -0.978228658146057, -0.978228658146057 }, {
//-0.887062599768095, -0.978228658146057, -0.887062599768095 }, {
//-0.887062599768095, -0.978228658146057, -0.730152005574049 }, {
//-0.887062599768095, -0.978228658146057, -0.519096129206812 }, {
//-0.887062599768095, -0.978228658146057, -0.269543155952345 }, {
//-0.887062599768095, -0.978228658146057, 0 }, { -0.887062599768095,
//-0.978228658146057, 0.269543155952345 }, { -0.887062599768095,
//-0.978228658146057, 0.519096129206812 }, { -0.887062599768095,
//-0.978228658146057, 0.730152005574049 }, { -0.887062599768095,
//-0.978228658146057, 0.887062599768095 }, { -0.887062599768095,
//-0.978228658146057, 0.978228658146057 }, { -0.887062599768095,
//-0.887062599768095, -0.978228658146057 }, { -0.887062599768095,
//-0.887062599768095, -0.887062599768095 }, { -0.887062599768095,
//-0.887062599768095, -0.730152005574049 }, { -0.887062599768095,
//-0.887062599768095, -0.519096129206812 }, { -0.887062599768095,
//-0.887062599768095, -0.269543155952345 }, { -0.887062599768095,
//-0.887062599768095, 0 }, { -0.887062599768095, -0.887062599768095,
// 0.269543155952345 }, { -0.887062599768095, -0.887062599768095,
// 0.519096129206812 }, { -0.887062599768095, -0.887062599768095,
// 0.730152005574049 }, { -0.887062599768095, -0.887062599768095,
// 0.887062599768095 }, { -0.887062599768095, -0.887062599768095,
// 0.978228658146057 }, { -0.887062599768095, -0.730152005574049,
//-0.978228658146057 }, { -0.887062599768095, -0.730152005574049,
//-0.887062599768095 }, { -0.887062599768095, -0.730152005574049,
//-0.730152005574049 }, { -0.887062599768095, -0.730152005574049,
//-0.519096129206812 }, { -0.887062599768095, -0.730152005574049,
//-0.269543155952345 }, { -0.887062599768095, -0.730152005574049, 0 }, {
//-0.887062599768095, -0.730152005574049, 0.269543155952345 }, {
// -0.887062599768095, -0.730152005574049, 0.519096129206812 }, {
// -0.887062599768095, -0.730152005574049, 0.730152005574049 }, {
// -0.887062599768095, -0.730152005574049, 0.887062599768095 }, {
// -0.887062599768095, -0.730152005574049, 0.978228658146057 }, {
// -0.887062599768095, -0.519096129206812, -0.978228658146057 }, {
//-0.887062599768095, -0.519096129206812, -0.887062599768095 }, {
//-0.887062599768095, -0.519096129206812, -0.730152005574049 }, {
//-0.887062599768095, -0.519096129206812, -0.519096129206812 }, {
//-0.887062599768095, -0.519096129206812, -0.269543155952345 }, {
//-0.887062599768095, -0.519096129206812, 0 }, { -0.887062599768095,
//-0.519096129206812, 0.269543155952345 }, { -0.887062599768095,
//-0.519096129206812, 0.519096129206812 }, { -0.887062599768095,
//-0.519096129206812, 0.730152005574049 }, { -0.887062599768095,
//-0.519096129206812, 0.887062599768095 }, { -0.887062599768095,
//-0.519096129206812, 0.978228658146057 }, { -0.887062599768095,
//-0.269543155952345, -0.978228658146057 }, { -0.887062599768095,
//-0.269543155952345, -0.887062599768095 }, { -0.887062599768095,
//-0.269543155952345, -0.730152005574049 }, { -0.887062599768095,
//-0.269543155952345, -0.519096129206812 }, { -0.887062599768095,
//-0.269543155952345, -0.269543155952345 }, { -0.887062599768095,
//-0.269543155952345, 0 }, { -0.887062599768095, -0.269543155952345,
// 0.269543155952345 }, { -0.887062599768095, -0.269543155952345,
// 0.519096129206812 }, { -0.887062599768095, -0.269543155952345,
// 0.730152005574049 }, { -0.887062599768095, -0.269543155952345,
// 0.887062599768095 }, { -0.887062599768095, -0.269543155952345,
// 0.978228658146057 }, { -0.887062599768095, 0, -0.978228658146057 }, {
//-0.887062599768095, 0, -0.887062599768095 }, { -0.887062599768095, 0,
//-0.730152005574049 }, { -0.887062599768095, 0, -0.519096129206812 }, {
//-0.887062599768095, 0, -0.269543155952345 }, { -0.887062599768095, 0, 0 }, {
//-0.887062599768095, 0, 0.269543155952345 }, { -0.887062599768095, 0,
// 0.519096129206812 }, { -0.887062599768095, 0, 0.730152005574049 }, {
// -0.887062599768095, 0, 0.887062599768095 }, { -0.887062599768095, 0,
// 0.978228658146057 }, { -0.887062599768095, 0.269543155952345,
// -0.978228658146057 }, { -0.887062599768095, 0.269543155952345,
// -0.887062599768095 }, { -0.887062599768095, 0.269543155952345,
// -0.730152005574049 }, { -0.887062599768095, 0.269543155952345,
// -0.519096129206812 }, { -0.887062599768095, 0.269543155952345,
// -0.269543155952345 }, { -0.887062599768095, 0.269543155952345, 0 }, {
// -0.887062599768095, 0.269543155952345, 0.269543155952345 }, {
// -0.887062599768095, 0.269543155952345, 0.519096129206812 }, {
// -0.887062599768095, 0.269543155952345, 0.730152005574049 }, {
// -0.887062599768095, 0.269543155952345, 0.887062599768095 }, {
// -0.887062599768095, 0.269543155952345, 0.978228658146057 }, {
// -0.887062599768095, 0.519096129206812, -0.978228658146057 }, {
//-0.887062599768095, 0.519096129206812, -0.887062599768095 }, {
//-0.887062599768095, 0.519096129206812, -0.730152005574049 }, {
//-0.887062599768095, 0.519096129206812, -0.519096129206812 }, {
//-0.887062599768095, 0.519096129206812, -0.269543155952345 }, {
//-0.887062599768095, 0.519096129206812, 0 }, { -0.887062599768095,
// 0.519096129206812, 0.269543155952345 }, { -0.887062599768095,
// 0.519096129206812, 0.519096129206812 }, { -0.887062599768095,
// 0.519096129206812, 0.730152005574049 }, { -0.887062599768095,
// 0.519096129206812, 0.887062599768095 }, { -0.887062599768095,
// 0.519096129206812, 0.978228658146057 }, { -0.887062599768095,
// 0.730152005574049, -0.978228658146057 }, { -0.887062599768095,
// 0.730152005574049, -0.887062599768095 }, { -0.887062599768095,
// 0.730152005574049, -0.730152005574049 }, { -0.887062599768095,
// 0.730152005574049, -0.519096129206812 }, { -0.887062599768095,
// 0.730152005574049, -0.269543155952345 }, { -0.887062599768095,
// 0.730152005574049, 0 }, { -0.887062599768095, 0.730152005574049,
// 0.269543155952345 }, { -0.887062599768095, 0.730152005574049,
// 0.519096129206812 }, { -0.887062599768095, 0.730152005574049,
// 0.730152005574049 }, { -0.887062599768095, 0.730152005574049,
// 0.887062599768095 }, { -0.887062599768095, 0.730152005574049,
// 0.978228658146057 }, { -0.887062599768095, 0.887062599768095,
// -0.978228658146057 }, { -0.887062599768095, 0.887062599768095,
// -0.887062599768095 }, { -0.887062599768095, 0.887062599768095,
// -0.730152005574049 }, { -0.887062599768095, 0.887062599768095,
// -0.519096129206812 }, { -0.887062599768095, 0.887062599768095,
// -0.269543155952345 }, { -0.887062599768095, 0.887062599768095, 0 }, {
// -0.887062599768095, 0.887062599768095, 0.269543155952345 }, {
// -0.887062599768095, 0.887062599768095, 0.519096129206812 }, {
// -0.887062599768095, 0.887062599768095, 0.730152005574049 }, {
// -0.887062599768095, 0.887062599768095, 0.887062599768095 }, {
// -0.887062599768095, 0.887062599768095, 0.978228658146057 }, {
// -0.887062599768095, 0.978228658146057, -0.978228658146057 }, {
//-0.887062599768095, 0.978228658146057, -0.887062599768095 }, {
//-0.887062599768095, 0.978228658146057, -0.730152005574049 }, {
//-0.887062599768095, 0.978228658146057, -0.519096129206812 }, {
//-0.887062599768095, 0.978228658146057, -0.269543155952345 }, {
//-0.887062599768095, 0.978228658146057, 0 }, { -0.887062599768095,
// 0.978228658146057, 0.269543155952345 }, { -0.887062599768095,
// 0.978228658146057, 0.519096129206812 }, { -0.887062599768095,
// 0.978228658146057, 0.730152005574049 }, { -0.887062599768095,
// 0.978228658146057, 0.887062599768095 }, { -0.887062599768095,
// 0.978228658146057, 0.978228658146057 }, { -0.730152005574049,
//-0.978228658146057, -0.978228658146057 }, { -0.730152005574049,
//-0.978228658146057, -0.887062599768095 }, { -0.730152005574049,
//-0.978228658146057, -0.730152005574049 }, { -0.730152005574049,
//-0.978228658146057, -0.519096129206812 }, { -0.730152005574049,
//-0.978228658146057, -0.269543155952345 }, { -0.730152005574049,
//-0.978228658146057, 0 }, { -0.730152005574049, -0.978228658146057,
// 0.269543155952345 }, { -0.730152005574049, -0.978228658146057,
// 0.519096129206812 }, { -0.730152005574049, -0.978228658146057,
// 0.730152005574049 }, { -0.730152005574049, -0.978228658146057,
// 0.887062599768095 }, { -0.730152005574049, -0.978228658146057,
// 0.978228658146057 }, { -0.730152005574049, -0.887062599768095,
//-0.978228658146057 }, { -0.730152005574049, -0.887062599768095,
//-0.887062599768095 }, { -0.730152005574049, -0.887062599768095,
//-0.730152005574049 }, { -0.730152005574049, -0.887062599768095,
//-0.519096129206812 }, { -0.730152005574049, -0.887062599768095,
//-0.269543155952345 }, { -0.730152005574049, -0.887062599768095, 0 }, {
//-0.730152005574049, -0.887062599768095, 0.269543155952345 }, {
// -0.730152005574049, -0.887062599768095, 0.519096129206812 }, {
// -0.730152005574049, -0.887062599768095, 0.730152005574049 }, {
// -0.730152005574049, -0.887062599768095, 0.887062599768095 }, {
// -0.730152005574049, -0.887062599768095, 0.978228658146057 }, {
// -0.730152005574049, -0.730152005574049, -0.978228658146057 }, {
//-0.730152005574049, -0.730152005574049, -0.887062599768095 }, {
//-0.730152005574049, -0.730152005574049, -0.730152005574049 }, {
//-0.730152005574049, -0.730152005574049, -0.519096129206812 }, {
//-0.730152005574049, -0.730152005574049, -0.269543155952345 }, {
//-0.730152005574049, -0.730152005574049, 0 }, { -0.730152005574049,
//-0.730152005574049, 0.269543155952345 }, { -0.730152005574049,
//-0.730152005574049, 0.519096129206812 }, { -0.730152005574049,
//-0.730152005574049, 0.730152005574049 }, { -0.730152005574049,
//-0.730152005574049, 0.887062599768095 }, { -0.730152005574049,
//-0.730152005574049, 0.978228658146057 }, { -0.730152005574049,
//-0.519096129206812, -0.978228658146057 }, { -0.730152005574049,
//-0.519096129206812, -0.887062599768095 }, { -0.730152005574049,
//-0.519096129206812, -0.730152005574049 }, { -0.730152005574049,
//-0.519096129206812, -0.519096129206812 }, { -0.730152005574049,
//-0.519096129206812, -0.269543155952345 }, { -0.730152005574049,
//-0.519096129206812, 0 }, { -0.730152005574049, -0.519096129206812,
// 0.269543155952345 }, { -0.730152005574049, -0.519096129206812,
// 0.519096129206812 }, { -0.730152005574049, -0.519096129206812,
// 0.730152005574049 }, { -0.730152005574049, -0.519096129206812,
// 0.887062599768095 }, { -0.730152005574049, -0.519096129206812,
// 0.978228658146057 }, { -0.730152005574049, -0.269543155952345,
//-0.978228658146057 }, { -0.730152005574049, -0.269543155952345,
//-0.887062599768095 }, { -0.730152005574049, -0.269543155952345,
//-0.730152005574049 }, { -0.730152005574049, -0.269543155952345,
//-0.519096129206812 }, { -0.730152005574049, -0.269543155952345,
//-0.269543155952345 }, { -0.730152005574049, -0.269543155952345, 0 }, {
//-0.730152005574049, -0.269543155952345, 0.269543155952345 }, {
// -0.730152005574049, -0.269543155952345, 0.519096129206812 }, {
// -0.730152005574049, -0.269543155952345, 0.730152005574049 }, {
// -0.730152005574049, -0.269543155952345, 0.887062599768095 }, {
// -0.730152005574049, -0.269543155952345, 0.978228658146057 }, {
// -0.730152005574049, 0, -0.978228658146057 }, { -0.730152005574049, 0,
//-0.887062599768095 }, { -0.730152005574049, 0, -0.730152005574049 }, {
//-0.730152005574049, 0, -0.519096129206812 }, { -0.730152005574049, 0,
//-0.269543155952345 }, { -0.730152005574049, 0, 0 }, { -0.730152005574049, 0,
// 0.269543155952345 }, { -0.730152005574049, 0, 0.519096129206812 }, {
// -0.730152005574049, 0, 0.730152005574049 }, { -0.730152005574049, 0,
// 0.887062599768095 }, { -0.730152005574049, 0, 0.978228658146057 }, {
// -0.730152005574049, 0.269543155952345, -0.978228658146057 }, {
//-0.730152005574049, 0.269543155952345, -0.887062599768095 }, {
//-0.730152005574049, 0.269543155952345, -0.730152005574049 }, {
//-0.730152005574049, 0.269543155952345, -0.519096129206812 }, {
//-0.730152005574049, 0.269543155952345, -0.269543155952345 }, {
//-0.730152005574049, 0.269543155952345, 0 }, { -0.730152005574049,
// 0.269543155952345, 0.269543155952345 }, { -0.730152005574049,
// 0.269543155952345, 0.519096129206812 }, { -0.730152005574049,
// 0.269543155952345, 0.730152005574049 }, { -0.730152005574049,
// 0.269543155952345, 0.887062599768095 }, { -0.730152005574049,
// 0.269543155952345, 0.978228658146057 }, { -0.730152005574049,
// 0.519096129206812, -0.978228658146057 }, { -0.730152005574049,
// 0.519096129206812, -0.887062599768095 }, { -0.730152005574049,
// 0.519096129206812, -0.730152005574049 }, { -0.730152005574049,
// 0.519096129206812, -0.519096129206812 }, { -0.730152005574049,
// 0.519096129206812, -0.269543155952345 }, { -0.730152005574049,
// 0.519096129206812, 0 }, { -0.730152005574049, 0.519096129206812,
// 0.269543155952345 }, { -0.730152005574049, 0.519096129206812,
// 0.519096129206812 }, { -0.730152005574049, 0.519096129206812,
// 0.730152005574049 }, { -0.730152005574049, 0.519096129206812,
// 0.887062599768095 }, { -0.730152005574049, 0.519096129206812,
// 0.978228658146057 }, { -0.730152005574049, 0.730152005574049,
// -0.978228658146057 }, { -0.730152005574049, 0.730152005574049,
// -0.887062599768095 }, { -0.730152005574049, 0.730152005574049,
// -0.730152005574049 }, { -0.730152005574049, 0.730152005574049,
// -0.519096129206812 }, { -0.730152005574049, 0.730152005574049,
// -0.269543155952345 }, { -0.730152005574049, 0.730152005574049, 0 }, {
// -0.730152005574049, 0.730152005574049, 0.269543155952345 }, {
// -0.730152005574049, 0.730152005574049, 0.519096129206812 }, {
// -0.730152005574049, 0.730152005574049, 0.730152005574049 }, {
// -0.730152005574049, 0.730152005574049, 0.887062599768095 }, {
// -0.730152005574049, 0.730152005574049, 0.978228658146057 }, {
// -0.730152005574049, 0.887062599768095, -0.978228658146057 }, {
//-0.730152005574049, 0.887062599768095, -0.887062599768095 }, {
//-0.730152005574049, 0.887062599768095, -0.730152005574049 }, {
//-0.730152005574049, 0.887062599768095, -0.519096129206812 }, {
//-0.730152005574049, 0.887062599768095, -0.269543155952345 }, {
//-0.730152005574049, 0.887062599768095, 0 }, { -0.730152005574049,
// 0.887062599768095, 0.269543155952345 }, { -0.730152005574049,
// 0.887062599768095, 0.519096129206812 }, { -0.730152005574049,
// 0.887062599768095, 0.730152005574049 }, { -0.730152005574049,
// 0.887062599768095, 0.887062599768095 }, { -0.730152005574049,
// 0.887062599768095, 0.978228658146057 }, { -0.730152005574049,
// 0.978228658146057, -0.978228658146057 }, { -0.730152005574049,
// 0.978228658146057, -0.887062599768095 }, { -0.730152005574049,
// 0.978228658146057, -0.730152005574049 }, { -0.730152005574049,
// 0.978228658146057, -0.519096129206812 }, { -0.730152005574049,
// 0.978228658146057, -0.269543155952345 }, { -0.730152005574049,
// 0.978228658146057, 0 }, { -0.730152005574049, 0.978228658146057,
// 0.269543155952345 }, { -0.730152005574049, 0.978228658146057,
// 0.519096129206812 }, { -0.730152005574049, 0.978228658146057,
// 0.730152005574049 }, { -0.730152005574049, 0.978228658146057,
// 0.887062599768095 }, { -0.730152005574049, 0.978228658146057,
// 0.978228658146057 }, { -0.519096129206812, -0.978228658146057,
//-0.978228658146057 }, { -0.519096129206812, -0.978228658146057,
//-0.887062599768095 }, { -0.519096129206812, -0.978228658146057,
//-0.730152005574049 }, { -0.519096129206812, -0.978228658146057,
//-0.519096129206812 }, { -0.519096129206812, -0.978228658146057,
//-0.269543155952345 }, { -0.519096129206812, -0.978228658146057, 0 }, {
//-0.519096129206812, -0.978228658146057, 0.269543155952345 }, {
// -0.519096129206812, -0.978228658146057, 0.519096129206812 }, {
// -0.519096129206812, -0.978228658146057, 0.730152005574049 }, {
// -0.519096129206812, -0.978228658146057, 0.887062599768095 }, {
// -0.519096129206812, -0.978228658146057, 0.978228658146057 }, {
// -0.519096129206812, -0.887062599768095, -0.978228658146057 }, {
//-0.519096129206812, -0.887062599768095, -0.887062599768095 }, {
//-0.519096129206812, -0.887062599768095, -0.730152005574049 }, {
//-0.519096129206812, -0.887062599768095, -0.519096129206812 }, {
//-0.519096129206812, -0.887062599768095, -0.269543155952345 }, {
//-0.519096129206812, -0.887062599768095, 0 }, { -0.519096129206812,
//-0.887062599768095, 0.269543155952345 }, { -0.519096129206812,
//-0.887062599768095, 0.519096129206812 }, { -0.519096129206812,
//-0.887062599768095, 0.730152005574049 }, { -0.519096129206812,
//-0.887062599768095, 0.887062599768095 }, { -0.519096129206812,
//-0.887062599768095, 0.978228658146057 }, { -0.519096129206812,
//-0.730152005574049, -0.978228658146057 }, { -0.519096129206812,
//-0.730152005574049, -0.887062599768095 }, { -0.519096129206812,
//-0.730152005574049, -0.730152005574049 }, { -0.519096129206812,
//-0.730152005574049, -0.519096129206812 }, { -0.519096129206812,
//-0.730152005574049, -0.269543155952345 }, { -0.519096129206812,
//-0.730152005574049, 0 }, { -0.519096129206812, -0.730152005574049,
// 0.269543155952345 }, { -0.519096129206812, -0.730152005574049,
// 0.519096129206812 }, { -0.519096129206812, -0.730152005574049,
// 0.730152005574049 }, { -0.519096129206812, -0.730152005574049,
// 0.887062599768095 }, { -0.519096129206812, -0.730152005574049,
// 0.978228658146057 }, { -0.519096129206812, -0.519096129206812,
//-0.978228658146057 }, { -0.519096129206812, -0.519096129206812,
//-0.887062599768095 }, { -0.519096129206812, -0.519096129206812,
//-0.730152005574049 }, { -0.519096129206812, -0.519096129206812,
//-0.519096129206812 }, { -0.519096129206812, -0.519096129206812,
//-0.269543155952345 }, { -0.519096129206812, -0.519096129206812, 0 }, {
//-0.519096129206812, -0.519096129206812, 0.269543155952345 }, {
// -0.519096129206812, -0.519096129206812, 0.519096129206812 }, {
// -0.519096129206812, -0.519096129206812, 0.730152005574049 }, {
// -0.519096129206812, -0.519096129206812, 0.887062599768095 }, {
// -0.519096129206812, -0.519096129206812, 0.978228658146057 }, {
// -0.519096129206812, -0.269543155952345, -0.978228658146057 }, {
//-0.519096129206812, -0.269543155952345, -0.887062599768095 }, {
//-0.519096129206812, -0.269543155952345, -0.730152005574049 }, {
//-0.519096129206812, -0.269543155952345, -0.519096129206812 }, {
//-0.519096129206812, -0.269543155952345, -0.269543155952345 }, {
//-0.519096129206812, -0.269543155952345, 0 }, { -0.519096129206812,
//-0.269543155952345, 0.269543155952345 }, { -0.519096129206812,
//-0.269543155952345, 0.519096129206812 }, { -0.519096129206812,
//-0.269543155952345, 0.730152005574049 }, { -0.519096129206812,
//-0.269543155952345, 0.887062599768095 }, { -0.519096129206812,
//-0.269543155952345, 0.978228658146057 }, { -0.519096129206812, 0,
//-0.978228658146057 }, { -0.519096129206812, 0, -0.887062599768095 }, {
//-0.519096129206812, 0, -0.730152005574049 }, { -0.519096129206812, 0,
//-0.519096129206812 }, { -0.519096129206812, 0, -0.269543155952345 }, {
//-0.519096129206812, 0, 0 }, { -0.519096129206812, 0, 0.269543155952345 }, {
// -0.519096129206812, 0, 0.519096129206812 }, { -0.519096129206812, 0,
// 0.730152005574049 }, { -0.519096129206812, 0, 0.887062599768095 }, {
// -0.519096129206812, 0, 0.978228658146057 }, { -0.519096129206812,
// 0.269543155952345, -0.978228658146057 }, { -0.519096129206812,
// 0.269543155952345, -0.887062599768095 }, { -0.519096129206812,
// 0.269543155952345, -0.730152005574049 }, { -0.519096129206812,
// 0.269543155952345, -0.519096129206812 }, { -0.519096129206812,
// 0.269543155952345, -0.269543155952345 }, { -0.519096129206812,
// 0.269543155952345, 0 }, { -0.519096129206812, 0.269543155952345,
// 0.269543155952345 }, { -0.519096129206812, 0.269543155952345,
// 0.519096129206812 }, { -0.519096129206812, 0.269543155952345,
// 0.730152005574049 }, { -0.519096129206812, 0.269543155952345,
// 0.887062599768095 }, { -0.519096129206812, 0.269543155952345,
// 0.978228658146057 }, { -0.519096129206812, 0.519096129206812,
// -0.978228658146057 }, { -0.519096129206812, 0.519096129206812,
// -0.887062599768095 }, { -0.519096129206812, 0.519096129206812,
// -0.730152005574049 }, { -0.519096129206812, 0.519096129206812,
// -0.519096129206812 }, { -0.519096129206812, 0.519096129206812,
// -0.269543155952345 }, { -0.519096129206812, 0.519096129206812, 0 }, {
// -0.519096129206812, 0.519096129206812, 0.269543155952345 }, {
// -0.519096129206812, 0.519096129206812, 0.519096129206812 }, {
// -0.519096129206812, 0.519096129206812, 0.730152005574049 }, {
// -0.519096129206812, 0.519096129206812, 0.887062599768095 }, {
// -0.519096129206812, 0.519096129206812, 0.978228658146057 }, {
// -0.519096129206812, 0.730152005574049, -0.978228658146057 }, {
//-0.519096129206812, 0.730152005574049, -0.887062599768095 }, {
//-0.519096129206812, 0.730152005574049, -0.730152005574049 }, {
//-0.519096129206812, 0.730152005574049, -0.519096129206812 }, {
//-0.519096129206812, 0.730152005574049, -0.269543155952345 }, {
//-0.519096129206812, 0.730152005574049, 0 }, { -0.519096129206812,
// 0.730152005574049, 0.269543155952345 }, { -0.519096129206812,
// 0.730152005574049, 0.519096129206812 }, { -0.519096129206812,
// 0.730152005574049, 0.730152005574049 }, { -0.519096129206812,
// 0.730152005574049, 0.887062599768095 }, { -0.519096129206812,
// 0.730152005574049, 0.978228658146057 }, { -0.519096129206812,
// 0.887062599768095, -0.978228658146057 }, { -0.519096129206812,
// 0.887062599768095, -0.887062599768095 }, { -0.519096129206812,
// 0.887062599768095, -0.730152005574049 }, { -0.519096129206812,
// 0.887062599768095, -0.519096129206812 }, { -0.519096129206812,
// 0.887062599768095, -0.269543155952345 }, { -0.519096129206812,
// 0.887062599768095, 0 }, { -0.519096129206812, 0.887062599768095,
// 0.269543155952345 }, { -0.519096129206812, 0.887062599768095,
// 0.519096129206812 }, { -0.519096129206812, 0.887062599768095,
// 0.730152005574049 }, { -0.519096129206812, 0.887062599768095,
// 0.887062599768095 }, { -0.519096129206812, 0.887062599768095,
// 0.978228658146057 }, { -0.519096129206812, 0.978228658146057,
// -0.978228658146057 }, { -0.519096129206812, 0.978228658146057,
// -0.887062599768095 }, { -0.519096129206812, 0.978228658146057,
// -0.730152005574049 }, { -0.519096129206812, 0.978228658146057,
// -0.519096129206812 }, { -0.519096129206812, 0.978228658146057,
// -0.269543155952345 }, { -0.519096129206812, 0.978228658146057, 0 }, {
// -0.519096129206812, 0.978228658146057, 0.269543155952345 }, {
// -0.519096129206812, 0.978228658146057, 0.519096129206812 }, {
// -0.519096129206812, 0.978228658146057, 0.730152005574049 }, {
// -0.519096129206812, 0.978228658146057, 0.887062599768095 }, {
// -0.519096129206812, 0.978228658146057, 0.978228658146057 }, {
// -0.269543155952345, -0.978228658146057, -0.978228658146057 }, {
//-0.269543155952345, -0.978228658146057, -0.887062599768095 }, {
//-0.269543155952345, -0.978228658146057, -0.730152005574049 }, {
//-0.269543155952345, -0.978228658146057, -0.519096129206812 }, {
//-0.269543155952345, -0.978228658146057, -0.269543155952345 }, {
//-0.269543155952345, -0.978228658146057, 0 }, { -0.269543155952345,
//-0.978228658146057, 0.269543155952345 }, { -0.269543155952345,
//-0.978228658146057, 0.519096129206812 }, { -0.269543155952345,
//-0.978228658146057, 0.730152005574049 }, { -0.269543155952345,
//-0.978228658146057, 0.887062599768095 }, { -0.269543155952345,
//-0.978228658146057, 0.978228658146057 }, { -0.269543155952345,
//-0.887062599768095, -0.978228658146057 }, { -0.269543155952345,
//-0.887062599768095, -0.887062599768095 }, { -0.269543155952345,
//-0.887062599768095, -0.730152005574049 }, { -0.269543155952345,
//-0.887062599768095, -0.519096129206812 }, { -0.269543155952345,
//-0.887062599768095, -0.269543155952345 }, { -0.269543155952345,
//-0.887062599768095, 0 }, { -0.269543155952345, -0.887062599768095,
// 0.269543155952345 }, { -0.269543155952345, -0.887062599768095,
// 0.519096129206812 }, { -0.269543155952345, -0.887062599768095,
// 0.730152005574049 }, { -0.269543155952345, -0.887062599768095,
// 0.887062599768095 }, { -0.269543155952345, -0.887062599768095,
// 0.978228658146057 }, { -0.269543155952345, -0.730152005574049,
//-0.978228658146057 }, { -0.269543155952345, -0.730152005574049,
//-0.887062599768095 }, { -0.269543155952345, -0.730152005574049,
//-0.730152005574049 }, { -0.269543155952345, -0.730152005574049,
//-0.519096129206812 }, { -0.269543155952345, -0.730152005574049,
//-0.269543155952345 }, { -0.269543155952345, -0.730152005574049, 0 }, {
//-0.269543155952345, -0.730152005574049, 0.269543155952345 }, {
// -0.269543155952345, -0.730152005574049, 0.519096129206812 }, {
// -0.269543155952345, -0.730152005574049, 0.730152005574049 }, {
// -0.269543155952345, -0.730152005574049, 0.887062599768095 }, {
// -0.269543155952345, -0.730152005574049, 0.978228658146057 }, {
// -0.269543155952345, -0.519096129206812, -0.978228658146057 }, {
//-0.269543155952345, -0.519096129206812, -0.887062599768095 }, {
//-0.269543155952345, -0.519096129206812, -0.730152005574049 }, {
//-0.269543155952345, -0.519096129206812, -0.519096129206812 }, {
//-0.269543155952345, -0.519096129206812, -0.269543155952345 }, {
//-0.269543155952345, -0.519096129206812, 0 }, { -0.269543155952345,
//-0.519096129206812, 0.269543155952345 }, { -0.269543155952345,
//-0.519096129206812, 0.519096129206812 }, { -0.269543155952345,
//-0.519096129206812, 0.730152005574049 }, { -0.269543155952345,
//-0.519096129206812, 0.887062599768095 }, { -0.269543155952345,
//-0.519096129206812, 0.978228658146057 }, { -0.269543155952345,
//-0.269543155952345, -0.978228658146057 }, { -0.269543155952345,
//-0.269543155952345, -0.887062599768095 }, { -0.269543155952345,
//-0.269543155952345, -0.730152005574049 }, { -0.269543155952345,
//-0.269543155952345, -0.519096129206812 }, { -0.269543155952345,
//-0.269543155952345, -0.269543155952345 }, { -0.269543155952345,
//-0.269543155952345, 0 }, { -0.269543155952345, -0.269543155952345,
// 0.269543155952345 }, { -0.269543155952345, -0.269543155952345,
// 0.519096129206812 }, { -0.269543155952345, -0.269543155952345,
// 0.730152005574049 }, { -0.269543155952345, -0.269543155952345,
// 0.887062599768095 }, { -0.269543155952345, -0.269543155952345,
// 0.978228658146057 }, { -0.269543155952345, 0, -0.978228658146057 }, {
//-0.269543155952345, 0, -0.887062599768095 }, { -0.269543155952345, 0,
//-0.730152005574049 }, { -0.269543155952345, 0, -0.519096129206812 }, {
//-0.269543155952345, 0, -0.269543155952345 }, { -0.269543155952345, 0, 0 }, {
//-0.269543155952345, 0, 0.269543155952345 }, { -0.269543155952345, 0,
// 0.519096129206812 }, { -0.269543155952345, 0, 0.730152005574049 }, {
// -0.269543155952345, 0, 0.887062599768095 }, { -0.269543155952345, 0,
// 0.978228658146057 }, { -0.269543155952345, 0.269543155952345,
// -0.978228658146057 }, { -0.269543155952345, 0.269543155952345,
// -0.887062599768095 }, { -0.269543155952345, 0.269543155952345,
// -0.730152005574049 }, { -0.269543155952345, 0.269543155952345,
// -0.519096129206812 }, { -0.269543155952345, 0.269543155952345,
// -0.269543155952345 }, { -0.269543155952345, 0.269543155952345, 0 }, {
// -0.269543155952345, 0.269543155952345, 0.269543155952345 }, {
// -0.269543155952345, 0.269543155952345, 0.519096129206812 }, {
// -0.269543155952345, 0.269543155952345, 0.730152005574049 }, {
// -0.269543155952345, 0.269543155952345, 0.887062599768095 }, {
// -0.269543155952345, 0.269543155952345, 0.978228658146057 }, {
// -0.269543155952345, 0.519096129206812, -0.978228658146057 }, {
//-0.269543155952345, 0.519096129206812, -0.887062599768095 }, {
//-0.269543155952345, 0.519096129206812, -0.730152005574049 }, {
//-0.269543155952345, 0.519096129206812, -0.519096129206812 }, {
//-0.269543155952345, 0.519096129206812, -0.269543155952345 }, {
//-0.269543155952345, 0.519096129206812, 0 }, { -0.269543155952345,
// 0.519096129206812, 0.269543155952345 }, { -0.269543155952345,
// 0.519096129206812, 0.519096129206812 }, { -0.269543155952345,
// 0.519096129206812, 0.730152005574049 }, { -0.269543155952345,
// 0.519096129206812, 0.887062599768095 }, { -0.269543155952345,
// 0.519096129206812, 0.978228658146057 }, { -0.269543155952345,
// 0.730152005574049, -0.978228658146057 }, { -0.269543155952345,
// 0.730152005574049, -0.887062599768095 }, { -0.269543155952345,
// 0.730152005574049, -0.730152005574049 }, { -0.269543155952345,
// 0.730152005574049, -0.519096129206812 }, { -0.269543155952345,
// 0.730152005574049, -0.269543155952345 }, { -0.269543155952345,
// 0.730152005574049, 0 }, { -0.269543155952345, 0.730152005574049,
// 0.269543155952345 }, { -0.269543155952345, 0.730152005574049,
// 0.519096129206812 }, { -0.269543155952345, 0.730152005574049,
// 0.730152005574049 }, { -0.269543155952345, 0.730152005574049,
// 0.887062599768095 }, { -0.269543155952345, 0.730152005574049,
// 0.978228658146057 }, { -0.269543155952345, 0.887062599768095,
// -0.978228658146057 }, { -0.269543155952345, 0.887062599768095,
// -0.887062599768095 }, { -0.269543155952345, 0.887062599768095,
// -0.730152005574049 }, { -0.269543155952345, 0.887062599768095,
// -0.519096129206812 }, { -0.269543155952345, 0.887062599768095,
// -0.269543155952345 }, { -0.269543155952345, 0.887062599768095, 0 }, {
// -0.269543155952345, 0.887062599768095, 0.269543155952345 }, {
// -0.269543155952345, 0.887062599768095, 0.519096129206812 }, {
// -0.269543155952345, 0.887062599768095, 0.730152005574049 }, {
// -0.269543155952345, 0.887062599768095, 0.887062599768095 }, {
// -0.269543155952345, 0.887062599768095, 0.978228658146057 }, {
// -0.269543155952345, 0.978228658146057, -0.978228658146057 }, {
//-0.269543155952345, 0.978228658146057, -0.887062599768095 }, {
//-0.269543155952345, 0.978228658146057, -0.730152005574049 }, {
//-0.269543155952345, 0.978228658146057, -0.519096129206812 }, {
//-0.269543155952345, 0.978228658146057, -0.269543155952345 }, {
//-0.269543155952345, 0.978228658146057, 0 }, { -0.269543155952345,
// 0.978228658146057, 0.269543155952345 }, { -0.269543155952345,
// 0.978228658146057, 0.519096129206812 }, { -0.269543155952345,
// 0.978228658146057, 0.730152005574049 }, { -0.269543155952345,
// 0.978228658146057, 0.887062599768095 }, { -0.269543155952345,
// 0.978228658146057, 0.978228658146057 }, { 0, -0.978228658146057,
//-0.978228658146057 }, { 0, -0.978228658146057, -0.887062599768095 }, { 0,
//-0.978228658146057, -0.730152005574049 }, { 0, -0.978228658146057,
//-0.519096129206812 }, { 0, -0.978228658146057, -0.269543155952345 }, { 0,
//-0.978228658146057, 0 }, { 0, -0.978228658146057, 0.269543155952345 }, { 0,
// -0.978228658146057, 0.519096129206812 }, { 0, -0.978228658146057,
// 0.730152005574049 }, { 0, -0.978228658146057, 0.887062599768095 }, { 0,
// -0.978228658146057, 0.978228658146057 }, { 0, -0.887062599768095,
//-0.978228658146057 }, { 0, -0.887062599768095, -0.887062599768095 }, { 0,
//-0.887062599768095, -0.730152005574049 }, { 0, -0.887062599768095,
//-0.519096129206812 }, { 0, -0.887062599768095, -0.269543155952345 }, { 0,
//-0.887062599768095, 0 }, { 0, -0.887062599768095, 0.269543155952345 }, { 0,
// -0.887062599768095, 0.519096129206812 }, { 0, -0.887062599768095,
// 0.730152005574049 }, { 0, -0.887062599768095, 0.887062599768095 }, { 0,
// -0.887062599768095, 0.978228658146057 }, { 0, -0.730152005574049,
//-0.978228658146057 }, { 0, -0.730152005574049, -0.887062599768095 }, { 0,
//-0.730152005574049, -0.730152005574049 }, { 0, -0.730152005574049,
//-0.519096129206812 }, { 0, -0.730152005574049, -0.269543155952345 }, { 0,
//-0.730152005574049, 0 }, { 0, -0.730152005574049, 0.269543155952345 }, { 0,
// -0.730152005574049, 0.519096129206812 }, { 0, -0.730152005574049,
// 0.730152005574049 }, { 0, -0.730152005574049, 0.887062599768095 }, { 0,
// -0.730152005574049, 0.978228658146057 }, { 0, -0.519096129206812,
//-0.978228658146057 }, { 0, -0.519096129206812, -0.887062599768095 }, { 0,
//-0.519096129206812, -0.730152005574049 }, { 0, -0.519096129206812,
//-0.519096129206812 }, { 0, -0.519096129206812, -0.269543155952345 }, { 0,
//-0.519096129206812, 0 }, { 0, -0.519096129206812, 0.269543155952345 }, { 0,
// -0.519096129206812, 0.519096129206812 }, { 0, -0.519096129206812,
// 0.730152005574049 }, { 0, -0.519096129206812, 0.887062599768095 }, { 0,
// -0.519096129206812, 0.978228658146057 }, { 0, -0.269543155952345,
//-0.978228658146057 }, { 0, -0.269543155952345, -0.887062599768095 }, { 0,
//-0.269543155952345, -0.730152005574049 }, { 0, -0.269543155952345,
//-0.519096129206812 }, { 0, -0.269543155952345, -0.269543155952345 }, { 0,
//-0.269543155952345, 0 }, { 0, -0.269543155952345, 0.269543155952345 }, { 0,
// -0.269543155952345, 0.519096129206812 }, { 0, -0.269543155952345,
// 0.730152005574049 }, { 0, -0.269543155952345, 0.887062599768095 }, { 0,
// -0.269543155952345, 0.978228658146057 }, { 0, 0, -0.978228658146057 }, { 0,
// 0, -0.887062599768095 }, { 0, 0, -0.730152005574049
//}, { 0, 0, -0.519096129206812 }, { 0, 0, -0.269543155952345 }, { 0, 0, 0 }, {
// 0, 0, 0.269543155952345 }, { 0, 0, 0.519096129206812 }, { 0, 0,
// 0.730152005574049 }, { 0, 0, 0.887062599768095 }, { 0, 0, 0.978228658146057
// }, { 0, 0.269543155952345, -0.978228658146057 }, { 0, 0.269543155952345,
//-0.887062599768095 }, { 0, 0.269543155952345, -0.730152005574049 }, { 0,
// 0.269543155952345, -0.519096129206812 }, { 0, 0.269543155952345,
//-0.269543155952345 }, { 0, 0.269543155952345, 0 }, { 0, 0.269543155952345,
// 0.269543155952345 }, { 0, 0.269543155952345, 0.519096129206812 }, { 0,
// 0.269543155952345, 0.730152005574049 }, { 0, 0.269543155952345,
// 0.887062599768095 }, { 0, 0.269543155952345, 0.978228658146057 }, { 0,
// 0.519096129206812, -0.978228658146057 }, { 0, 0.519096129206812,
//-0.887062599768095 }, { 0, 0.519096129206812, -0.730152005574049 }, { 0,
// 0.519096129206812, -0.519096129206812 }, { 0, 0.519096129206812,
//-0.269543155952345 }, { 0, 0.519096129206812, 0 }, { 0, 0.519096129206812,
// 0.269543155952345 }, { 0, 0.519096129206812, 0.519096129206812 }, { 0,
// 0.519096129206812, 0.730152005574049 }, { 0, 0.519096129206812,
// 0.887062599768095 }, { 0, 0.519096129206812, 0.978228658146057 }, { 0,
// 0.730152005574049, -0.978228658146057 }, { 0, 0.730152005574049,
//-0.887062599768095 }, { 0, 0.730152005574049, -0.730152005574049 }, { 0,
// 0.730152005574049, -0.519096129206812 }, { 0, 0.730152005574049,
//-0.269543155952345 }, { 0, 0.730152005574049, 0 }, { 0, 0.730152005574049,
// 0.269543155952345 }, { 0, 0.730152005574049, 0.519096129206812 }, { 0,
// 0.730152005574049, 0.730152005574049 }, { 0, 0.730152005574049,
// 0.887062599768095 }, { 0, 0.730152005574049, 0.978228658146057 }, { 0,
// 0.887062599768095, -0.978228658146057 }, { 0, 0.887062599768095,
//-0.887062599768095 }, { 0, 0.887062599768095, -0.730152005574049 }, { 0,
// 0.887062599768095, -0.519096129206812 }, { 0, 0.887062599768095,
//-0.269543155952345 }, { 0, 0.887062599768095, 0 }, { 0, 0.887062599768095,
// 0.269543155952345 }, { 0, 0.887062599768095, 0.519096129206812 }, { 0,
// 0.887062599768095, 0.730152005574049 }, { 0, 0.887062599768095,
// 0.887062599768095 }, { 0, 0.887062599768095, 0.978228658146057 }, { 0,
// 0.978228658146057, -0.978228658146057 }, { 0, 0.978228658146057,
//-0.887062599768095 }, { 0, 0.978228658146057, -0.730152005574049 }, { 0,
// 0.978228658146057, -0.519096129206812 }, { 0, 0.978228658146057,
//-0.269543155952345 }, { 0, 0.978228658146057, 0 }, { 0, 0.978228658146057,
// 0.269543155952345 }, { 0, 0.978228658146057, 0.519096129206812 }, { 0,
// 0.978228658146057, 0.730152005574049 }, { 0, 0.978228658146057,
// 0.887062599768095 }, { 0, 0.978228658146057, 0.978228658146057 }, {
// 0.269543155952345, -0.978228658146057, -0.978228658146057 }, {
// 0.269543155952345, -0.978228658146057, -0.887062599768095 }, {
// 0.269543155952345, -0.978228658146057, -0.730152005574049 }, {
// 0.269543155952345, -0.978228658146057, -0.519096129206812 }, {
// 0.269543155952345, -0.978228658146057, -0.269543155952345 }, {
// 0.269543155952345, -0.978228658146057, 0 }, { 0.269543155952345,
//-0.978228658146057, 0.269543155952345 }, { 0.269543155952345,
//-0.978228658146057, 0.519096129206812 }, { 0.269543155952345,
//-0.978228658146057, 0.730152005574049 }, { 0.269543155952345,
//-0.978228658146057, 0.887062599768095 }, { 0.269543155952345,
//-0.978228658146057, 0.978228658146057 }, { 0.269543155952345,
//-0.887062599768095, -0.978228658146057 }, { 0.269543155952345,
//-0.887062599768095, -0.887062599768095 }, { 0.269543155952345,
//-0.887062599768095, -0.730152005574049 }, { 0.269543155952345,
//-0.887062599768095, -0.519096129206812 }, { 0.269543155952345,
//-0.887062599768095, -0.269543155952345 }, { 0.269543155952345,
//-0.887062599768095, 0 }, { 0.269543155952345, -0.887062599768095,
// 0.269543155952345 }, { 0.269543155952345, -0.887062599768095,
// 0.519096129206812 }, { 0.269543155952345, -0.887062599768095,
// 0.730152005574049 }, { 0.269543155952345, -0.887062599768095,
// 0.887062599768095 }, { 0.269543155952345, -0.887062599768095,
// 0.978228658146057 }, { 0.269543155952345, -0.730152005574049,
//-0.978228658146057 }, { 0.269543155952345, -0.730152005574049,
//-0.887062599768095 }, { 0.269543155952345, -0.730152005574049,
//-0.730152005574049 }, { 0.269543155952345, -0.730152005574049,
//-0.519096129206812 }, { 0.269543155952345, -0.730152005574049,
//-0.269543155952345 }, { 0.269543155952345, -0.730152005574049, 0 }, {
// 0.269543155952345, -0.730152005574049, 0.269543155952345 }, {
// 0.269543155952345, -0.730152005574049, 0.519096129206812 }, {
// 0.269543155952345, -0.730152005574049, 0.730152005574049 }, {
// 0.269543155952345, -0.730152005574049, 0.887062599768095 }, {
// 0.269543155952345, -0.730152005574049, 0.978228658146057 }, {
// 0.269543155952345, -0.519096129206812, -0.978228658146057 }, {
// 0.269543155952345, -0.519096129206812, -0.887062599768095 }, {
// 0.269543155952345, -0.519096129206812, -0.730152005574049 }, {
// 0.269543155952345, -0.519096129206812, -0.519096129206812 }, {
// 0.269543155952345, -0.519096129206812, -0.269543155952345 }, {
// 0.269543155952345, -0.519096129206812, 0 }, { 0.269543155952345,
//-0.519096129206812, 0.269543155952345 }, { 0.269543155952345,
//-0.519096129206812, 0.519096129206812 }, { 0.269543155952345,
//-0.519096129206812, 0.730152005574049 }, { 0.269543155952345,
//-0.519096129206812, 0.887062599768095 }, { 0.269543155952345,
//-0.519096129206812, 0.978228658146057 }, { 0.269543155952345,
//-0.269543155952345, -0.978228658146057 }, { 0.269543155952345,
//-0.269543155952345, -0.887062599768095 }, { 0.269543155952345,
//-0.269543155952345, -0.730152005574049 }, { 0.269543155952345,
//-0.269543155952345, -0.519096129206812 }, { 0.269543155952345,
//-0.269543155952345, -0.269543155952345 }, { 0.269543155952345,
//-0.269543155952345, 0 }, { 0.269543155952345, -0.269543155952345,
// 0.269543155952345 }, { 0.269543155952345, -0.269543155952345,
// 0.519096129206812 }, { 0.269543155952345, -0.269543155952345,
// 0.730152005574049 }, { 0.269543155952345, -0.269543155952345,
// 0.887062599768095 }, { 0.269543155952345, -0.269543155952345,
// 0.978228658146057 }, { 0.269543155952345, 0, -0.978228658146057 }, {
// 0.269543155952345, 0, -0.887062599768095 }, { 0.269543155952345, 0,
//-0.730152005574049 }, { 0.269543155952345, 0, -0.519096129206812 }, {
// 0.269543155952345, 0, -0.269543155952345 }, { 0.269543155952345, 0, 0 }, {
// 0.269543155952345, 0, 0.269543155952345 }, { 0.269543155952345, 0,
// 0.519096129206812 }, { 0.269543155952345, 0, 0.730152005574049 }, {
// 0.269543155952345, 0, 0.887062599768095 }, { 0.269543155952345, 0,
// 0.978228658146057 }, { 0.269543155952345, 0.269543155952345,
//-0.978228658146057 }, { 0.269543155952345, 0.269543155952345,
//-0.887062599768095 }, { 0.269543155952345, 0.269543155952345,
//-0.730152005574049 }, { 0.269543155952345, 0.269543155952345,
//-0.519096129206812 }, { 0.269543155952345, 0.269543155952345,
//-0.269543155952345 }, { 0.269543155952345, 0.269543155952345, 0 }, {
// 0.269543155952345, 0.269543155952345, 0.269543155952345 }, {
// 0.269543155952345, 0.269543155952345, 0.519096129206812 }, {
// 0.269543155952345, 0.269543155952345, 0.730152005574049 }, {
// 0.269543155952345, 0.269543155952345, 0.887062599768095 }, {
// 0.269543155952345, 0.269543155952345, 0.978228658146057 }, {
// 0.269543155952345, 0.519096129206812, -0.978228658146057 }, {
// 0.269543155952345, 0.519096129206812, -0.887062599768095 }, {
// 0.269543155952345, 0.519096129206812, -0.730152005574049 }, {
// 0.269543155952345, 0.519096129206812, -0.519096129206812 }, {
// 0.269543155952345, 0.519096129206812, -0.269543155952345 }, {
// 0.269543155952345, 0.519096129206812, 0 }, { 0.269543155952345,
// 0.519096129206812, 0.269543155952345 }, { 0.269543155952345,
// 0.519096129206812, 0.519096129206812 }, { 0.269543155952345,
// 0.519096129206812, 0.730152005574049 }, { 0.269543155952345,
// 0.519096129206812, 0.887062599768095 }, { 0.269543155952345,
// 0.519096129206812, 0.978228658146057 }, { 0.269543155952345,
// 0.730152005574049, -0.978228658146057 }, { 0.269543155952345,
// 0.730152005574049, -0.887062599768095 }, { 0.269543155952345,
// 0.730152005574049, -0.730152005574049 }, { 0.269543155952345,
// 0.730152005574049, -0.519096129206812 }, { 0.269543155952345,
// 0.730152005574049, -0.269543155952345 }, { 0.269543155952345,
// 0.730152005574049, 0 }, { 0.269543155952345, 0.730152005574049,
// 0.269543155952345 }, { 0.269543155952345, 0.730152005574049,
// 0.519096129206812
//}, { 0.269543155952345, 0.730152005574049, 0.730152005574049 }, {
// 0.269543155952345, 0.730152005574049, 0.887062599768095 }, {
// 0.269543155952345, 0.730152005574049, 0.978228658146057 }, {
// 0.269543155952345, 0.887062599768095, -0.978228658146057 }, {
// 0.269543155952345, 0.887062599768095, -0.887062599768095 }, {
// 0.269543155952345, 0.887062599768095, -0.730152005574049 }, {
// 0.269543155952345, 0.887062599768095, -0.519096129206812 }, {
// 0.269543155952345, 0.887062599768095, -0.269543155952345 }, {
// 0.269543155952345, 0.887062599768095, 0 }, { 0.269543155952345,
// 0.887062599768095, 0.269543155952345 }, { 0.269543155952345,
// 0.887062599768095, 0.519096129206812 }, { 0.269543155952345,
// 0.887062599768095, 0.730152005574049 }, { 0.269543155952345,
// 0.887062599768095, 0.887062599768095 }, { 0.269543155952345,
// 0.887062599768095, 0.978228658146057 }, { 0.269543155952345,
// 0.978228658146057, -0.978228658146057 }, { 0.269543155952345,
// 0.978228658146057, -0.887062599768095 }, { 0.269543155952345,
// 0.978228658146057, -0.730152005574049 }, { 0.269543155952345,
// 0.978228658146057, -0.519096129206812 }, { 0.269543155952345,
// 0.978228658146057, -0.269543155952345 }, { 0.269543155952345,
// 0.978228658146057, 0 }, { 0.269543155952345, 0.978228658146057,
// 0.269543155952345 }, { 0.269543155952345, 0.978228658146057,
// 0.519096129206812
//}, { 0.269543155952345, 0.978228658146057, 0.730152005574049 }, {
// 0.269543155952345, 0.978228658146057, 0.887062599768095 }, {
// 0.269543155952345, 0.978228658146057, 0.978228658146057 }, {
// 0.519096129206812, -0.978228658146057, -0.978228658146057 }, {
// 0.519096129206812, -0.978228658146057, -0.887062599768095 }, {
// 0.519096129206812, -0.978228658146057, -0.730152005574049 }, {
// 0.519096129206812, -0.978228658146057, -0.519096129206812 }, {
// 0.519096129206812, -0.978228658146057, -0.269543155952345 }, {
// 0.519096129206812, -0.978228658146057, 0 }, { 0.519096129206812,
//-0.978228658146057, 0.269543155952345 }, { 0.519096129206812,
//-0.978228658146057, 0.519096129206812 }, { 0.519096129206812,
//-0.978228658146057, 0.730152005574049 }, { 0.519096129206812,
//-0.978228658146057, 0.887062599768095 }, { 0.519096129206812,
//-0.978228658146057, 0.978228658146057 }, { 0.519096129206812,
//-0.887062599768095, -0.978228658146057 }, { 0.519096129206812,
//-0.887062599768095, -0.887062599768095 }, { 0.519096129206812,
//-0.887062599768095, -0.730152005574049 }, { 0.519096129206812,
//-0.887062599768095, -0.519096129206812 }, { 0.519096129206812,
//-0.887062599768095, -0.269543155952345 }, { 0.519096129206812,
//-0.887062599768095, 0 }, { 0.519096129206812, -0.887062599768095,
// 0.269543155952345 }, { 0.519096129206812, -0.887062599768095,
// 0.519096129206812 }, { 0.519096129206812, -0.887062599768095,
// 0.730152005574049 }, { 0.519096129206812, -0.887062599768095,
// 0.887062599768095 }, { 0.519096129206812, -0.887062599768095,
// 0.978228658146057 }, { 0.519096129206812, -0.730152005574049,
//-0.978228658146057 }, { 0.519096129206812, -0.730152005574049,
//-0.887062599768095 }, { 0.519096129206812, -0.730152005574049,
//-0.730152005574049 }, { 0.519096129206812, -0.730152005574049,
//-0.519096129206812 }, { 0.519096129206812, -0.730152005574049,
//-0.269543155952345 }, { 0.519096129206812, -0.730152005574049, 0 }, {
// 0.519096129206812, -0.730152005574049, 0.269543155952345 }, {
// 0.519096129206812, -0.730152005574049, 0.519096129206812 }, {
// 0.519096129206812, -0.730152005574049, 0.730152005574049 }, {
// 0.519096129206812, -0.730152005574049, 0.887062599768095 }, {
// 0.519096129206812, -0.730152005574049, 0.978228658146057 }, {
// 0.519096129206812, -0.519096129206812, -0.978228658146057 }, {
// 0.519096129206812, -0.519096129206812, -0.887062599768095 }, {
// 0.519096129206812, -0.519096129206812, -0.730152005574049 }, {
// 0.519096129206812, -0.519096129206812, -0.519096129206812 }, {
// 0.519096129206812, -0.519096129206812, -0.269543155952345 }, {
// 0.519096129206812, -0.519096129206812, 0 }, { 0.519096129206812,
//-0.519096129206812, 0.269543155952345 }, { 0.519096129206812,
//-0.519096129206812, 0.519096129206812 }, { 0.519096129206812,
//-0.519096129206812, 0.730152005574049 }, { 0.519096129206812,
//-0.519096129206812, 0.887062599768095 }, { 0.519096129206812,
//-0.519096129206812, 0.978228658146057 }, { 0.519096129206812,
//-0.269543155952345, -0.978228658146057 }, { 0.519096129206812,
//-0.269543155952345, -0.887062599768095 }, { 0.519096129206812,
//-0.269543155952345, -0.730152005574049 }, { 0.519096129206812,
//-0.269543155952345, -0.519096129206812 }, { 0.519096129206812,
//-0.269543155952345, -0.269543155952345 }, { 0.519096129206812,
//-0.269543155952345, 0 }, { 0.519096129206812, -0.269543155952345,
// 0.269543155952345 }, { 0.519096129206812, -0.269543155952345,
// 0.519096129206812 }, { 0.519096129206812, -0.269543155952345,
// 0.730152005574049 }, { 0.519096129206812, -0.269543155952345,
// 0.887062599768095 }, { 0.519096129206812, -0.269543155952345,
// 0.978228658146057 }, { 0.519096129206812, 0, -0.978228658146057 }, {
// 0.519096129206812, 0, -0.887062599768095 }, { 0.519096129206812, 0,
//-0.730152005574049 }, { 0.519096129206812, 0, -0.519096129206812 }, {
// 0.519096129206812, 0, -0.269543155952345 }, { 0.519096129206812, 0, 0 }, {
// 0.519096129206812, 0, 0.269543155952345 }, { 0.519096129206812, 0,
// 0.519096129206812 }, { 0.519096129206812, 0, 0.730152005574049 }, {
// 0.519096129206812, 0, 0.887062599768095 }, { 0.519096129206812, 0,
// 0.978228658146057 }, { 0.519096129206812, 0.269543155952345,
//-0.978228658146057 }, { 0.519096129206812, 0.269543155952345,
//-0.887062599768095 }, { 0.519096129206812, 0.269543155952345,
//-0.730152005574049 }, { 0.519096129206812, 0.269543155952345,
//-0.519096129206812 }, { 0.519096129206812, 0.269543155952345,
//-0.269543155952345 }, { 0.519096129206812, 0.269543155952345, 0 }, {
// 0.519096129206812, 0.269543155952345, 0.269543155952345 }, {
// 0.519096129206812, 0.269543155952345, 0.519096129206812 }, {
// 0.519096129206812, 0.269543155952345, 0.730152005574049 }, {
// 0.519096129206812, 0.269543155952345, 0.887062599768095 }, {
// 0.519096129206812, 0.269543155952345, 0.978228658146057 }, {
// 0.519096129206812, 0.519096129206812, -0.978228658146057 }, {
// 0.519096129206812, 0.519096129206812, -0.887062599768095 }, {
// 0.519096129206812, 0.519096129206812, -0.730152005574049 }, {
// 0.519096129206812, 0.519096129206812, -0.519096129206812 }, {
// 0.519096129206812, 0.519096129206812, -0.269543155952345 }, {
// 0.519096129206812, 0.519096129206812, 0 }, { 0.519096129206812,
// 0.519096129206812, 0.269543155952345 }, { 0.519096129206812,
// 0.519096129206812, 0.519096129206812 }, { 0.519096129206812,
// 0.519096129206812, 0.730152005574049 }, { 0.519096129206812,
// 0.519096129206812, 0.887062599768095 }, { 0.519096129206812,
// 0.519096129206812, 0.978228658146057 }, { 0.519096129206812,
// 0.730152005574049, -0.978228658146057 }, { 0.519096129206812,
// 0.730152005574049, -0.887062599768095 }, { 0.519096129206812,
// 0.730152005574049, -0.730152005574049 }, { 0.519096129206812,
// 0.730152005574049, -0.519096129206812 }, { 0.519096129206812,
// 0.730152005574049, -0.269543155952345 }, { 0.519096129206812,
// 0.730152005574049, 0 }, { 0.519096129206812, 0.730152005574049,
// 0.269543155952345 }, { 0.519096129206812, 0.730152005574049,
// 0.519096129206812
//}, { 0.519096129206812, 0.730152005574049, 0.730152005574049 }, {
// 0.519096129206812, 0.730152005574049, 0.887062599768095 }, {
// 0.519096129206812, 0.730152005574049, 0.978228658146057 }, {
// 0.519096129206812, 0.887062599768095, -0.978228658146057 }, {
// 0.519096129206812, 0.887062599768095, -0.887062599768095 }, {
// 0.519096129206812, 0.887062599768095, -0.730152005574049 }, {
// 0.519096129206812, 0.887062599768095, -0.519096129206812 }, {
// 0.519096129206812, 0.887062599768095, -0.269543155952345 }, {
// 0.519096129206812, 0.887062599768095, 0 }, { 0.519096129206812,
// 0.887062599768095, 0.269543155952345 }, { 0.519096129206812,
// 0.887062599768095, 0.519096129206812 }, { 0.519096129206812,
// 0.887062599768095, 0.730152005574049 }, { 0.519096129206812,
// 0.887062599768095, 0.887062599768095 }, { 0.519096129206812,
// 0.887062599768095, 0.978228658146057 }, { 0.519096129206812,
// 0.978228658146057, -0.978228658146057 }, { 0.519096129206812,
// 0.978228658146057, -0.887062599768095 }, { 0.519096129206812,
// 0.978228658146057, -0.730152005574049 }, { 0.519096129206812,
// 0.978228658146057, -0.519096129206812 }, { 0.519096129206812,
// 0.978228658146057, -0.269543155952345 }, { 0.519096129206812,
// 0.978228658146057, 0 }, { 0.519096129206812, 0.978228658146057,
// 0.269543155952345 }, { 0.519096129206812, 0.978228658146057,
// 0.519096129206812
//}, { 0.519096129206812, 0.978228658146057, 0.730152005574049 }, {
// 0.519096129206812, 0.978228658146057, 0.887062599768095 }, {
// 0.519096129206812, 0.978228658146057, 0.978228658146057 }, {
// 0.730152005574049, -0.978228658146057, -0.978228658146057 }, {
// 0.730152005574049, -0.978228658146057, -0.887062599768095 }, {
// 0.730152005574049, -0.978228658146057, -0.730152005574049 }, {
// 0.730152005574049, -0.978228658146057, -0.519096129206812 }, {
// 0.730152005574049, -0.978228658146057, -0.269543155952345 }, {
// 0.730152005574049, -0.978228658146057, 0 }, { 0.730152005574049,
//-0.978228658146057, 0.269543155952345 }, { 0.730152005574049,
//-0.978228658146057, 0.519096129206812 }, { 0.730152005574049,
//-0.978228658146057, 0.730152005574049 }, { 0.730152005574049,
//-0.978228658146057, 0.887062599768095 }, { 0.730152005574049,
//-0.978228658146057, 0.978228658146057 }, { 0.730152005574049,
//-0.887062599768095, -0.978228658146057 }, { 0.730152005574049,
//-0.887062599768095, -0.887062599768095 }, { 0.730152005574049,
//-0.887062599768095, -0.730152005574049 }, { 0.730152005574049,
//-0.887062599768095, -0.519096129206812 }, { 0.730152005574049,
//-0.887062599768095, -0.269543155952345 }, { 0.730152005574049,
//-0.887062599768095, 0 }, { 0.730152005574049, -0.887062599768095,
// 0.269543155952345 }, { 0.730152005574049, -0.887062599768095,
// 0.519096129206812 }, { 0.730152005574049, -0.887062599768095,
// 0.730152005574049 }, { 0.730152005574049, -0.887062599768095,
// 0.887062599768095 }, { 0.730152005574049, -0.887062599768095,
// 0.978228658146057 }, { 0.730152005574049, -0.730152005574049,
//-0.978228658146057 }, { 0.730152005574049, -0.730152005574049,
//-0.887062599768095 }, { 0.730152005574049, -0.730152005574049,
//-0.730152005574049 }, { 0.730152005574049, -0.730152005574049,
//-0.519096129206812 }, { 0.730152005574049, -0.730152005574049,
//-0.269543155952345 }, { 0.730152005574049, -0.730152005574049, 0 }, {
// 0.730152005574049, -0.730152005574049, 0.269543155952345 }, {
// 0.730152005574049, -0.730152005574049, 0.519096129206812 }, {
// 0.730152005574049, -0.730152005574049, 0.730152005574049 }, {
// 0.730152005574049, -0.730152005574049, 0.887062599768095 }, {
// 0.730152005574049, -0.730152005574049, 0.978228658146057 }, {
// 0.730152005574049, -0.519096129206812, -0.978228658146057 }, {
// 0.730152005574049, -0.519096129206812, -0.887062599768095 }, {
// 0.730152005574049, -0.519096129206812, -0.730152005574049 }, {
// 0.730152005574049, -0.519096129206812, -0.519096129206812 }, {
// 0.730152005574049, -0.519096129206812, -0.269543155952345 }, {
// 0.730152005574049, -0.519096129206812, 0 }, { 0.730152005574049,
//-0.519096129206812, 0.269543155952345 }, { 0.730152005574049,
//-0.519096129206812, 0.519096129206812 }, { 0.730152005574049,
//-0.519096129206812, 0.730152005574049 }, { 0.730152005574049,
//-0.519096129206812, 0.887062599768095 }, { 0.730152005574049,
//-0.519096129206812, 0.978228658146057 }, { 0.730152005574049,
//-0.269543155952345, -0.978228658146057 }, { 0.730152005574049,
//-0.269543155952345, -0.887062599768095 }, { 0.730152005574049,
//-0.269543155952345, -0.730152005574049 }, { 0.730152005574049,
//-0.269543155952345, -0.519096129206812 }, { 0.730152005574049,
//-0.269543155952345, -0.269543155952345 }, { 0.730152005574049,
//-0.269543155952345, 0 }, { 0.730152005574049, -0.269543155952345,
// 0.269543155952345 }, { 0.730152005574049, -0.269543155952345,
// 0.519096129206812 }, { 0.730152005574049, -0.269543155952345,
// 0.730152005574049 }, { 0.730152005574049, -0.269543155952345,
// 0.887062599768095 }, { 0.730152005574049, -0.269543155952345,
// 0.978228658146057 }, { 0.730152005574049, 0, -0.978228658146057 }, {
// 0.730152005574049, 0, -0.887062599768095 }, { 0.730152005574049, 0,
//-0.730152005574049 }, { 0.730152005574049, 0, -0.519096129206812 }, {
// 0.730152005574049, 0, -0.269543155952345 }, { 0.730152005574049, 0, 0 }, {
// 0.730152005574049, 0, 0.269543155952345 }, { 0.730152005574049, 0,
// 0.519096129206812 }, { 0.730152005574049, 0, 0.730152005574049 }, {
// 0.730152005574049, 0, 0.887062599768095 }, { 0.730152005574049, 0,
// 0.978228658146057 }, { 0.730152005574049, 0.269543155952345,
//-0.978228658146057 }, { 0.730152005574049, 0.269543155952345,
//-0.887062599768095 }, { 0.730152005574049, 0.269543155952345,
//-0.730152005574049 }, { 0.730152005574049, 0.269543155952345,
//-0.519096129206812 }, { 0.730152005574049, 0.269543155952345,
//-0.269543155952345 }, { 0.730152005574049, 0.269543155952345, 0 }, {
// 0.730152005574049, 0.269543155952345, 0.269543155952345 }, {
// 0.730152005574049, 0.269543155952345, 0.519096129206812 }, {
// 0.730152005574049, 0.269543155952345, 0.730152005574049 }, {
// 0.730152005574049, 0.269543155952345, 0.887062599768095 }, {
// 0.730152005574049, 0.269543155952345, 0.978228658146057 }, {
// 0.730152005574049, 0.519096129206812, -0.978228658146057 }, {
// 0.730152005574049, 0.519096129206812, -0.887062599768095 }, {
// 0.730152005574049, 0.519096129206812, -0.730152005574049 }, {
// 0.730152005574049, 0.519096129206812, -0.519096129206812 }, {
// 0.730152005574049, 0.519096129206812, -0.269543155952345 }, {
// 0.730152005574049, 0.519096129206812, 0 }, { 0.730152005574049,
// 0.519096129206812, 0.269543155952345 }, { 0.730152005574049,
// 0.519096129206812, 0.519096129206812 }, { 0.730152005574049,
// 0.519096129206812, 0.730152005574049 }, { 0.730152005574049,
// 0.519096129206812, 0.887062599768095 }, { 0.730152005574049,
// 0.519096129206812, 0.978228658146057 }, { 0.730152005574049,
// 0.730152005574049, -0.978228658146057 }, { 0.730152005574049,
// 0.730152005574049, -0.887062599768095 }, { 0.730152005574049,
// 0.730152005574049, -0.730152005574049 }, { 0.730152005574049,
// 0.730152005574049, -0.519096129206812 }, { 0.730152005574049,
// 0.730152005574049, -0.269543155952345 }, { 0.730152005574049,
// 0.730152005574049, 0 }, { 0.730152005574049, 0.730152005574049,
// 0.269543155952345 }, { 0.730152005574049, 0.730152005574049,
// 0.519096129206812
//}, { 0.730152005574049, 0.730152005574049, 0.730152005574049 }, {
// 0.730152005574049, 0.730152005574049, 0.887062599768095 }, {
// 0.730152005574049, 0.730152005574049, 0.978228658146057 }, {
// 0.730152005574049, 0.887062599768095, -0.978228658146057 }, {
// 0.730152005574049, 0.887062599768095, -0.887062599768095 }, {
// 0.730152005574049, 0.887062599768095, -0.730152005574049 }, {
// 0.730152005574049, 0.887062599768095, -0.519096129206812 }, {
// 0.730152005574049, 0.887062599768095, -0.269543155952345 }, {
// 0.730152005574049, 0.887062599768095, 0 }, { 0.730152005574049,
// 0.887062599768095, 0.269543155952345 }, { 0.730152005574049,
// 0.887062599768095, 0.519096129206812 }, { 0.730152005574049,
// 0.887062599768095, 0.730152005574049 }, { 0.730152005574049,
// 0.887062599768095, 0.887062599768095 }, { 0.730152005574049,
// 0.887062599768095, 0.978228658146057 }, { 0.730152005574049,
// 0.978228658146057, -0.978228658146057 }, { 0.730152005574049,
// 0.978228658146057, -0.887062599768095 }, { 0.730152005574049,
// 0.978228658146057, -0.730152005574049 }, { 0.730152005574049,
// 0.978228658146057, -0.519096129206812 }, { 0.730152005574049,
// 0.978228658146057, -0.269543155952345 }, { 0.730152005574049,
// 0.978228658146057, 0 }, { 0.730152005574049, 0.978228658146057,
// 0.269543155952345 }, { 0.730152005574049, 0.978228658146057,
// 0.519096129206812
//}, { 0.730152005574049, 0.978228658146057, 0.730152005574049 }, {
// 0.730152005574049, 0.978228658146057, 0.887062599768095 }, {
// 0.730152005574049, 0.978228658146057, 0.978228658146057 }, {
// 0.887062599768095, -0.978228658146057, -0.978228658146057 }, {
// 0.887062599768095, -0.978228658146057, -0.887062599768095 }, {
// 0.887062599768095, -0.978228658146057, -0.730152005574049 }, {
// 0.887062599768095, -0.978228658146057, -0.519096129206812 }, {
// 0.887062599768095, -0.978228658146057, -0.269543155952345 }, {
// 0.887062599768095, -0.978228658146057, 0 }, { 0.887062599768095,
//-0.978228658146057, 0.269543155952345 }, { 0.887062599768095,
//-0.978228658146057, 0.519096129206812 }, { 0.887062599768095,
//-0.978228658146057, 0.730152005574049 }, { 0.887062599768095,
//-0.978228658146057, 0.887062599768095 }, { 0.887062599768095,
//-0.978228658146057, 0.978228658146057 }, { 0.887062599768095,
//-0.887062599768095, -0.978228658146057 }, { 0.887062599768095,
//-0.887062599768095, -0.887062599768095 }, { 0.887062599768095,
//-0.887062599768095, -0.730152005574049 }, { 0.887062599768095,
//-0.887062599768095, -0.519096129206812 }, { 0.887062599768095,
//-0.887062599768095, -0.269543155952345 }, { 0.887062599768095,
//-0.887062599768095, 0 }, { 0.887062599768095, -0.887062599768095,
// 0.269543155952345 }, { 0.887062599768095, -0.887062599768095,
// 0.519096129206812 }, { 0.887062599768095, -0.887062599768095,
// 0.730152005574049 }, { 0.887062599768095, -0.887062599768095,
// 0.887062599768095 }, { 0.887062599768095, -0.887062599768095,
// 0.978228658146057 }, { 0.887062599768095, -0.730152005574049,
//-0.978228658146057 }, { 0.887062599768095, -0.730152005574049,
//-0.887062599768095 }, { 0.887062599768095, -0.730152005574049,
//-0.730152005574049 }, { 0.887062599768095, -0.730152005574049,
//-0.519096129206812 }, { 0.887062599768095, -0.730152005574049,
//-0.269543155952345 }, { 0.887062599768095, -0.730152005574049, 0 }, {
// 0.887062599768095, -0.730152005574049, 0.269543155952345 }, {
// 0.887062599768095, -0.730152005574049, 0.519096129206812 }, {
// 0.887062599768095, -0.730152005574049, 0.730152005574049 }, {
// 0.887062599768095, -0.730152005574049, 0.887062599768095 }, {
// 0.887062599768095, -0.730152005574049, 0.978228658146057 }, {
// 0.887062599768095, -0.519096129206812, -0.978228658146057 }, {
// 0.887062599768095, -0.519096129206812, -0.887062599768095 }, {
// 0.887062599768095, -0.519096129206812, -0.730152005574049 }, {
// 0.887062599768095, -0.519096129206812, -0.519096129206812 }, {
// 0.887062599768095, -0.519096129206812, -0.269543155952345 }, {
// 0.887062599768095, -0.519096129206812, 0 }, { 0.887062599768095,
//-0.519096129206812, 0.269543155952345 }, { 0.887062599768095,
//-0.519096129206812, 0.519096129206812 }, { 0.887062599768095,
//-0.519096129206812, 0.730152005574049 }, { 0.887062599768095,
//-0.519096129206812, 0.887062599768095 }, { 0.887062599768095,
//-0.519096129206812, 0.978228658146057 }, { 0.887062599768095,
//-0.269543155952345, -0.978228658146057 }, { 0.887062599768095,
//-0.269543155952345, -0.887062599768095 }, { 0.887062599768095,
//-0.269543155952345, -0.730152005574049 }, { 0.887062599768095,
//-0.269543155952345, -0.519096129206812 }, { 0.887062599768095,
//-0.269543155952345, -0.269543155952345 }, { 0.887062599768095,
//-0.269543155952345, 0 }, { 0.887062599768095, -0.269543155952345,
// 0.269543155952345 }, { 0.887062599768095, -0.269543155952345,
// 0.519096129206812 }, { 0.887062599768095, -0.269543155952345,
// 0.730152005574049 }, { 0.887062599768095, -0.269543155952345,
// 0.887062599768095 }, { 0.887062599768095, -0.269543155952345,
// 0.978228658146057 }, { 0.887062599768095, 0, -0.978228658146057 }, {
// 0.887062599768095, 0, -0.887062599768095 }, { 0.887062599768095, 0,
//-0.730152005574049 }, { 0.887062599768095, 0, -0.519096129206812 }, {
// 0.887062599768095, 0, -0.269543155952345 }, { 0.887062599768095, 0, 0 }, {
// 0.887062599768095, 0, 0.269543155952345 }, { 0.887062599768095, 0,
// 0.519096129206812 }, { 0.887062599768095, 0, 0.730152005574049 }, {
// 0.887062599768095, 0, 0.887062599768095 }, { 0.887062599768095, 0,
// 0.978228658146057 }, { 0.887062599768095, 0.269543155952345,
//-0.978228658146057 }, { 0.887062599768095, 0.269543155952345,
//-0.887062599768095 }, { 0.887062599768095, 0.269543155952345,
//-0.730152005574049 }, { 0.887062599768095, 0.269543155952345,
//-0.519096129206812 }, { 0.887062599768095, 0.269543155952345,
//-0.269543155952345 }, { 0.887062599768095, 0.269543155952345, 0 }, {
// 0.887062599768095, 0.269543155952345, 0.269543155952345 }, {
// 0.887062599768095, 0.269543155952345, 0.519096129206812 }, {
// 0.887062599768095, 0.269543155952345, 0.730152005574049 }, {
// 0.887062599768095, 0.269543155952345, 0.887062599768095 }, {
// 0.887062599768095, 0.269543155952345, 0.978228658146057 }, {
// 0.887062599768095, 0.519096129206812, -0.978228658146057 }, {
// 0.887062599768095, 0.519096129206812, -0.887062599768095 }, {
// 0.887062599768095, 0.519096129206812, -0.730152005574049 }, {
// 0.887062599768095, 0.519096129206812, -0.519096129206812 }, {
// 0.887062599768095, 0.519096129206812, -0.269543155952345 }, {
// 0.887062599768095, 0.519096129206812, 0 }, { 0.887062599768095,
// 0.519096129206812, 0.269543155952345 }, { 0.887062599768095,
// 0.519096129206812, 0.519096129206812 }, { 0.887062599768095,
// 0.519096129206812, 0.730152005574049 }, { 0.887062599768095,
// 0.519096129206812, 0.887062599768095 }, { 0.887062599768095,
// 0.519096129206812, 0.978228658146057 }, { 0.887062599768095,
// 0.730152005574049, -0.978228658146057 }, { 0.887062599768095,
// 0.730152005574049, -0.887062599768095 }, { 0.887062599768095,
// 0.730152005574049, -0.730152005574049 }, { 0.887062599768095,
// 0.730152005574049, -0.519096129206812 }, { 0.887062599768095,
// 0.730152005574049, -0.269543155952345 }, { 0.887062599768095,
// 0.730152005574049, 0 }, { 0.887062599768095, 0.730152005574049,
// 0.269543155952345 }, { 0.887062599768095, 0.730152005574049,
// 0.519096129206812
//}, { 0.887062599768095, 0.730152005574049, 0.730152005574049 }, {
// 0.887062599768095, 0.730152005574049, 0.887062599768095 }, {
// 0.887062599768095, 0.730152005574049, 0.978228658146057 }, {
// 0.887062599768095, 0.887062599768095, -0.978228658146057 }, {
// 0.887062599768095, 0.887062599768095, -0.887062599768095 }, {
// 0.887062599768095, 0.887062599768095, -0.730152005574049 }, {
// 0.887062599768095, 0.887062599768095, -0.519096129206812 }, {
// 0.887062599768095, 0.887062599768095, -0.269543155952345 }, {
// 0.887062599768095, 0.887062599768095, 0 }, { 0.887062599768095,
// 0.887062599768095, 0.269543155952345 }, { 0.887062599768095,
// 0.887062599768095, 0.519096129206812 }, { 0.887062599768095,
// 0.887062599768095, 0.730152005574049 }, { 0.887062599768095,
// 0.887062599768095, 0.887062599768095 }, { 0.887062599768095,
// 0.887062599768095, 0.978228658146057 }, { 0.887062599768095,
// 0.978228658146057, -0.978228658146057 }, { 0.887062599768095,
// 0.978228658146057, -0.887062599768095 }, { 0.887062599768095,
// 0.978228658146057, -0.730152005574049 }, { 0.887062599768095,
// 0.978228658146057, -0.519096129206812 }, { 0.887062599768095,
// 0.978228658146057, -0.269543155952345 }, { 0.887062599768095,
// 0.978228658146057, 0 }, { 0.887062599768095, 0.978228658146057,
// 0.269543155952345 }, { 0.887062599768095, 0.978228658146057,
// 0.519096129206812
//}, { 0.887062599768095, 0.978228658146057, 0.730152005574049 }, {
// 0.887062599768095, 0.978228658146057, 0.887062599768095 }, {
// 0.887062599768095, 0.978228658146057, 0.978228658146057 }, {
// 0.978228658146057, -0.978228658146057, -0.978228658146057 }, {
// 0.978228658146057, -0.978228658146057, -0.887062599768095 }, {
// 0.978228658146057, -0.978228658146057, -0.730152005574049 }, {
// 0.978228658146057, -0.978228658146057, -0.519096129206812 }, {
// 0.978228658146057, -0.978228658146057, -0.269543155952345 }, {
// 0.978228658146057, -0.978228658146057, 0 }, { 0.978228658146057,
//-0.978228658146057, 0.269543155952345 }, { 0.978228658146057,
//-0.978228658146057, 0.519096129206812 }, { 0.978228658146057,
//-0.978228658146057, 0.730152005574049 }, { 0.978228658146057,
//-0.978228658146057, 0.887062599768095 }, { 0.978228658146057,
//-0.978228658146057, 0.978228658146057 }, { 0.978228658146057,
//-0.887062599768095, -0.978228658146057 }, { 0.978228658146057,
//-0.887062599768095, -0.887062599768095 }, { 0.978228658146057,
//-0.887062599768095, -0.730152005574049 }, { 0.978228658146057,
//-0.887062599768095, -0.519096129206812 }, { 0.978228658146057,
//-0.887062599768095, -0.269543155952345 }, { 0.978228658146057,
//-0.887062599768095, 0 }, { 0.978228658146057, -0.887062599768095,
// 0.269543155952345 }, { 0.978228658146057, -0.887062599768095,
// 0.519096129206812 }, { 0.978228658146057, -0.887062599768095,
// 0.730152005574049 }, { 0.978228658146057, -0.887062599768095,
// 0.887062599768095 }, { 0.978228658146057, -0.887062599768095,
// 0.978228658146057 }, { 0.978228658146057, -0.730152005574049,
//-0.978228658146057 }, { 0.978228658146057, -0.730152005574049,
//-0.887062599768095 }, { 0.978228658146057, -0.730152005574049,
//-0.730152005574049 }, { 0.978228658146057, -0.730152005574049,
//-0.519096129206812 }, { 0.978228658146057, -0.730152005574049,
//-0.269543155952345 }, { 0.978228658146057, -0.730152005574049, 0 }, {
// 0.978228658146057, -0.730152005574049, 0.269543155952345 }, {
// 0.978228658146057, -0.730152005574049, 0.519096129206812 }, {
// 0.978228658146057, -0.730152005574049, 0.730152005574049 }, {
// 0.978228658146057, -0.730152005574049, 0.887062599768095 }, {
// 0.978228658146057, -0.730152005574049, 0.978228658146057 }, {
// 0.978228658146057, -0.519096129206812, -0.978228658146057 }, {
// 0.978228658146057, -0.519096129206812, -0.887062599768095 }, {
// 0.978228658146057, -0.519096129206812, -0.730152005574049 }, {
// 0.978228658146057, -0.519096129206812, -0.519096129206812 }, {
// 0.978228658146057, -0.519096129206812, -0.269543155952345 }, {
// 0.978228658146057, -0.519096129206812, 0 }, { 0.978228658146057,
//-0.519096129206812, 0.269543155952345 }, { 0.978228658146057,
//-0.519096129206812, 0.519096129206812 }, { 0.978228658146057,
//-0.519096129206812, 0.730152005574049 }, { 0.978228658146057,
//-0.519096129206812, 0.887062599768095 }, { 0.978228658146057,
//-0.519096129206812, 0.978228658146057 }, { 0.978228658146057,
//-0.269543155952345, -0.978228658146057 }, { 0.978228658146057,
//-0.269543155952345, -0.887062599768095 }, { 0.978228658146057,
//-0.269543155952345, -0.730152005574049 }, { 0.978228658146057,
//-0.269543155952345, -0.519096129206812 }, { 0.978228658146057,
//-0.269543155952345, -0.269543155952345 }, { 0.978228658146057,
//-0.269543155952345, 0 }, { 0.978228658146057, -0.269543155952345,
// 0.269543155952345 }, { 0.978228658146057, -0.269543155952345,
// 0.519096129206812 }, { 0.978228658146057, -0.269543155952345,
// 0.730152005574049 }, { 0.978228658146057, -0.269543155952345,
// 0.887062599768095 }, { 0.978228658146057, -0.269543155952345,
// 0.978228658146057 }, { 0.978228658146057, 0, -0.978228658146057 }, {
// 0.978228658146057, 0, -0.887062599768095 }, { 0.978228658146057, 0,
//-0.730152005574049 }, { 0.978228658146057, 0, -0.519096129206812 }, {
// 0.978228658146057, 0, -0.269543155952345 }, { 0.978228658146057, 0, 0 }, {
// 0.978228658146057, 0, 0.269543155952345 }, { 0.978228658146057, 0,
// 0.519096129206812 }, { 0.978228658146057, 0, 0.730152005574049 }, {
// 0.978228658146057, 0, 0.887062599768095 }, { 0.978228658146057, 0,
// 0.978228658146057 }, { 0.978228658146057, 0.269543155952345,
//-0.978228658146057 }, { 0.978228658146057, 0.269543155952345,
//-0.887062599768095 }, { 0.978228658146057, 0.269543155952345,
//-0.730152005574049 }, { 0.978228658146057, 0.269543155952345,
//-0.519096129206812 }, { 0.978228658146057, 0.269543155952345,
//-0.269543155952345 }, { 0.978228658146057, 0.269543155952345, 0 }, {
// 0.978228658146057, 0.269543155952345, 0.269543155952345 }, {
// 0.978228658146057, 0.269543155952345, 0.519096129206812 }, {
// 0.978228658146057, 0.269543155952345, 0.730152005574049 }, {
// 0.978228658146057, 0.269543155952345, 0.887062599768095 }, {
// 0.978228658146057, 0.269543155952345, 0.978228658146057 }, {
// 0.978228658146057, 0.519096129206812, -0.978228658146057 }, {
// 0.978228658146057, 0.519096129206812, -0.887062599768095 }, {
// 0.978228658146057, 0.519096129206812, -0.730152005574049 }, {
// 0.978228658146057, 0.519096129206812, -0.519096129206812 }, {
// 0.978228658146057, 0.519096129206812, -0.269543155952345 }, {
// 0.978228658146057, 0.519096129206812, 0 }, { 0.978228658146057,
// 0.519096129206812, 0.269543155952345 }, { 0.978228658146057,
// 0.519096129206812, 0.519096129206812 }, { 0.978228658146057,
// 0.519096129206812, 0.730152005574049 }, { 0.978228658146057,
// 0.519096129206812, 0.887062599768095 }, { 0.978228658146057,
// 0.519096129206812, 0.978228658146057 }, { 0.978228658146057,
// 0.730152005574049, -0.978228658146057 }, { 0.978228658146057,
// 0.730152005574049, -0.887062599768095 }, { 0.978228658146057,
// 0.730152005574049, -0.730152005574049 }, { 0.978228658146057,
// 0.730152005574049, -0.519096129206812 }, { 0.978228658146057,
// 0.730152005574049, -0.269543155952345 }, { 0.978228658146057,
// 0.730152005574049, 0 }, { 0.978228658146057, 0.730152005574049,
// 0.269543155952345 }, { 0.978228658146057, 0.730152005574049,
// 0.519096129206812
//}, { 0.978228658146057, 0.730152005574049, 0.730152005574049 }, {
// 0.978228658146057, 0.730152005574049, 0.887062599768095 }, {
// 0.978228658146057, 0.730152005574049, 0.978228658146057 }, {
// 0.978228658146057, 0.887062599768095, -0.978228658146057 }, {
// 0.978228658146057, 0.887062599768095, -0.887062599768095 }, {
// 0.978228658146057, 0.887062599768095, -0.730152005574049 }, {
// 0.978228658146057, 0.887062599768095, -0.519096129206812 }, {
// 0.978228658146057, 0.887062599768095, -0.269543155952345 }, {
// 0.978228658146057, 0.887062599768095, 0 }, { 0.978228658146057,
// 0.887062599768095, 0.269543155952345 }, { 0.978228658146057,
// 0.887062599768095, 0.519096129206812 }, { 0.978228658146057,
// 0.887062599768095, 0.730152005574049 }, { 0.978228658146057,
// 0.887062599768095, 0.887062599768095 }, { 0.978228658146057,
// 0.887062599768095, 0.978228658146057 }, { 0.978228658146057,
// 0.978228658146057, -0.978228658146057 }, { 0.978228658146057,
// 0.978228658146057, -0.887062599768095 }, { 0.978228658146057,
// 0.978228658146057, -0.730152005574049 }, { 0.978228658146057,
// 0.978228658146057, -0.519096129206812 }, { 0.978228658146057,
// 0.978228658146057, -0.269543155952345 }, { 0.978228658146057,
// 0.978228658146057, 0 }, { 0.978228658146057, 0.978228658146057,
// 0.269543155952345 }, { 0.978228658146057, 0.978228658146057,
// 0.519096129206812
//}, { 0.978228658146057, 0.978228658146057, 0.730152005574049 }, {
// 0.978228658146057, 0.978228658146057, 0.887062599768095 }, {
// 0.978228658146057, 0.978228658146057, 0.978228658146057 } } , {
// 0.0001725162974448952, 0.0003891722293953698, 0.0005773113824254295,
// 0.0007226649964007826, 0.0008144284884499459, 0.0008457919413030904,
// 0.0008144284884499459, 0.0007226649964007826, 0.0005773113824254295,
// 0.0003891722293953698, 0.0001725162974448952, 0.0003891722293953698,
// 0.0008779171960894872, 0.001302332365587617, 0.001630229444526092,
// 0.001837234830723208, 0.001907986319418319, 0.001837234830723208,
// 0.001630229444526092, 0.001302332365587617, 0.0008779171960894872,
// 0.0003891722293953698, 0.0005773113824254295, 0.001302332365587617,
// 0.001931924329551641, 0.002418338060123659, 0.002725416923023606,
// 0.002830372098809633, 0.002725416923023606, 0.002418338060123659,
// 0.001931924329551641, 0.001302332365587617, 0.0005773113824254295,
// 0.0007226649964007826, 0.001630229444526091, 0.002418338060123659,
// 0.003027219484522948, 0.003411613681671851, 0.003542994136034277,
// 0.003411613681671851, 0.003027219484522948, 0.002418338060123659,
// 0.001630229444526091, 0.0007226649964007826, 0.0008144284884499458,
// 0.001837234830723208, 0.002725416923023605, 0.003411613681671851,
// 0.003844817983128415, 0.003992881035014378, 0.003844817983128415,
// 0.003411613681671851, 0.002725416923023605, 0.001837234830723208,
// 0.0008144284884499458, 0.0008457919413030904, 0.001907986319418319,
// 0.002830372098809633, 0.003542994136034277, 0.003992881035014378,
// 0.004146645960806982, 0.003992881035014378, 0.003542994136034277,
// 0.002830372098809633, 0.001907986319418319, 0.0008457919413030904,
// 0.0008144284884499458, 0.001837234830723208, 0.002725416923023605,
// 0.003411613681671851, 0.003844817983128415, 0.003992881035014378,
// 0.003844817983128415, 0.003411613681671851, 0.002725416923023605,
// 0.001837234830723208, 0.0008144284884499458, 0.0007226649964007826,
// 0.001630229444526091, 0.002418338060123659, 0.003027219484522948,
// 0.003411613681671851, 0.003542994136034277, 0.003411613681671851,
// 0.003027219484522948, 0.002418338060123659, 0.001630229444526091,
// 0.0007226649964007826, 0.0005773113824254295, 0.001302332365587617,
// 0.001931924329551641, 0.002418338060123659, 0.002725416923023606,
// 0.002830372098809633, 0.002725416923023606, 0.002418338060123659,
// 0.001931924329551641, 0.001302332365587617, 0.0005773113824254295,
// 0.0003891722293953698, 0.0008779171960894872, 0.001302332365587617,
// 0.001630229444526092, 0.001837234830723208, 0.001907986319418319,
// 0.001837234830723208, 0.001630229444526092, 0.001302332365587617,
// 0.0008779171960894872, 0.0003891722293953698, 0.0001725162974448952,
// 0.0003891722293953698, 0.0005773113824254295, 0.0007226649964007826,
// 0.0008144284884499459, 0.0008457919413030904, 0.0008144284884499459,
// 0.0007226649964007826, 0.0005773113824254295, 0.0003891722293953698,
// 0.0001725162974448952, 0.0003891722293953698, 0.0008779171960894872,
// 0.001302332365587617, 0.001630229444526092, 0.001837234830723208,
// 0.001907986319418319, 0.001837234830723208, 0.001630229444526092,
// 0.001302332365587617, 0.0008779171960894872, 0.0003891722293953698,
// 0.0008779171960894872, 0.001980456324920899, 0.002937876581146601,
// 0.003677565753199903, 0.00414454046130778, 0.004304145756554233,
// 0.00414454046130778, 0.003677565753199903, 0.002937876581146601,
// 0.001980456324920899, 0.0008779171960894872, 0.001302332365587617,
// 0.002937876581146601, 0.004358146502622003, 0.00545542669434242,
// 0.006148152932065786, 0.006384916880470592, 0.006148152932065786,
// 0.00545542669434242, 0.004358146502622003, 0.002937876581146601,
// 0.001302332365587617, 0.001630229444526091, 0.003677565753199903,
// 0.005455426694342419, 0.006828976584297535, 0.007696115219236448,
// 0.007992490837543083, 0.007696115219236448, 0.006828976584297535,
// 0.005455426694342419, 0.003677565753199903, 0.001630229444526091,
// 0.001837234830723208, 0.004144540461307779, 0.006148152932065785,
// 0.007696115219236448, 0.008673362507049162, 0.00900737169254045,
// 0.008673362507049162, 0.007696115219236448, 0.006148152932065785,
// 0.004144540461307779, 0.001837234830723208, 0.001907986319418319,
// 0.004304145756554232, 0.006384916880470592, 0.007992490837543083,
// 0.00900737169254045, 0.0093542434945662, 0.00900737169254045,
// 0.007992490837543083, 0.006384916880470592, 0.004304145756554232,
// 0.001907986319418319, 0.001837234830723208, 0.004144540461307779,
// 0.006148152932065785, 0.007696115219236448, 0.008673362507049162,
// 0.00900737169254045, 0.008673362507049162, 0.007696115219236448,
// 0.006148152932065785, 0.004144540461307779, 0.001837234830723208,
// 0.001630229444526091, 0.003677565753199903, 0.005455426694342419,
// 0.006828976584297535, 0.007696115219236448, 0.007992490837543083,
// 0.007696115219236448, 0.006828976584297535, 0.005455426694342419,
// 0.003677565753199903, 0.001630229444526091, 0.001302332365587617,
// 0.002937876581146601, 0.004358146502622003, 0.00545542669434242,
// 0.006148152932065786, 0.006384916880470592, 0.006148152932065786,
// 0.00545542669434242, 0.004358146502622003, 0.002937876581146601,
// 0.001302332365587617, 0.0008779171960894872, 0.001980456324920899,
// 0.002937876581146601, 0.003677565753199903, 0.00414454046130778,
// 0.004304145756554233, 0.00414454046130778, 0.003677565753199903,
// 0.002937876581146601, 0.001980456324920899, 0.0008779171960894872,
// 0.0003891722293953698, 0.0008779171960894872, 0.001302332365587617,
// 0.001630229444526092, 0.001837234830723208, 0.001907986319418319,
// 0.001837234830723208, 0.001630229444526092, 0.001302332365587617,
// 0.0008779171960894872, 0.0003891722293953698, 0.0005773113824254295,
// 0.001302332365587617, 0.001931924329551641, 0.002418338060123659,
// 0.002725416923023606, 0.002830372098809633, 0.002725416923023606,
// 0.002418338060123659, 0.001931924329551641, 0.001302332365587617,
// 0.0005773113824254295, 0.001302332365587617, 0.002937876581146601,
// 0.004358146502622003, 0.00545542669434242, 0.006148152932065786,
// 0.006384916880470592, 0.006148152932065786, 0.00545542669434242,
// 0.004358146502622003, 0.002937876581146601, 0.001302332365587617,
// 0.001931924329551641, 0.004358146502622003, 0.006465023432299391,
// 0.008092766360859161, 0.009120380131152508, 0.009471603862029818,
// 0.009120380131152508, 0.008092766360859161, 0.006465023432299391,
// 0.004358146502622003, 0.001931924329551641, 0.002418338060123659,
// 0.00545542669434242, 0.008092766360859161, 0.01013033719943688,
// 0.01141668027912899, 0.01185633399796608, 0.01141668027912899,
// 0.01013033719943688, 0.008092766360859161, 0.00545542669434242,
// 0.002418338060123659, 0.002725416923023606, 0.006148152932065786,
// 0.009120380131152508, 0.01141668027912899, 0.01286636229671586,
// 0.0133618429349884, 0.01286636229671586, 0.01141668027912899,
// 0.009120380131152508, 0.006148152932065786, 0.002725416923023606,
// 0.002830372098809633, 0.006384916880470591, 0.009471603862029816,
// 0.01185633399796608, 0.0133618429349884, 0.01387640441812148,
// 0.0133618429349884, 0.01185633399796608, 0.009471603862029816,
// 0.006384916880470591, 0.002830372098809633, 0.002725416923023606,
// 0.006148152932065786, 0.009120380131152508, 0.01141668027912899,
// 0.01286636229671586, 0.0133618429349884, 0.01286636229671586,
// 0.01141668027912899, 0.009120380131152508, 0.006148152932065786,
// 0.002725416923023606, 0.002418338060123659, 0.00545542669434242,
// 0.008092766360859161, 0.01013033719943688, 0.01141668027912899,
// 0.01185633399796608, 0.01141668027912899, 0.01013033719943688,
// 0.008092766360859161, 0.00545542669434242, 0.002418338060123659,
// 0.001931924329551641, 0.004358146502622003, 0.006465023432299391,
// 0.008092766360859161, 0.009120380131152508, 0.009471603862029818,
// 0.009120380131152508, 0.008092766360859161, 0.006465023432299391,
// 0.004358146502622003, 0.001931924329551641, 0.001302332365587617,
// 0.002937876581146601, 0.004358146502622003, 0.00545542669434242,
// 0.006148152932065786, 0.006384916880470592, 0.006148152932065786,
// 0.00545542669434242, 0.004358146502622003, 0.002937876581146601,
// 0.001302332365587617, 0.0005773113824254295, 0.001302332365587617,
// 0.001931924329551641, 0.002418338060123659, 0.002725416923023606,
// 0.002830372098809633, 0.002725416923023606, 0.002418338060123659,
// 0.001931924329551641, 0.001302332365587617, 0.0005773113824254295,
// 0.0007226649964007826, 0.001630229444526091, 0.002418338060123659,
// 0.003027219484522948, 0.003411613681671851, 0.003542994136034277,
// 0.003411613681671851, 0.003027219484522948, 0.002418338060123659,
// 0.001630229444526091, 0.0007226649964007826, 0.001630229444526091,
// 0.003677565753199903, 0.005455426694342419, 0.006828976584297535,
// 0.007696115219236448, 0.007992490837543083, 0.007696115219236448,
// 0.006828976584297535, 0.005455426694342419, 0.003677565753199903,
// 0.001630229444526091, 0.002418338060123659, 0.00545542669434242,
// 0.008092766360859161, 0.01013033719943688, 0.01141668027912899,
// 0.01185633399796608, 0.01141668027912899, 0.01013033719943688,
// 0.008092766360859161, 0.00545542669434242, 0.002418338060123659,
// 0.003027219484522948, 0.006828976584297535, 0.01013033719943688,
// 0.01268092110883573, 0.01429113553618759, 0.01484148386260749,
// 0.01429113553618759, 0.01268092110883573, 0.01013033719943688,
// 0.006828976584297535, 0.003027219484522948, 0.003411613681671851,
// 0.007696115219236449, 0.01141668027912899, 0.01429113553618759,
// 0.01610581385695847, 0.01672604502608865, 0.01610581385695847,
// 0.01429113553618759, 0.01141668027912899, 0.007696115219236449,
// 0.003411613681671851, 0.003542994136034277, 0.007992490837543083,
// 0.01185633399796608, 0.01484148386260749, 0.01672604502608865,
// 0.01737016115418938, 0.01672604502608865, 0.01484148386260749,
// 0.01185633399796608, 0.007992490837543083, 0.003542994136034277,
// 0.003411613681671851, 0.007696115219236449, 0.01141668027912899,
// 0.01429113553618759, 0.01610581385695847, 0.01672604502608865,
// 0.01610581385695847, 0.01429113553618759, 0.01141668027912899,
// 0.007696115219236449, 0.003411613681671851, 0.003027219484522948,
// 0.006828976584297535, 0.01013033719943688, 0.01268092110883573,
// 0.01429113553618759, 0.01484148386260749, 0.01429113553618759,
// 0.01268092110883573, 0.01013033719943688, 0.006828976584297535,
// 0.003027219484522948, 0.002418338060123659, 0.00545542669434242,
// 0.008092766360859161, 0.01013033719943688, 0.01141668027912899,
// 0.01185633399796608, 0.01141668027912899, 0.01013033719943688,
// 0.008092766360859161, 0.00545542669434242, 0.002418338060123659,
// 0.001630229444526091, 0.003677565753199903, 0.005455426694342419,
// 0.006828976584297535, 0.007696115219236448, 0.007992490837543083,
// 0.007696115219236448, 0.006828976584297535, 0.005455426694342419,
// 0.003677565753199903, 0.001630229444526091, 0.0007226649964007826,
// 0.001630229444526091, 0.002418338060123659, 0.003027219484522948,
// 0.003411613681671851, 0.003542994136034277, 0.003411613681671851,
// 0.003027219484522948, 0.002418338060123659, 0.001630229444526091,
// 0.0007226649964007826, 0.0008144284884499458, 0.001837234830723208,
// 0.002725416923023605, 0.003411613681671851, 0.003844817983128415,
// 0.003992881035014378, 0.003844817983128415, 0.003411613681671851,
// 0.002725416923023605, 0.001837234830723208, 0.0008144284884499458,
// 0.001837234830723208, 0.004144540461307779, 0.006148152932065785,
// 0.007696115219236448, 0.008673362507049162, 0.00900737169254045,
// 0.008673362507049162, 0.007696115219236448, 0.006148152932065785,
// 0.004144540461307779, 0.001837234830723208, 0.002725416923023606,
// 0.006148152932065786, 0.009120380131152508, 0.01141668027912899,
// 0.01286636229671586, 0.0133618429349884, 0.01286636229671586,
// 0.01141668027912899, 0.009120380131152508, 0.006148152932065786,
// 0.002725416923023606, 0.003411613681671851, 0.007696115219236449,
// 0.01141668027912899, 0.01429113553618759, 0.01610581385695847,
// 0.01672604502608865, 0.01610581385695847, 0.01429113553618759,
// 0.01141668027912899, 0.007696115219236449, 0.003411613681671851,
// 0.003844817983128416, 0.008673362507049164, 0.01286636229671586,
// 0.01610581385695847, 0.01815091875226831, 0.01884990643823629,
// 0.01815091875226831, 0.01610581385695847, 0.01286636229671586,
// 0.008673362507049164, 0.003844817983128416, 0.003992881035014378,
// 0.00900737169254045, 0.0133618429349884, 0.01672604502608865,
// 0.01884990643823628, 0.0195758119784354, 0.01884990643823628,
// 0.01672604502608865, 0.0133618429349884, 0.00900737169254045,
// 0.003992881035014378, 0.003844817983128416, 0.008673362507049164,
// 0.01286636229671586, 0.01610581385695847, 0.01815091875226831,
// 0.01884990643823629, 0.01815091875226831, 0.01610581385695847,
// 0.01286636229671586, 0.008673362507049164, 0.003844817983128416,
// 0.003411613681671851, 0.007696115219236449, 0.01141668027912899,
// 0.01429113553618759, 0.01610581385695847, 0.01672604502608865,
// 0.01610581385695847, 0.01429113553618759, 0.01141668027912899,
// 0.007696115219236449, 0.003411613681671851, 0.002725416923023606,
// 0.006148152932065786, 0.009120380131152508, 0.01141668027912899,
// 0.01286636229671586, 0.0133618429349884, 0.01286636229671586,
// 0.01141668027912899, 0.009120380131152508, 0.006148152932065786,
// 0.002725416923023606, 0.001837234830723208, 0.004144540461307779,
// 0.006148152932065785, 0.007696115219236448, 0.008673362507049162,
// 0.00900737169254045, 0.008673362507049162, 0.007696115219236448,
// 0.006148152932065785, 0.004144540461307779, 0.001837234830723208,
// 0.0008144284884499458, 0.001837234830723208, 0.002725416923023605,
// 0.003411613681671851, 0.003844817983128415, 0.003992881035014378,
// 0.003844817983128415, 0.003411613681671851, 0.002725416923023605,
// 0.001837234830723208, 0.0008144284884499458, 0.0008457919413030904,
// 0.001907986319418319, 0.002830372098809633, 0.003542994136034277,
// 0.003992881035014378, 0.004146645960806982, 0.003992881035014378,
// 0.003542994136034277, 0.002830372098809633, 0.001907986319418319,
// 0.0008457919413030904, 0.001907986319418319, 0.004304145756554232,
// 0.006384916880470592, 0.007992490837543083, 0.00900737169254045,
// 0.0093542434945662, 0.00900737169254045, 0.007992490837543083,
// 0.006384916880470592, 0.004304145756554232, 0.001907986319418319,
// 0.002830372098809633, 0.006384916880470591, 0.009471603862029816,
// 0.01185633399796608, 0.0133618429349884, 0.01387640441812148,
// 0.0133618429349884, 0.01185633399796608, 0.009471603862029816,
// 0.006384916880470591, 0.002830372098809633, 0.003542994136034277,
// 0.007992490837543083, 0.01185633399796608, 0.01484148386260749,
// 0.01672604502608865, 0.01737016115418938, 0.01672604502608865,
// 0.01484148386260749, 0.01185633399796608, 0.007992490837543083,
// 0.003542994136034277, 0.003992881035014378, 0.00900737169254045,
// 0.0133618429349884, 0.01672604502608865, 0.01884990643823628,
// 0.0195758119784354, 0.01884990643823628, 0.01672604502608865,
// 0.0133618429349884, 0.00900737169254045, 0.003992881035014378,
// 0.004146645960806982, 0.0093542434945662, 0.01387640441812148,
// 0.01737016115418937, 0.0195758119784354, 0.02032967197321064,
// 0.0195758119784354, 0.01737016115418937, 0.01387640441812148,
// 0.0093542434945662, 0.004146645960806982, 0.003992881035014378,
// 0.00900737169254045, 0.0133618429349884, 0.01672604502608865,
// 0.01884990643823628, 0.0195758119784354, 0.01884990643823628,
// 0.01672604502608865, 0.0133618429349884, 0.00900737169254045,
// 0.003992881035014378, 0.003542994136034277, 0.007992490837543083,
// 0.01185633399796608, 0.01484148386260749, 0.01672604502608865,
// 0.01737016115418938, 0.01672604502608865, 0.01484148386260749,
// 0.01185633399796608, 0.007992490837543083, 0.003542994136034277,
// 0.002830372098809633, 0.006384916880470591, 0.009471603862029816,
// 0.01185633399796608, 0.0133618429349884, 0.01387640441812148,
// 0.0133618429349884, 0.01185633399796608, 0.009471603862029816,
// 0.006384916880470591, 0.002830372098809633, 0.001907986319418319,
// 0.004304145756554232, 0.006384916880470592, 0.007992490837543083,
// 0.00900737169254045, 0.0093542434945662, 0.00900737169254045,
// 0.007992490837543083, 0.006384916880470592, 0.004304145756554232,
// 0.001907986319418319, 0.0008457919413030904, 0.001907986319418319,
// 0.002830372098809633, 0.003542994136034277, 0.003992881035014378,
// 0.004146645960806982, 0.003992881035014378, 0.003542994136034277,
// 0.002830372098809633, 0.001907986319418319, 0.0008457919413030904,
// 0.0008144284884499458, 0.001837234830723208, 0.002725416923023605,
// 0.003411613681671851, 0.003844817983128415, 0.003992881035014378,
// 0.003844817983128415, 0.003411613681671851, 0.002725416923023605,
// 0.001837234830723208, 0.0008144284884499458, 0.001837234830723208,
// 0.004144540461307779, 0.006148152932065785, 0.007696115219236448,
// 0.008673362507049162, 0.00900737169254045, 0.008673362507049162,
// 0.007696115219236448, 0.006148152932065785, 0.004144540461307779,
// 0.001837234830723208, 0.002725416923023606, 0.006148152932065786,
// 0.009120380131152508, 0.01141668027912899, 0.01286636229671586,
// 0.0133618429349884, 0.01286636229671586, 0.01141668027912899,
// 0.009120380131152508, 0.006148152932065786, 0.002725416923023606,
// 0.003411613681671851, 0.007696115219236449, 0.01141668027912899,
// 0.01429113553618759, 0.01610581385695847, 0.01672604502608865,
// 0.01610581385695847, 0.01429113553618759, 0.01141668027912899,
// 0.007696115219236449, 0.003411613681671851, 0.003844817983128416,
// 0.008673362507049164, 0.01286636229671586, 0.01610581385695847,
// 0.01815091875226831, 0.01884990643823629, 0.01815091875226831,
// 0.01610581385695847, 0.01286636229671586, 0.008673362507049164,
// 0.003844817983128416, 0.003992881035014378, 0.00900737169254045,
// 0.0133618429349884, 0.01672604502608865, 0.01884990643823628,
// 0.0195758119784354, 0.01884990643823628, 0.01672604502608865,
// 0.0133618429349884, 0.00900737169254045, 0.003992881035014378,
// 0.003844817983128416, 0.008673362507049164, 0.01286636229671586,
// 0.01610581385695847, 0.01815091875226831, 0.01884990643823629,
// 0.01815091875226831, 0.01610581385695847, 0.01286636229671586,
// 0.008673362507049164, 0.003844817983128416, 0.003411613681671851,
// 0.007696115219236449, 0.01141668027912899, 0.01429113553618759,
// 0.01610581385695847, 0.01672604502608865, 0.01610581385695847,
// 0.01429113553618759, 0.01141668027912899, 0.007696115219236449,
// 0.003411613681671851, 0.002725416923023606, 0.006148152932065786,
// 0.009120380131152508, 0.01141668027912899, 0.01286636229671586,
// 0.0133618429349884, 0.01286636229671586, 0.01141668027912899,
// 0.009120380131152508, 0.006148152932065786, 0.002725416923023606,
// 0.001837234830723208, 0.004144540461307779, 0.006148152932065785,
// 0.007696115219236448, 0.008673362507049162, 0.00900737169254045,
// 0.008673362507049162, 0.007696115219236448, 0.006148152932065785,
// 0.004144540461307779, 0.001837234830723208, 0.0008144284884499458,
// 0.001837234830723208, 0.002725416923023605, 0.003411613681671851,
// 0.003844817983128415, 0.003992881035014378, 0.003844817983128415,
// 0.003411613681671851, 0.002725416923023605, 0.001837234830723208,
// 0.0008144284884499458, 0.0007226649964007826, 0.001630229444526091,
// 0.002418338060123659, 0.003027219484522948, 0.003411613681671851,
// 0.003542994136034277, 0.003411613681671851, 0.003027219484522948,
// 0.002418338060123659, 0.001630229444526091, 0.0007226649964007826,
// 0.001630229444526091, 0.003677565753199903, 0.005455426694342419,
// 0.006828976584297535, 0.007696115219236448, 0.007992490837543083,
// 0.007696115219236448, 0.006828976584297535, 0.005455426694342419,
// 0.003677565753199903, 0.001630229444526091, 0.002418338060123659,
// 0.00545542669434242, 0.008092766360859161, 0.01013033719943688,
// 0.01141668027912899, 0.01185633399796608, 0.01141668027912899,
// 0.01013033719943688, 0.008092766360859161, 0.00545542669434242,
// 0.002418338060123659, 0.003027219484522948, 0.006828976584297535,
// 0.01013033719943688, 0.01268092110883573, 0.01429113553618759,
// 0.01484148386260749, 0.01429113553618759, 0.01268092110883573,
// 0.01013033719943688, 0.006828976584297535, 0.003027219484522948,
// 0.003411613681671851, 0.007696115219236449, 0.01141668027912899,
// 0.01429113553618759, 0.01610581385695847, 0.01672604502608865,
// 0.01610581385695847, 0.01429113553618759, 0.01141668027912899,
// 0.007696115219236449, 0.003411613681671851, 0.003542994136034277,
// 0.007992490837543083, 0.01185633399796608, 0.01484148386260749,
// 0.01672604502608865, 0.01737016115418938, 0.01672604502608865,
// 0.01484148386260749, 0.01185633399796608, 0.007992490837543083,
// 0.003542994136034277, 0.003411613681671851, 0.007696115219236449,
// 0.01141668027912899, 0.01429113553618759, 0.01610581385695847,
// 0.01672604502608865, 0.01610581385695847, 0.01429113553618759,
// 0.01141668027912899, 0.007696115219236449, 0.003411613681671851,
// 0.003027219484522948, 0.006828976584297535, 0.01013033719943688,
// 0.01268092110883573, 0.01429113553618759, 0.01484148386260749,
// 0.01429113553618759, 0.01268092110883573, 0.01013033719943688,
// 0.006828976584297535, 0.003027219484522948, 0.002418338060123659,
// 0.00545542669434242, 0.008092766360859161, 0.01013033719943688,
// 0.01141668027912899, 0.01185633399796608, 0.01141668027912899,
// 0.01013033719943688, 0.008092766360859161, 0.00545542669434242,
// 0.002418338060123659, 0.001630229444526091, 0.003677565753199903,
// 0.005455426694342419, 0.006828976584297535, 0.007696115219236448,
// 0.007992490837543083, 0.007696115219236448, 0.006828976584297535,
// 0.005455426694342419, 0.003677565753199903, 0.001630229444526091,
// 0.0007226649964007826, 0.001630229444526091, 0.002418338060123659,
// 0.003027219484522948, 0.003411613681671851, 0.003542994136034277,
// 0.003411613681671851, 0.003027219484522948, 0.002418338060123659,
// 0.001630229444526091, 0.0007226649964007826, 0.0005773113824254295,
// 0.001302332365587617, 0.001931924329551641, 0.002418338060123659,
// 0.002725416923023606, 0.002830372098809633, 0.002725416923023606,
// 0.002418338060123659, 0.001931924329551641, 0.001302332365587617,
// 0.0005773113824254295, 0.001302332365587617, 0.002937876581146601,
// 0.004358146502622003, 0.00545542669434242, 0.006148152932065786,
// 0.006384916880470592, 0.006148152932065786, 0.00545542669434242,
// 0.004358146502622003, 0.002937876581146601, 0.001302332365587617,
// 0.001931924329551641, 0.004358146502622003, 0.006465023432299391,
// 0.008092766360859161, 0.009120380131152508, 0.009471603862029818,
// 0.009120380131152508, 0.008092766360859161, 0.006465023432299391,
// 0.004358146502622003, 0.001931924329551641, 0.002418338060123659,
// 0.00545542669434242, 0.008092766360859161, 0.01013033719943688,
// 0.01141668027912899, 0.01185633399796608, 0.01141668027912899,
// 0.01013033719943688, 0.008092766360859161, 0.00545542669434242,
// 0.002418338060123659, 0.002725416923023606, 0.006148152932065786,
// 0.009120380131152508, 0.01141668027912899, 0.01286636229671586,
// 0.0133618429349884, 0.01286636229671586, 0.01141668027912899,
// 0.009120380131152508, 0.006148152932065786, 0.002725416923023606,
// 0.002830372098809633, 0.006384916880470591, 0.009471603862029816,
// 0.01185633399796608, 0.0133618429349884, 0.01387640441812148,
// 0.0133618429349884, 0.01185633399796608, 0.009471603862029816,
// 0.006384916880470591, 0.002830372098809633, 0.002725416923023606,
// 0.006148152932065786, 0.009120380131152508, 0.01141668027912899,
// 0.01286636229671586, 0.0133618429349884, 0.01286636229671586,
// 0.01141668027912899, 0.009120380131152508, 0.006148152932065786,
// 0.002725416923023606, 0.002418338060123659, 0.00545542669434242,
// 0.008092766360859161, 0.01013033719943688, 0.01141668027912899,
// 0.01185633399796608, 0.01141668027912899, 0.01013033719943688,
// 0.008092766360859161, 0.00545542669434242, 0.002418338060123659,
// 0.001931924329551641, 0.004358146502622003, 0.006465023432299391,
// 0.008092766360859161, 0.009120380131152508, 0.009471603862029818,
// 0.009120380131152508, 0.008092766360859161, 0.006465023432299391,
// 0.004358146502622003, 0.001931924329551641, 0.001302332365587617,
// 0.002937876581146601, 0.004358146502622003, 0.00545542669434242,
// 0.006148152932065786, 0.006384916880470592, 0.006148152932065786,
// 0.00545542669434242, 0.004358146502622003, 0.002937876581146601,
// 0.001302332365587617, 0.0005773113824254295, 0.001302332365587617,
// 0.001931924329551641, 0.002418338060123659, 0.002725416923023606,
// 0.002830372098809633, 0.002725416923023606, 0.002418338060123659,
// 0.001931924329551641, 0.001302332365587617, 0.0005773113824254295,
// 0.0003891722293953698, 0.0008779171960894872, 0.001302332365587617,
// 0.001630229444526092, 0.001837234830723208, 0.001907986319418319,
// 0.001837234830723208, 0.001630229444526092, 0.001302332365587617,
// 0.0008779171960894872, 0.0003891722293953698, 0.0008779171960894872,
// 0.001980456324920899, 0.002937876581146601, 0.003677565753199903,
// 0.00414454046130778, 0.004304145756554233, 0.00414454046130778,
// 0.003677565753199903, 0.002937876581146601, 0.001980456324920899,
// 0.0008779171960894872, 0.001302332365587617, 0.002937876581146601,
// 0.004358146502622003, 0.00545542669434242, 0.006148152932065786,
// 0.006384916880470592, 0.006148152932065786, 0.00545542669434242,
// 0.004358146502622003, 0.002937876581146601, 0.001302332365587617,
// 0.001630229444526091, 0.003677565753199903, 0.005455426694342419,
// 0.006828976584297535, 0.007696115219236448, 0.007992490837543083,
// 0.007696115219236448, 0.006828976584297535, 0.005455426694342419,
// 0.003677565753199903, 0.001630229444526091, 0.001837234830723208,
// 0.004144540461307779, 0.006148152932065785, 0.007696115219236448,
// 0.008673362507049162, 0.00900737169254045, 0.008673362507049162,
// 0.007696115219236448, 0.006148152932065785, 0.004144540461307779,
// 0.001837234830723208, 0.001907986319418319, 0.004304145756554232,
// 0.006384916880470592, 0.007992490837543083, 0.00900737169254045,
// 0.0093542434945662, 0.00900737169254045, 0.007992490837543083,
// 0.006384916880470592, 0.004304145756554232, 0.001907986319418319,
// 0.001837234830723208, 0.004144540461307779, 0.006148152932065785,
// 0.007696115219236448, 0.008673362507049162, 0.00900737169254045,
// 0.008673362507049162, 0.007696115219236448, 0.006148152932065785,
// 0.004144540461307779, 0.001837234830723208, 0.001630229444526091,
// 0.003677565753199903, 0.005455426694342419, 0.006828976584297535,
// 0.007696115219236448, 0.007992490837543083, 0.007696115219236448,
// 0.006828976584297535, 0.005455426694342419, 0.003677565753199903,
// 0.001630229444526091, 0.001302332365587617, 0.002937876581146601,
// 0.004358146502622003, 0.00545542669434242, 0.006148152932065786,
// 0.006384916880470592, 0.006148152932065786, 0.00545542669434242,
// 0.004358146502622003, 0.002937876581146601, 0.001302332365587617,
// 0.0008779171960894872, 0.001980456324920899, 0.002937876581146601,
// 0.003677565753199903, 0.00414454046130778, 0.004304145756554233,
// 0.00414454046130778, 0.003677565753199903, 0.002937876581146601,
// 0.001980456324920899, 0.0008779171960894872, 0.0003891722293953698,
// 0.0008779171960894872, 0.001302332365587617, 0.001630229444526092,
// 0.001837234830723208, 0.001907986319418319, 0.001837234830723208,
// 0.001630229444526092, 0.001302332365587617, 0.0008779171960894872,
// 0.0003891722293953698, 0.0001725162974448952, 0.0003891722293953698,
// 0.0005773113824254295, 0.0007226649964007826, 0.0008144284884499459,
// 0.0008457919413030904, 0.0008144284884499459, 0.0007226649964007826,
// 0.0005773113824254295, 0.0003891722293953698, 0.0001725162974448952,
// 0.0003891722293953698, 0.0008779171960894872, 0.001302332365587617,
// 0.001630229444526092, 0.001837234830723208, 0.001907986319418319,
// 0.001837234830723208, 0.001630229444526092, 0.001302332365587617,
// 0.0008779171960894872, 0.0003891722293953698, 0.0005773113824254295,
// 0.001302332365587617, 0.001931924329551641, 0.002418338060123659,
// 0.002725416923023606, 0.002830372098809633, 0.002725416923023606,
// 0.002418338060123659, 0.001931924329551641, 0.001302332365587617,
// 0.0005773113824254295, 0.0007226649964007826, 0.001630229444526091,
// 0.002418338060123659, 0.003027219484522948, 0.003411613681671851,
// 0.003542994136034277, 0.003411613681671851, 0.003027219484522948,
// 0.002418338060123659, 0.001630229444526091, 0.0007226649964007826,
// 0.0008144284884499458, 0.001837234830723208, 0.002725416923023605,
// 0.003411613681671851, 0.003844817983128415, 0.003992881035014378,
// 0.003844817983128415, 0.003411613681671851, 0.002725416923023605,
// 0.001837234830723208, 0.0008144284884499458, 0.0008457919413030904,
// 0.001907986319418319, 0.002830372098809633, 0.003542994136034277,
// 0.003992881035014378, 0.004146645960806982, 0.003992881035014378,
// 0.003542994136034277, 0.002830372098809633, 0.001907986319418319,
// 0.0008457919413030904, 0.0008144284884499458, 0.001837234830723208,
// 0.002725416923023605, 0.003411613681671851, 0.003844817983128415,
// 0.003992881035014378, 0.003844817983128415, 0.003411613681671851,
// 0.002725416923023605, 0.001837234830723208, 0.0008144284884499458,
// 0.0007226649964007826, 0.001630229444526091, 0.002418338060123659,
// 0.003027219484522948, 0.003411613681671851, 0.003542994136034277,
// 0.003411613681671851, 0.003027219484522948, 0.002418338060123659,
// 0.001630229444526091, 0.0007226649964007826, 0.0005773113824254295,
// 0.001302332365587617, 0.001931924329551641, 0.002418338060123659,
// 0.002725416923023606, 0.002830372098809633, 0.002725416923023606,
// 0.002418338060123659, 0.001931924329551641, 0.001302332365587617,
// 0.0005773113824254295, 0.0003891722293953698, 0.0008779171960894872,
// 0.001302332365587617, 0.001630229444526092, 0.001837234830723208,
// 0.001907986319418319, 0.001837234830723208, 0.001630229444526092,
// 0.001302332365587617, 0.0008779171960894872, 0.0003891722293953698,
// 0.0001725162974448952, 0.0003891722293953698, 0.0005773113824254295,
// 0.0007226649964007826, 0.0008144284884499459, 0.0008457919413030904,
// 0.0008144284884499459, 0.0007226649964007826, 0.0005773113824254295,
// 0.0003891722293953698, 0.0001725162974448952 } };
// default: 						throw
// std::runtime_error("not supported order");
//				}
//}

// std::vector<Euclidean_Vector<space_dimension>>
// ReferenceGeometry<space_dimension>::reference_post_nodes(const int
// post_order) const { 	const auto n = post_order; 			case
// Figure::tetrahedral: {
//				const auto num_reference_post_point = (n + 2) *
//(n
//+
// 3) * (n + 4) / 6.0; 				const auto delta = 2.0 / (n +
// 1);
//
//				const auto x0_start_coord = -1.0;
//				const auto x1_start_coord = -1.0;
//				const auto x2_start_coord = -1.0;
//
//				for (int i = 0; i <= n + 1; ++i) {
//					for (int j = 0; j <= n + 1; ++j) {
//						for (int k = 0; k <= n + 1 - i -
// j;
//++k) { 							const double
// x0_coord = x0_start_coord + delta * k;
// const double x1_coord = x1_start_coord + delta * j;
// const double x2_coord = x2_start_coord
//+ delta * i;
//
//							reference_post_nodes.push_back({
// x0_coord, x1_coord, x2_coord });
//						}
//					}
//				}
//				break;
//			}
//			case Figure::hexahedral: {
//				const auto num_reference_post_point = (n + 2) *
//(n
//+ 2) * (n + 2); 				const auto delta = 2.0 / (n +
// 1);
//
//				const auto x0_start_coord = -1.0;
//				const auto x1_start_coord = -1.0;
//				const auto x2_start_coord = -1.0;
//
//				for (int i = 0; i <= n + 1; ++i) {
//					for (int j = 0; j <= n + 1; ++j) {
//						for (int k = 0; k <= n + 1; ++k)
//{ 							const double x0_coord =
// x0_start_coord + delta * k; const double x1_coord = x1_start_coord + delta *
// j;
// const double x2_coord = x2_start_coord + delta * i;
//
//							reference_post_nodes.push_back({
// x0_coord, x1_coord, x2_coord });
//						}
//					}
//				}
//
//	return reference_post_nodes;
// }

// std::vector<std::vector<uint>>
// ReferenceGeometry<space_dimension>::reference_connectivity(const int
// post_order) const { 	const auto n = post_order;

//		case Figure::tetrahedral: {
//			const auto num_simplex = (n + 1) * (n + 1) * (n + 1);
//			reference_connectivities.reserve(num_simplex);
//
//			for (int iz = 0; iz <= n; iz++) {
//				for (int iy = 0; iy <= n - iz; iy++) {
//					for (int ix = 0; ix <= n - iz - iy;
// ix++) { 						const int sum = ix + iy
//+ iz;
//						//	     I3��������I4
//						//      /��    /��
//						//     I3+1����I4+1
//						//     �� I1��������I2
//						//	   ��/    ��/
//						//     I1+1����I2+1
//
//						const uint I1 = ((n + 1) - iz +
// 2)
//* iy
//- iy * (iy + 1) / 2 + ix + iz * (3 * (n + 1) * (n + 1) - 3 * (n + 1) * iz + 12
// *
//(n + 1) + iz * iz - 6 * iz + 11) / 6;
// const uint I2 = ((n + 1) - iz + 2) * (iy
//+ 1) - (iy + 1) * (iy + 2) / 2 + ix + iz * (3 * (n + 1) * (n + 1) - 3 * (n +
// 1) * iz + 12 * (n + 1) + iz * iz - 6 * iz + 11) / 6;
// const uint I3 = ((n + 1)
//- iz + 1) * iy - iy * (iy + 1) / 2 + ix + (iz + 1) * (3 * (n + 1) * (n + 1) -
// 3 * (n + 1) * iz + 9 * (n + 1) + iz * iz - 4 * iz + 6) / 6;
// const uint I4 =
//((n + 1) - iz + 1) * (iy + 1) - (iy + 1) * (iy + 2) / 2 + ix + (iz + 1) * (3 *
//(n + 1) * (n + 1) - 3 * (n + 1) * iz + 9 * (n + 1) + iz * iz - 4 * iz + 6) /
// 6;
//
//						if (sum == n)
//							reference_connectivities.push_back({
// I1, I1 + 1, I2, I3 }); 						else if
// (sum
// ==
// n
// - 1) { auto sliced_hexahedral_connectivities =
// this->sliced_hexahedral_connectivities({ I1, I1 + 1, I2, I2 + 1, I3, I3 + 1,
// I4 }); 							for (auto&
// connectivity : sliced_hexahedral_connectivities)
//								reference_connectivities.push_back(std::move(connectivity));
//						}
//						else {
//							auto
// hexahedral_connectivities = this->hexahedral_connectivities({ I1, I1 + 1, I2,
// I2 + 1, I3, I3 + 1, I4, I4
// +
// 1 }); 							for (auto&
// connectivity : hexahedral_connectivities)
//								reference_connectivities.push_back(std::move(connectivity));
//						}
//					}
//				}
//			}
//
//			return reference_connectivities;
//		}
//		case Figure::hexahedral: {
//			const auto num_simplex = 6 * (n + 1) * (n + 1) * (n +
// 1); 			reference_connectivities.reserve(num_simplex);
//
//			for (int iz = 0; iz <= n; iz++) {
//				for (int iy = 0; iy <= n; iy++) {
//					for (int ix = 0; ix <= n; ix++) {
//						//	     I3��������I4
//						//      /��    /��
//						//     I3+1����I4+1
//						//     �� I1��������I2
//						//	   ��/    ��/
//						//     I1+1����I2+1
//						const uint I1 = ((n + 1) + 1) *
//((n
//+ 1)
//+
// 1) * iz + ((n + 1) + 1) * iy + ix; const uint I2 = ((n + 1) + 1) * ((n + 1) +
// 1) * iz + ((n + 1) + 1) * (iy + 1) + ix; const uint I3 = ((n + 1) + 1) * ((n
// + 1) + 1) * (iz + 1) + ((n + 1) + 1) * iy
// + ix; 						const uint I4 = ((n + 1)
// + 1)
// *
//((n + 1) + 1) * (iz + 1) + ((n + 1) + 1) * (iy + 1) + ix;
//
//						auto hexahedral_connectivities =
// this->hexahedral_connectivities({ I1, I1 + 1, I2, I2 + 1, I3, I3 + 1, I4, I4
// +
// 1 }); 						for (auto& connectivity
// : hexahedral_connectivities)
//							reference_connectivities.push_back(std::move(connectivity));
//					}
//				}
//			}
//
//			return reference_connectivities;
//		}
//}

// int ReferenceGeometry<space_dimension>::scale_function_order(void) const {
//		case Figure::tetrahedral:	return 0;
//		case Figure::hexahedral:	return 2;
// }

// std::vector<std::vector<uint>>
// ReferenceGeometry<space_dimension>::sliced_hexahedral_connectivities(const
// std::array<uint, 7>& node_indexes) const {
//	//	     4����������6
//	//      /��     ��
//  	//     5 ��   X ��
//	//     �� 0����������2
//	//	   ��/     /
//	//     1����������3
//	return { { node_indexes[0], node_indexes[1], node_indexes[2],
// node_indexes[4] },
//       		 { node_indexes[1], node_indexes[2], node_indexes[4],
//       node_indexes[5] }, 		 { node_indexes[2], node_indexes[4],
//       node_indexes[5], node_indexes[6] }, 		 { node_indexes[1],
//       node_indexes[3], node_indexes[5], node_indexes[2] }, 		 {
//       node_indexes[3], node_indexes[2], node_indexes[5], node_indexes[6] } };
// }

// std::vector<std::vector<uint>>
// ReferenceGeometry<space_dimension>::hexahedral_connectivities(const
// std::array<uint, 8>& node_indexes) const {
//	//	     4����������6
//	//      /��    /��
//	//     5����������7 ��
//	//     �� 0����������2
//	//	   ��/    ��/
//	//     1����������3
//	return { { node_indexes[0], node_indexes[1], node_indexes[2],
// node_indexes[4] }, 			 { node_indexes[1], node_indexes[2],
// node_indexes[4], node_indexes[5] }, 			 { node_indexes[2],
// node_indexes[4], node_indexes[5],
// node_indexes[6] }, 			 { node_indexes[1], node_indexes[3],
// node_indexes[5], node_indexes[2] }, 			 { node_indexes[3],
// node_indexes[2], node_indexes[5],
// node_indexes[6] }, 			 { node_indexes[3], node_indexes[6],
// node_indexes[5], node_indexes[7] } };
// }
//