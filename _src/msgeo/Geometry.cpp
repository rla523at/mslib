#include "msgeo/Geometry.h"

#include "msexception/Exception.h"
#include "msgeo/Reference_Geometry_Container.h"
#include "msmath/Vector.h"

namespace ms::geo
{

Geometry::Geometry(const Figure figure, const std::vector<Node_View>& consisting_nodes)
    : _reference_geometry(Reference_Geometry_Container::get(figure)),
      _nodes(consisting_nodes)
{
  this->_parametric_functions = this->_reference_geometry.cal_parametric_functions(this->_nodes);
}

Geometry::Geometry(const Figure figure, std::vector<Node_View>&& consisting_nodes)
    : _reference_geometry(Reference_Geometry_Container::get(figure)),
      _nodes(std::move(consisting_nodes))
{
  this->_parametric_functions = this->_reference_geometry.cal_parametric_functions(this->_nodes);
}

void Geometry::change_nodes(std::vector<Node_View>&& new_nodes)
{
  REQUIRE(this->_reference_geometry.is_valid_num_points(static_cast<int>(new_nodes.size())), "The number of nodes is not valid for the current reference geometry.");

  this->_nodes                = std::move(new_nodes);
  this->_parametric_functions = this->_reference_geometry.cal_parametric_functions(this->_nodes);

  if (this->_is_scale_function_initialized)
  {
    this->_scale_function = this->_reference_geometry.cal_scale_function(this->_parametric_functions);
  }
  if (this->_is_normal_functions_initialized)
  {
    this->_normal_functions = this->_reference_geometry.cal_normal_functions(this->_parametric_functions);
  }
}

void Geometry::cal_normal(double* normal, const Node_View node) const
{
  if (!this->_is_normal_functions_initialized)
  {
    this->_normal_functions                = this->_reference_geometry.cal_normal_functions(this->_parametric_functions);
    this->_is_normal_functions_initialized = true;
  }

  const auto dimension = this->dimension();
  for (int i = 0; i < dimension; ++i)
  {
    normal[i] = this->_normal_functions[i](node.to_vector_view());
  }
}

void Geometry::cal_projected_volumes(double* projected_volumes) const
{
  // calculating approximately

  const auto ref_geo_dim = this->_reference_geometry.dimension();
  if (ref_geo_dim <= 2)
  {
    const auto dimension = this->dimension();
    const auto num_nodes = this->_nodes.size();

    std::vector<double> coordinate_mins(dimension, std::numeric_limits<double>::max());
    std::vector<double> coordinate_maxs(dimension, std::numeric_limits<double>::min());

    for (int i = 0; i < num_nodes; ++i)
    {
      const auto& node = this->_nodes[i];
      for (int j = 0; j < dimension; ++j)
      {
        const auto coordinate = node[j];

        if (coordinate < coordinate_mins[j]) coordinate_mins[j] = coordinate;

        if (coordinate_maxs[j] < coordinate) coordinate_maxs[j] = coordinate;
      }
    }

    for (int i = 0; i < dimension; ++i)
    {
      projected_volumes[i] = coordinate_maxs[i] - coordinate_mins[i];
    }
  }
  else if (ref_geo_dim == 3)
  {
    EXCEPTION("not supproted space dimension");
    // double yz_projected_volume = 0.0;
    // double xz_projected_volume = 0.0;
    // double xy_projected_volume = 0.0;

    // const auto face_geometries = this->face_geometries();
    // for (const auto& geometry : face_geometries)
    //{
    //   const auto normal_v = geometry.normalized_normal_vector(geometry.center_point());

    //  Euclidean_Vector yz_plane_normalized_normal_vector = {1, 0, 0};
    //  Euclidean_Vector xz_plane_normalized_normal_vector = {0, 1, 0};
    //  Euclidean_Vector xy_plane_normalized_normal_vector = {0, 0, 1};

    //  const auto volume = geometry.volume();

    //  yz_projected_volume += volume * std::abs(normal_v.inner_product(yz_plane_normalized_normal_vector));
    //  xz_projected_volume += volume * std::abs(normal_v.inner_product(xz_plane_normalized_normal_vector));
    //  xy_projected_volume += volume * std::abs(normal_v.inner_product(xy_plane_normalized_normal_vector));
    //}

    // return {0.5 * yz_projected_volume, 0.5 * xz_projected_volume, 0.5 * xy_projected_volume};
  }
  else
  {
    EXCEPTION("not supproted space dimension");
  }
}

double Geometry::cal_volume(const int expected_scale_function_order) const
{
  const auto& quadrature_rule = this->get_quadrature_rule(expected_scale_function_order);

  auto volume = 0.0;
  for (const auto weight : quadrature_rule.weights)
  {
    volume += weight;
  }

  return volume;
}

Node Geometry::center(void) const
{
  const auto ref_center = this->_reference_geometry.center_point();

  const auto          dimension = this->dimension();
  std::vector<double> coordinates(dimension);

  for (int i = 0; i < dimension; ++i)
  {
    coordinates[i] = this->_parametric_functions[i](ref_center.to_vector_view());
  }

  return coordinates;
}

void Geometry::center(double* coordinates) const
{
  const auto ref_center = this->_reference_geometry.center_point();

  const auto dimension = this->dimension();
  for (int i = 0; i < dimension; ++i)
  {
    coordinates[i] = this->_parametric_functions[i](ref_center.to_vector_view());
  }
}

int Geometry::dimension(void) const
{
  return static_cast<int>(this->_parametric_functions.size());
}

Figure Geometry::face_figure(const int face_index) const
{
  return this->_reference_geometry.face_figure(face_index);
}

const std::vector<std::vector<int>>& Geometry::get_face_vnode_indexes_s(void) const
{
  return this->_reference_geometry.get_face_vnode_indexes_s();
}

std::vector<int> Geometry::face_node_indexes(const int face_index) const
{
  const auto num_nodes       = static_cast<int>(this->_nodes.size());
  const auto parameter_order = this->_reference_geometry.cal_parameter_order(num_nodes);

  const auto  face_figure  = this->_reference_geometry.face_figure(face_index);
  const auto& face_ref_geo = Reference_Geometry_Container::get(face_figure);

  return face_ref_geo.node_indexes(parameter_order);
}

const Quadrature_Rule& Geometry::get_quadrature_rule(const int integrand_degree) const
{
  REQUIRE(0 <= integrand_degree, "integrand degree should be positive");

  if (!this->_degree_to_quadrature_rule.contains(integrand_degree))
  {
    this->create_and_store_quadrature_rule(integrand_degree);
  }

  return this->_degree_to_quadrature_rule.at(integrand_degree);
}

bool Geometry::is_point(void) const
{
  return this->_reference_geometry.is_point();
}

bool Geometry::is_line(void) const
{
  return this->_reference_geometry.is_line();
}

Geometry_Consisting_Nodes_Info Geometry::make_partitioned_geometry_node_info(const int partition_order) const
{
  REQUIRE(0 <= partition_order, "partition order should not be negative");

  const auto& ref_pg_nodes_info = this->_reference_geometry.get_partition_geometry_nodes_info(partition_order);
  const auto& ref_nodes         = ref_pg_nodes_info.nodes;

  auto  pg_nodes_info = ref_pg_nodes_info;
  auto& nodes         = pg_nodes_info.nodes;

  const auto num_nodes = ref_nodes.num_nodes();

  for (int i = 0; i < num_nodes; ++i)
  {
    auto       node_wrap            = nodes[i];
    const auto ref_node_view        = ref_nodes[i];
    auto       node_vector_wrap     = node_wrap.to_vector_wrap();
    const auto ref_node_vector_view = ref_node_view.to_vector_view();

    this->_parametric_functions.calculate(node_vector_wrap, ref_node_vector_view);
  }

  // const auto& ref_numbered_nodes = ref_pg_nodes_info.numbered_nodes;

  // const auto num_new_nodes = ref_numbered_nodes.size();

  // auto partition_geometry_nodes_info = ref_pg_nodes_info;

  // for (int i = 0; i < num_new_nodes; ++i)
  //{
  //   auto&       node          = partition_geometry_nodes_info.numbered_nodes[i].node;
  //   const auto& ref_node      = ref_numbered_nodes[i].node;
  //   auto        new_node_wrap = node.to_vector_wrap();
  //   const auto  ref_node_view = ref_node.to_vector_view();

  //  this->_parametric_functions.calculate(new_node_wrap, ref_node_view);
  //}

  return pg_nodes_info;
}

int Geometry::num_faces(void) const
{
  return this->_reference_geometry.num_faces();
}

int Geometry::num_nodes(void) const
{
  return static_cast<int>(this->_nodes.size());
}

int Geometry::num_vertices(void) const
{
  return this->_reference_geometry.num_vertices();
}

// std::vector<Euclidean_Vector> Geometry::post_points(const ushort post_order) const
//{
//   const auto& ref_post_points = this->reference_geometry_->get_post_points(post_order);
//   const auto  num_post_points = ref_post_points.size();
//
//   std::vector<Euclidean_Vector> post_points;
//   post_points.reserve(num_post_points);
//
//   for (const auto& point : ref_post_points)
//   {
//     post_points.push_back(this->mapping_vf_(point));
//   }
//
//   return post_points;
// }
//
// std::vector<Euclidean_Vector> Geometry::post_element_centers(const ushort post_order) const
//{
//   const auto  post_points        = this->post_points(post_order);
//   const auto& ref_connectivities = this->reference_geometry_->get_connectivities(post_order);
//
//   const auto                    num_post_element = ref_connectivities.size();
//   std::vector<Euclidean_Vector> post_element_centers;
//   post_element_centers.reserve(num_post_element);
//
//   for (const auto& connectivity : ref_connectivities)
//   {
//     Euclidean_Vector center(this->space_dimension_);
//
//     for (const auto index : connectivity)
//     {
//       center += post_points[index];
//     }
//
//     center *= 1.0 / connectivity.size();
//
//     post_element_centers.push_back(center);
//   }
//
//   return post_element_centers;
// }
//
// std::vector<std::vector<int>> Geometry::post_connectivities(const ushort post_order, const size_t connectivity_start_index) const
//{
//   const auto& ref_connectivities = this->reference_geometry_->get_connectivities(post_order);
//
//   const auto                    num_connectivity = ref_connectivities.size();
//   std::vector<std::vector<int>> connectivities(num_connectivity);
//
//   for (ushort i = 0; i < num_connectivity; ++i)
//   {
//     auto&       connectivity      = connectivities[i];
//     const auto& ref_connecitivity = ref_connectivities[i];
//
//     const auto num_index = ref_connecitivity.size();
//     connectivity.resize(num_index);
//
//     for (ushort j = 0; j < num_index; ++j)
//     {
//       const auto new_index = static_cast<int>(ref_connecitivity[j] + connectivity_start_index);
//       connectivity[j]      = new_index;
//     }
//   }
//
//   return connectivities;
// }

void Geometry::create_and_store_quadrature_rule(const int integrand_degree) const
{
  if (!this->_is_scale_function_initialized)
  {
    this->_scale_function                = this->_reference_geometry.cal_scale_function(this->_parametric_functions);
    this->_is_scale_function_initialized = true;
  }

  REQUIRE(0 <= integrand_degree, "integrand degree should be positive");

  const auto  ref_quadrature_points  = this->_reference_geometry.quadrature_points(integrand_degree);
  const auto& ref_quadrature_weights = this->_reference_geometry.get_quadrature_weights(integrand_degree);

  const auto num_QP    = ref_quadrature_points.num_nodes();
  const auto dimension = this->dimension();

  std::vector<double> transformed_coordiantes(num_QP * dimension);
  std::vector<double> transformed_QW(num_QP);

  auto ptr = transformed_coordiantes.data();
  for (int i = 0; i < num_QP; ++i)
  {
    const auto ref_QP     = ref_quadrature_points[i];
    const auto ref_weight = ref_quadrature_weights[i];

    for (int j = 0; j < dimension; ++j)
    {
      ptr[j] = this->_parametric_functions[j](ref_QP.to_vector_view());
    }
    ptr += dimension;

    transformed_QW[i] = this->_scale_function(ref_QP.to_vector_view()) * ref_weight;
  }

  auto points          = Nodes(num_QP, dimension, std::move(transformed_coordiantes));
  auto quadrature_rule = Quadrature_Rule(std::move(points), std::move(transformed_QW));
  this->_degree_to_quadrature_rule.emplace(integrand_degree, std::move(quadrature_rule));
}

} // namespace ms::geo

// #include "Exception.h"
// #include "Figure.h"
// #include "Point_Const_Wrapper.h"
// #include "Reference_Geometry_Container.h"
//
// Parametric_Domain::Parametric_Domain(const Figure figure, const Points_CW points)
//     : reference_geometry_(RGeo_Container::get(figure))
//{
//     this->coordiantes_.resize(points.num_coordinates());
//     points.copy_coordinates(this->coordiantes_.data());
//
//     this->parametric_functions_ = this->reference_geometry_.cal_parametric_functions(points);
//     this->normal_functions_ = this->reference_geometry_.cal_normal_functions(this->parametric_functions_);
// }
//
// void Parametric_Domain::change_points(const Points_CW points)
//{
//     this->coordiantes_.resize(points.num_coordinates());
//     this->coordiantes_.shrink_to_fit();
//     points.copy_coordinates(this->coordiantes_.data());
//
//     this->parametric_functions_ = this->reference_geometry_.cal_parametric_functions(points);
// }
//
// Geometry::Geometry(const std::shared_ptr<const Reference_Geometry>&
//                        reference_goemetry,
//     std::vector<Euclidean_Vector>&& consisting_points)
//     : reference_geometry_(reference_goemetry)
//     , points_(std::move(consisting_points))
//{
//     this->space_dimension_ = this->reference_geometry_->check_space_dimension(this->points_);
//     this->mapping_vf_ = this->reference_geometry_->make_mapping_function(this->points_);
//     this->normal_vf_ = this->reference_geometry_->make_normal_vector_function(this->mapping_vf_);
//     this->scale_f_ = this->reference_geometry_->scale_function(this->mapping_vf_);
// };
//
// void Geometry::change_points(std::vector<Euclidean_Vector>&& points)
//{
//         REQUIRE(this->space_dimension_ ==
//  this->reference_geometry_->check_space_dimension(points), "space dimension
//  should be matched");
//
//	this->points_ = std::move(points);
//	this->mapping_vf_ =
//  this->reference_geometry_->make_mapping_function(this->points_);
//	this->normal_vf_ =
//  this->reference_geometry_->make_normal_vector_function(this->mapping_vf_);
//	this->scale_f_ =
//  this->reference_geometry_->scale_function(this->mapping_vf_);
// }
//
// bool Geometry::operator==(const Geometry& other) const
//{
//         return this->points_ == other.points_ && this->reference_geometry_ == other.reference_geometry_;
// }
//
// std::vector<Geometry> Geometry::face_geometries(void) const
//{
//         auto set_of_face_points = this->set_of_face_points();
//         auto face_reference_geometries = this->reference_geometry_->face_reference_geometries();
//         const auto
//             num_face
//             = face_reference_geometries.size();
//
//         std::vector<Geometry> face_geometries;
//         face_geometries.reserve(num_face);
//
//         for (ushort i = 0; i < num_face; ++i) {
//             face_geometries.push_back({ face_reference_geometries[i],
//                 std::move(set_of_face_points[i]) });
//         }
//
//         return face_geometries;
// }
//
// ushort Geometry::num_post_nodes(const ushort post_order) const
//{
//         return this->reference_geometry_->num_post_nodes(post_order);
// }
//
// ushort Geometry::num_post_elements(const ushort post_order) const
//{
//         return this->reference_geometry_->num_post_elements(post_order);
// }
//
// ushort Geometry::num_vertices(void) const
//{
//         return this->reference_geometry_->num_vertices();
// }
//
// Euclidean_Vector Geometry::normalized_normal_vector(const Euclidean_Vector&
//         point) const
//{
//         auto normal_v = this->normal_vf_(point);
//         normal_v.normalize();
//         return normal_v;
// }
//
// std::vector<Euclidean_Vector> Geometry::normalized_normal_vectors(const std::vector<Euclidean_Vector>& points) const
//{
//         const auto num_points = points.size();
//
//         std::vector<Euclidean_Vector> normalized_normal_vectors(num_points);
//         for (ushort i = 0; i < num_points; ++i) {
//             normalized_normal_vectors[i] = this->normalized_normal_vector(points[i]);
//         }
//
//         return normalized_normal_vectors;
// }
//

//

// Vector_Function<Polynomial> Geometry::orthonormal_basis_vector_function(const ushort solution_order) const
//{
//         const auto initial_basis_vector_function = this->initial_basis_vector_function(solution_order);
//         return ms::Gram_Schmidt_process(initial_basis_vector_function, *this);
// }
//
// std::vector<double> Geometry::projected_volumes(void) const
//{
//         // This only work for linear mesh and convex geometry and cell
//
//         if (this->space_dimension_ == 2) {
//             double x_projected_volume = 0.0;
//             double y_projected_volume = 0.0;
//
//             const auto set_of_face_points = this->set_of_face_points();
//             for (const auto& points : set_of_face_points) {
//                 const auto& start_point_v = points[0];
//                 const auto& end_point_v = points[1];
//                 const auto distance_v = end_point_v - start_point_v;
//
//                 x_projected_volume += std::abs(distance_v.at(0));
//                 y_projected_volume += std::abs(distance_v.at(1));
//             }
//
//             return { 0.5 * y_projected_volume, 0.5 * x_projected_volume };
//         } else if (this->space_dimension_ == 3) {
//             double yz_projected_volume = 0.0;
//             double xz_projected_volume = 0.0;
//             double xy_projected_volume = 0.0;
//
//             const auto face_geometries = this->face_geometries();
//             for (const auto& geometry : face_geometries) {
//                 const auto normal_v = geometry.normalized_normal_vector(geometry.center_point());
//
//                 Euclidean_Vector yz_plane_normalized_normal_vector = {
//                     1, 0, 0
//                 };
//                 Euclidean_Vector xz_plane_normalized_normal_vector = {
//                     0, 1, 0
//                 };
//                 Euclidean_Vector
//                     xy_plane_normalized_normal_vector
//                     = { 0, 0, 1 };
//
//                 const auto volume = geometry.volume();
//
//                 yz_projected_volume += volume * std::abs(normal_v.inner_product(yz_plane_normalized_normal_vector));
//                 xz_projected_volume += volume * std::abs(normal_v.inner_product(xz_plane_normalized_normal_vector));
//                 xy_projected_volume += volume * std::abs(normal_v.inner_product(xy_plane_normalized_normal_vector));
//             }
//
//             return { 0.5 * yz_projected_volume, 0.5 * xz_projected_volume,
//                 0.5
//                     * xy_projected_volume };
//         } else {
//             EXCEPTION("not supproted space dimension");
//             return {};
//         }
// }
//
// std::vector<std::vector<Euclidean_Vector>> Geometry::set_of_face_points(void)
//     const
//{
//         const auto set_of_face_point_index_sequences = this->reference_geometry_->set_of_face_node_index_sequences();
//         const auto num_face = set_of_face_point_index_sequences.size();
//
//         std::vector<std::vector<Euclidean_Vector>> set_of_face_points;
//         set_of_face_points.reserve(num_face);
//
//         for (size_t i = 0; i < num_face; ++i) {
//             set_of_face_points.push_back(ms::extract_by_index(this->points_,
//                 set_of_face_point_index_sequences[i]));
//         }
//
//         return set_of_face_points;
// }
//
// std::vector<std::vector<Euclidean_Vector>>
// Geometry::set_of_sub_simplex_vertices(void) const
//{
//         const auto set_of_sub_simplex_vertex_index_sequences = this->reference_geometry_->set_of_sub_simplex_vertex_index_sequences();
//         const auto num_sub_simplex = set_of_sub_simplex_vertex_index_sequences.size();
//
//         std::vector<std::vector<Euclidean_Vector>> set_of_simplex_vertices;
//         set_of_simplex_vertices.reserve(num_sub_simplex);
//
//         for (size_t i = 0; i < num_sub_simplex; ++i) {
//             set_of_simplex_vertices.push_back(ms::extract_by_index(this->points_,
//                 set_of_sub_simplex_vertex_index_sequences[i]));
//         }
//
//         return set_of_simplex_vertices;
// }
//
// std::vector<Geometry> Geometry::sub_simplex_geometries(void) const
//{
//         auto sub_simplex_reference_geometries = this->reference_geometry_->sub_simplex_reference_geometries();
//         auto set_of_sub_simplex_vertices = this->set_of_sub_simplex_vertices();
//
//         const auto num_sub_simplex = sub_simplex_reference_geometries.size();
//
//         std::vector<Geometry> sub_simplex_geometries;
//         sub_simplex_geometries.reserve(num_sub_simplex);
//
//         for (ushort i = 0; i < num_sub_simplex; ++i) {
//             sub_simplex_geometries.push_back({ sub_simplex_reference_geometries[i], std::move(set_of_sub_simplex_vertices[i]) });
//         }
//
//         return sub_simplex_geometries;
// }
//
// double Geometry::volume(void) const
//{
//         const auto& quadrature_rule = this->get_quadrature_rule(0);
//
//         auto volume = 0.0;
//         for (const auto weight : quadrature_rule.weights) {
//             volume += weight;
//         }
//
//         return volume;
// }
//
// std::vector<Euclidean_Vector> Geometry::vertices(void) const
//{
//         const auto num_vertices = this->reference_geometry_->num_vertices();
//         return { this->points_.begin(), this->points_.begin() + num_vertices };
// }
//
// bool Geometry::is_line(void) const
//{
//         return this->reference_geometry_->is_line();
// }
//
// bool Geometry::is_simplex(void) const
//{
//         return this->reference_geometry_->is_simplex();
// }
//
// const Quadrature_Rule& Geometry::get_quadrature_rule(const ushort
//         integrand_degree) const
//{
//         if (!this->degree_to_quadrature_rule_.contains(integrand_degree))
//             this->degree_to_quadrature_rule_.emplace(integrand_degree,
//                 this->make_quadrature_rule(integrand_degree));
//
//         return this->degree_to_quadrature_rule_.at(integrand_degree);
// }
//
// Vector_Function<Polynomial> Geometry::initial_basis_vector_function(const ushort solution_order) const
//{
//         const auto num_basis = ms::combination_with_repetition(1 + this->space_dimension_, solution_order);
//
//         std::vector<Polynomial> initial_basis_functions(num_basis);
//
//         ushort index = 0;
//         if (this->space_dimension_ == 2) {
//             Polynomial x("x0");
//             Polynomial y("x1");
//
//             const auto center_point = this->center_point();
//             const auto x_c = center_point.at(0);
//             const auto y_c = center_point.at(1);
//
//             // 1 (x - x_c) (y - y_c)  ...
//             for (ushort a = 0; a <= solution_order; ++a) {
//                 for (ushort b = 0; b <= a; ++b) {
//                     initial_basis_functions[index++] = ((x - x_c) ^ (a - b)) * ((y - y_c) ^ b);
//                 }
//             }
//         } else if (this->space_dimension_ == 3) {
//             Polynomial x("x0");
//             Polynomial y("x1");
//             Polynomial z("x2");
//
//             const auto center_point = this->center_point();
//             const auto x_c = center_point.at(0);
//             const auto y_c = center_point.at(1);
//             const auto z_c = center_point.at(2);
//
//             // 1 (x - x_c) (y - y_c) (z - z_c) ...
//             for (ushort a = 0; a <= solution_order; ++a) {
//                 for (ushort b = 0; b <= a; ++b) {
//                     for (ushort c = 0; c <= b; ++c) {
//                         initial_basis_functions[index++] = ((x - x_c) ^ (a - b)) * ((y - y_c) ^ (b - c)) * ((z - z_c) ^ c);
//                     }
//                 }
//             }
//         } else
//             EXCEPTION("not supported space dimension");
//
//         return initial_basis_functions;
// }
//
// bool Geometry::is_axis_parallel_node(const Euclidean_Vector& node) const
//{
//         for (const auto& my_node : this->points_) {
//             if (my_node.is_axis_translation(node)) {
//                 return true;
//             }
//         }
//         return false;
// }
//
// bool Geometry::is_on_axis_plane(const ushort axis_tag) const
//{
//         if (this->space_dimension_ <= axis_tag) {
//             return false;
//         }
//
//         const auto& ref_node = this->points_.front();
//
//         constexpr auto epsilon = 1.0E-10;
//         for (ushort j = 1; j < this->points_.size(); ++j) {
//             if (std::abs(this->points_[j][axis_tag] - ref_node[axis_tag]) > epsilon) {
//                 return false;
//             }
//         }
//         return true;
// }
//
// bool Geometry::is_on_this_axis_plane(const ushort axis_tag, const double reference_value) const
//{
//         if (this->space_dimension_ <= axis_tag) {
//             return false;
//         }
//
//         constexpr auto epsilon = 1.0E-10;
//         for (const auto& node : this->points_) {
//             if (std::abs(node[axis_tag] - reference_value) > epsilon) {
//                 return false;
//             }
//         }
//
//         return true;
// }
//
// bool Geometry::is_on_same_axis_plane(const Geometry& other) const
//{
//         for (ushort i = 0; i < this->space_dimension_; ++i) {
//             if (this->is_on_axis_plane(i)) {
//                 const auto reference_value = this->points_.front()[i];
//
//                 if (other.is_on_this_axis_plane(i, reference_value)) {
//                     return true;
//                 }
//             }
//         }
//
//         return false;
// }
//

// namespace ms {
//// double integrate(const Polynomial& integrand, const Quadrature_Rule&
// quadrature_rule)
//	//{
//	//	const auto& QP_set = quadrature_rule.points;
//	//	const auto& QW_set = quadrature_rule.weights;
//
//	//	double result = 0.0;
//	//	for (ushort i = 0; i < QP_set.size(); ++i)
//	//		result += integrand(QP_set[i]) * QW_set[i];
//
//	//	return result;
//	//}
//
//	//double integrate(const Polynomial& integrand, const Geometry&
// geometry)
//	//{
//	//	const auto quadrature_rule =
// geometry.get_quadrature_rule(integrand.degree());
// //	return ms::integrate(integrand, quadrature_rule);
// //}
//
// double inner_product(const Polynomial& f1, const Polynomial& f2, const Geometry& geometry)
// {
//        const auto integrand_degree = f1.degree() + f2.degree();
//
//        const auto quadrature_rule = geometry.get_quadrature_rule(integrand_degree);
//        const auto&
//            QP_set
//            = quadrature_rule.points;
//        const auto& QW_set = quadrature_rule.weights;
//
//        double result = 0.0;
//        for (ushort i = 0; i < QP_set.size(); ++i)
//            result += f1(QP_set[i]) * f2(QP_set[i]) * QW_set[i];
//
//        return result;
// }
//
// double L2_Norm(const Polynomial& function, const Geometry& geometry)
// {
//        return std::sqrt(ms::inner_product(function, function,
//            geometry));
// }
//
// Vector_Function<Polynomial> Gram_Schmidt_process(const Vector_Function<Polynomial>& vf, const Geometry& geometry)
// {
//        const auto range_dimension = vf.size();
//
//        std::vector<Polynomial> normalized_functions(range_dimension);
//
//        for (ushort i = 0; i < range_dimension; ++i) {
//            normalized_functions[i] = vf[i];
//
//            for (ushort j = 0; j < i; ++j)
//                normalized_functions[i] -= ms::inner_product(normalized_functions[i], normalized_functions[j], geometry)
//                    * normalized_functions[j];
//
//            normalized_functions[i] *= 1.0 / ms::L2_Norm(normalized_functions[i], geometry);
//        }
//
//        return normalized_functions;
// }
// };