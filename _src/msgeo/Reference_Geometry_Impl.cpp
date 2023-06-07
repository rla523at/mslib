#include "Reference_Geometry_Impl.h"

#include "Figure.h"
#include "Node.h"

#include "msexception/Exception.h"
#include "msmath/Matrix.h"

namespace ms::geo
{

ms::sym::Polynomials Reference_Point::cal_normal_functions(const ms::sym::Polynomials& parametric_functions) const
{
  ms::sym::Polynomials result;
  result[0] = 1.0;
  return result;
}

ms::sym::Polynomials Reference_Point::cal_parametric_functions(const std::vector<Node_Const_Wrapper>& consisting_node_wraps) const
{
  REQUIRE(consisting_node_wraps.size() == 1, "Point must consist of only one single point.");

  const auto& node_cwrap = consisting_node_wraps.front();
  const auto  dimension = node_cwrap.dimension();

  ms::sym::Polynomials result(dimension);
  for (int i = 0; i < dimension; ++i)
  {
    result[i] = node_cwrap[i];
  }

  return result;
}

int Reference_Point::cal_parameter_order(const int num_points) const
{
  return 0;
}

ms::sym::Symbol Reference_Point::cal_scale_function(const ms::sym::Polynomials& parametric_functions) const
{
  EXCEPTION("Point doesn't have scale function");
  return {};
}

Node_Const_Wrapper Reference_Point::center_point(void) const
{
  return {1, this->_center_coords.data()};
}

int Reference_Point::dimension(void) const
{
  return 0;
}

Figure Reference_Point::face_figure(const int face_index) const
{
  EXCEPTION("Point doesn't have any face");
  return Figure::NOT_FIGURE;
}

std::vector<std::vector<int>> Reference_Point::face_index_to_face_vnode_indexes(void) const
{
  EXCEPTION("Point doesn't have any faces");
  return {};
}

bool Reference_Point::is_valid_num_points(const int num_points) const
{
  return num_points == 1;
}

bool Reference_Point::is_point(void) const
{
  return true;
}

bool Reference_Point::is_line(void) const
{
  return false;
}

std::vector<int> Reference_Point::node_indexes(const int parameter_order) const
{
  return {0};
}

int Reference_Point::num_faces(void) const
{
  EXCEPTION("Point doesn't have any vertices");
  return -1;
}

int Reference_Point::num_vertices(void) const
{
  return 1;
}

Nodes_Const_Wrapper Reference_Point::quadrature_points(const int integrand_degree) const
{
  EXCEPTION("Integration is not possible over a point.");
  return {Coordinates_Type::NOT_SUPPROTED, -1, -1, nullptr};
}

const std::vector<double>& Reference_Point::get_quadrature_weights(const int integrand_degree) const
{
  EXCEPTION("Integration is not possible over a point.");
  return {};
}

const Partition_Data& Reference_Point::get_partition_data(const int partition_order) const
{
  EXCEPTION("Partition is not possible over a point.");
  return {Nodes{Coordinates_Type::NOT_SUPPROTED, -1, -1}, {}};
}

ms::sym::Polynomials Reference_Geometry_Common::cal_parametric_functions(const std::vector<Node_Const_Wrapper>& consisting_nodes) const
{
  const auto  num_nodes       = static_cast<int>(consisting_nodes.size());
  const auto  param_order     = this->cal_parameter_order(num_nodes);
  const auto& shape_functions = this->get_shape_functions(param_order);

  const auto           dim = consisting_nodes.front().dimension();
  ms::sym::Polynomials parametric_functions(dim);

  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < num_nodes; ++j)
    {
      const auto& node = consisting_nodes[j];
      parametric_functions[i] += node[i] * shape_functions[j];
    }
  }

  return parametric_functions;
}

bool Reference_Geometry_Common::is_point(void) const
{
  return false;
}

Nodes_Const_Wrapper Reference_Geometry_Common::quadrature_points(const int integrand_degree) const
{
  REQUIRE(0 <= integrand_degree, "integrand degree should be positive");

  const auto tag = this->cal_quadrature_rule_tag(integrand_degree);

  if (!this->_tag_to_quadrature_points.contains(tag))
  {
    this->create_and_store_quadrature_points(tag);
  }

  return this->_tag_to_quadrature_points.at(tag);
}

const std::vector<double>& Reference_Geometry_Common::get_quadrature_weights(const int integrand_degree) const
{
  const auto tag = this->cal_quadrature_rule_tag(integrand_degree);

  if (!this->_tag_to_quadrature_weights.contains(tag))
  {
    this->create_and_store_quadrature_weights(tag);
  }

  return this->_tag_to_quadrature_weights.at(tag);
}

const Partition_Data& Reference_Geometry_Common::get_partition_data(const int partition_order) const
{
  if (!this->_order_to_partition_data.contains(partition_order))
  {
    this->create_and_store_partition_data(partition_order);
  }

  return this->_order_to_partition_data.at(partition_order);
}

const ms::sym::Polynomials& Reference_Geometry_Common::get_shape_functions(const int porder) const
{
  if (!this->_parameter_order_to_shape_functions.contains(porder))
  {
    this->_parameter_order_to_shape_functions.emplace(porder, this->make_shape_functions(porder));
  }

  return this->_parameter_order_to_shape_functions[porder];
}

ms::sym::Polynomials Reference_Geometry_Common::make_shape_functions(const int parameter_order) const
{
  constexpr auto coordinate_type = Coordinates_Type::NODAL;

  const auto          num_ref_points = this->num_parametric_function_reference_points(parameter_order);
  const auto          dim            = this->dimension();
  const auto          coords         = this->make_parametric_functions_reference_coords(parameter_order);
  Nodes_Const_Wrapper ref_points(coordinate_type, num_ref_points, dim, coords.data());

  const auto bases     = this->make_parametric_function_bases(parameter_order);
  const auto num_bases = bases.size();
  REQUIRE(num_ref_points == num_bases, "it should be same for invertibility");

  const auto       n = static_cast<int>(num_bases);
  ms::math::Matrix coefficients(n, n);

  for (int i = 0; i < n; ++i)
  {
    const auto& basis = bases[i];

    for (int j = 0; j < n; ++j)
    {
      const auto& point = ref_points[j];

      coefficients.at(i, j) = basis(point);
    }
  }

  coefficients.inverse();

  ms::sym::Polynomials shape_functions(n);
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      shape_functions[i] += coefficients.at(i, j) * bases[j];
    }
  }

  return shape_functions;
}

Node_Const_Wrapper Reference_Line::center_point(void) const
{
  return Node_Const_Wrapper(this->dimension(), this->center_coords_.data());
};

ms::sym::Polynomials Reference_Line::cal_normal_functions(const ms::sym::Polynomials& curve) const
{
  return ms::geo::cal_paramteric_curve_normal_functions(curve);
}

int Reference_Line::cal_parameter_order(const int num_points) const
{
  switch (num_points)
  {
  case 2:
    return 0;
  case 3:
    return 1;
  default:
    EXCEPTION("given num nodes does not match any order of parameter");
    return -1;
    break;
  }
}

Figure Reference_Line::face_figure(const int face_index) const
{
  REQUIRE(face_index < this->num_faces(), "face index can't exceed num face");

  return Figure::POINT;
}

std::vector<std::vector<int>> Reference_Line::face_index_to_face_vnode_indexes(void) const
{
  // 0 �������� 1
  constexpr auto num_faces = 2;

  std::vector<std::vector<int>> result(num_faces);
  result[0] = {0};
  result[1] = {1};

  return result;
}

bool Reference_Line::is_valid_num_points(const int num_points) const
{
  switch (num_points)
  {
  case 2:
  case 3:
    return true;
  default:
    return false;
  }
}

bool Reference_Line::is_line(void) const
{
  return true;
}

std::vector<int> Reference_Line::node_indexes(const int parameter_order) const
{
  // 0 ������ 2 ������ ������ ������ order+1 ������ 1

  std::vector<int> node_indexes(parameter_order + 1);

  for (auto i = 0; i < parameter_order + 1; ++i)
  {
    node_indexes[i] = i;
  }

  return node_indexes;
}

int Reference_Line::num_faces(void) const
{
  return 2;
}

int Reference_Line::num_vertices(void) const
{
  return 2;
}

ms::sym::Symbol Reference_Line::cal_scale_function(const ms::sym::Polynomials& parametric_functions) const
{
  constexpr int r = 0;

  const auto df_dr = ms::sym::get_differentiate(parametric_functions, r);
  return ms::sym::cal_L2_norm(df_dr);
}

int Reference_Line::dimension(void) const
{
  return 1;
}

int Reference_Line::num_quadrature_points(const int tag) const
{
  return tag + 1;
}

int Reference_Line::num_parametric_function_reference_points(const int param_order) const
{
  return param_order + 2;
}

int Reference_Line::cal_quadrature_rule_tag(const int integrand_degree) const
{
  return integrand_degree / 2;
}

void Reference_Line::create_and_store_partition_data(const int partition_order) const
{
  constexpr auto dimension      = 1;
  constexpr auto X0_start_coord = -1.0;

  const auto num_consisting_nodes    = partition_order + 2;
  const auto num_consisting_elements = partition_order + 1;

  Nodes consisting_nodes(Coordinates_Type::NODAL, num_consisting_nodes, dimension);

  const auto delta = 2.0 / num_consisting_elements;

  for (int i = 0; i < num_consisting_nodes; ++i)
  {
    auto node_wrap = consisting_nodes[i];

    const auto X0_coord = X0_start_coord + delta * i;
    node_wrap[0]        = X0_coord;
  }

  std::vector<std::vector<int>> connectivities(num_consisting_elements);

  for (int i = 0; i < num_consisting_elements; i++)
  {
    auto& connectivity = connectivities[i];

    //   i �������� i+1
    connectivity = {i, i + 1};
  }

  Partition_Data data(std::move(consisting_nodes), std::move(connectivities));
  this->_order_to_partition_data.emplace(partition_order, std::move(data));
}

void Reference_Line::create_and_store_quadrature_points(const int tag) const
{
  REQUIRE(0 <= tag, "param_order can not be negative");

  std::vector<double> coordinates;

  switch (tag)
  {
  case 0:
    coordinates = {0.000000000000000};
    break;
  case 1:
    coordinates = {-0.577350269189626, 0.577350269189626};
    break;
  case 2:
    coordinates = {-0.774596669241483, 0.000000000000000, 0.774596669241483};
    break;
  case 3:
    coordinates = {-0.861136311594052, -0.339981043584856, 0.339981043584856, 0.861136311594052};
    break;
  case 4:
    coordinates = {-0.906179845938664, -0.538469310105683, 0.000000000000000, 0.538469310105683, 0.906179845938664};
    break;
  case 5:
    coordinates = {-0.932469514203152, -0.661209386466264, -0.238619186083197, 0.238619186083197, 0.661209386466264, 0.932469514203152};
    break;
  case 6:
    coordinates = {-0.949107912342758, -0.741531185599394, -0.405845151377397, 0.000000000000000, 0.405845151377397, 0.741531185599394, 0.949107912342758};
    break;
  case 7:
    coordinates = {-0.960289856497536, -0.796666477413627, -0.525532409916329, -0.183434642495650, 0.183434642495650, 0.525532409916329, 0.796666477413627, 0.960289856497536};
    break;
  case 8:
    coordinates = {-0.968160239507626, -0.836031107326636, -0.613371432700590, -0.324253423403809, 0.000000000000000, 0.324253423403809, 0.613371432700590, 0.836031107326636, 0.968160239507626};
    break;
  case 9:
    coordinates = {-0.973906528517172, -0.865063366688985, -0.679409568299024, -0.433395394129247, -0.148874338981631, 0.148874338981631, 0.433395394129247, 0.679409568299024, 0.865063366688985, 0.973906528517172};
    break;
  case 10:
    coordinates = {-0.978228658146057, -0.887062599768095, -0.730152005574049, -0.519096129206812, -0.269543155952345, 0.000000000000000, 0.269543155952345, 0.519096129206812, 0.730152005574049, 0.887062599768095, 0.978228658146057};
    break;

  default:
    EXCEPTION("unsupported param_order");
    break;
  }

  constexpr auto type      = Coordinates_Type::NODAL;
  const auto     num_nodes = tag + 1;
  constexpr auto dimension = 1;

  this->_tag_to_quadrature_points.emplace(tag, Nodes(type, num_nodes, dimension, std::move(coordinates)));
}

void Reference_Line::create_and_store_quadrature_weights(const int tag) const
{
  REQUIRE(0 <= tag, "param_order can not be negative");

  std::vector<double> weights;

  switch (tag)
  {
  case 0:
    weights = {2.000000000000000};
    break;
  case 1:
    weights = {1.000000000000000, 1.000000000000000};
    break;
  case 2:
    weights = {0.555555555555554, 0.888888888888889, 0.555555555555554};
    break;
  case 3:
    weights = {0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454};
    break;
  case 4:
    weights = {0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189};
    break;
  case 5:
    weights = {0.171324492379171, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.171324492379171};
    break;
  case 6:
    weights = {0.129484966168870, 0.279705391489277, 0.381830050505119, 0.417959183673469, 0.381830050505119, 0.279705391489277, 0.129484966168870};
    break;
  case 7:
    weights = {0.101228536290377, 0.222381034453374, 0.313706645877887, 0.362683783378362, 0.362683783378362, 0.313706645877887, 0.222381034453374, 0.101228536290377};
    break;
  case 8:
    weights = {0.081274388361575, 0.180648160694857, 0.260610696402936, 0.312347077040003, 0.330239355001260, 0.312347077040003, 0.260610696402936, 0.180648160694857, 0.081274388361575};
    break;
  case 9:
    weights = {0.066671344308688, 0.149451349150581, 0.219086362515982, 0.269266719309996, 0.295524224714753, 0.295524224714753, 0.269266719309996, 0.219086362515982, 0.149451349150581, 0.066671344308688};
    break;
  case 10:
    weights = {0.055668567116174, 0.125580369464904, 0.186290210927734, 0.233193764591990, 0.262804544510247, 0.272925086777901, 0.262804544510247, 0.233193764591990, 0.186290210927734, 0.125580369464904, 0.055668567116174};
    break;
  default:
    EXCEPTION("unsupported param_order");
    break;
  }

  this->_tag_to_quadrature_weights.emplace(tag, std::move(weights));
}

std::vector<double> Reference_Line::make_parametric_functions_reference_coords(const int tag) const
{
  REQUIRE(0 <= tag, "param_order can not be negative");

  std::vector<double> VCP_coords;

  switch (tag)
  {
  case 0:
    VCP_coords = {
        -1.0,
        1.0};
    break;
  case 1:
    VCP_coords = {
        -1.0,
        1.0,
        0.0};
    break;
  default:
    EXCEPTION("unsupported param_order");
    break;
  }

  return VCP_coords;
}

ms::sym::Polynomials Reference_Line::make_parametric_function_bases(const int tag) const
{
  REQUIRE(0 <= tag, "param_order can not be negative");

  const auto           num_bases = tag + 2;
  ms::sym::Polynomials VCF_bases(num_bases);

  ms::sym::Polynomial r("x0");
  for (int i = 0; i < num_bases; ++i)
  {
    VCF_bases[i] = (r ^ i);
  }

  // 1 r r^2 ...
  return VCF_bases;
}

} // namespace ms::geo

namespace ms::geo
{

ms::sym::Symbols cal_principal_normal(const ms::sym::Polynomials& curve)
{
  const auto var_index = 0;

  auto tangent      = ms::sym::get_differentiate(curve, var_index);
  auto unit_tangent = ms::sym::get_normalize(tangent);

  return ms::sym::get_differentiate(unit_tangent, var_index);
}

ms::sym::Symbol cal_curvature(const ms::sym::Polynomials& curve)
{
  const auto var_index = 0;

  auto       tangent      = ms::sym::get_differentiate(curve, var_index);
  auto       unit_tangent = ms::sym::get_normalize(tangent);
  const auto diff_ut      = ms::sym::get_differentiate(unit_tangent, var_index);

  auto curvature = ms::sym::cal_L2_norm(diff_ut) / ms::sym::cal_L2_norm(tangent);
  return curvature;
}

ms::sym::Polynomials cal_paramteric_curve_normal_functions(const ms::sym::Polynomials& parameteric_curve)
{
  const auto range_dimension = parameteric_curve.size();
  REQUIRE(2 <= range_dimension, "1D can't have normal");

  ms::sym::Polynomials normal_functions(range_dimension);

  const auto& x = parameteric_curve[0];
  const auto& y = parameteric_curve[1];

  constexpr auto var_index = 0;
  const auto     dx_dr     = x.get_diff_polynomial(var_index);
  const auto     dy_dr     = y.get_diff_polynomial(var_index);

  normal_functions[0] = -1 * dy_dr;
  normal_functions[1] = dx_dr;

  return normal_functions;
}

} // namespace ms::geo

// ms::sym::Polynomials Reference_Line::cal_normal_functions(
//     const ms::sym::Polynomials& parameteric_curve) const {
//   return {};
// }

// ms::sym::Polynomials cal_paramteric_curve_normal_functions(
//     const ms::sym::Polynomials& parameteric_curve) {
//   const auto codomain_dim = parameteric_curve.size();
//   ms::sym::Polynomials tangent(codomain_dim);
//
//   // variable index
//   constexpr auto r = 0;
//   for (int i = 0; i < codomain_dim; ++i) {
//     tangent[i] = parameteric_curve[i].get_differentiate(r);
//   }
//
//   Polynomial A;
//   for (int i = 0; i < codomain_dim; ++i) {
//     A += tangent[i] * tangent[i];
//   }
//
// }

//
// std::vector<std::shared_ptr<const Reference_Geometry>>
// Reference_Line::face_reference_geometries(void) const {
//  EXCEPTION("Reference Node class is not supported yet");
//  return {};
//}
//
// Figure Reference_Line::figure(void) const { return Figure::line; }
//
// ushort Reference_Line::num_vertices(void) const { return 2; }
//

//
// Quadrature_Rule Reference_Line::make_quadrature_rule(
//    const ushort integrand_order) const {
//  switch (integrand_order) {
//    case 0:
//    case 1:
//      return {{{0.000000000000000}}, {2.000000000000000}};
//    case 2:
//    case 3:
//      return {{{-0.577350269189626}, {0.577350269189626}},
//              {1.000000000000000, 1.000000000000000}};
//    case 4:
//    case 5:
//      return {{{-0.774596669241483}, {0.000000000000000},
//      {0.774596669241483}},
//              {0.555555555555554, 0.888888888888889, 0.555555555555554}};
//    case 6:
//    case 7:
//      return {{{-0.861136311594052},
//               {-0.339981043584856},
//               {0.339981043584856},
//               {0.861136311594052}},
//              {0.347854845137454, 0.652145154862546, 0.652145154862546,
//               0.347854845137454}};
//    case 8:
//    case 9:
//      return {{{-0.906179845938664},
//               {-0.538469310105683},
//               {0.000000000000000},
//               {0.538469310105683},
//               {0.906179845938664}},
//              {0.236926885056189, 0.478628670499366, 0.568888888888889,
//               0.478628670499366, 0.236926885056189}};
//    case 10:
//    case 11:
//      return {{{-0.932469514203152},
//               {-0.661209386466264},
//               {-0.238619186083197},
//               {0.238619186083197},
//               {0.661209386466264},
//               {0.932469514203152}},
//              {0.171324492379171, 0.360761573048139, 0.467913934572691,
//               0.467913934572691, 0.360761573048139, 0.171324492379171}};
//    case 12:
//    case 13:
//      return {{{-0.949107912342758},
//               {-0.741531185599394},
//               {-0.405845151377397},
//               {0.000000000000000},
//               {0.405845151377397},
//               {0.741531185599394},
//               {0.949107912342758}},
//              {0.129484966168870, 0.279705391489277, 0.381830050505119,
//               0.417959183673469, 0.381830050505119, 0.279705391489277,
//               0.129484966168870}};
//    case 14:
//    case 15:
//      return {{{-0.960289856497536},
//               {-0.796666477413627},
//               {-0.525532409916329},
//               {-0.183434642495650},
//               {0.183434642495650},
//               {0.525532409916329},
//               {0.796666477413627},
//               {0.960289856497536}},
//              {0.101228536290377, 0.222381034453374, 0.313706645877887,
//               0.362683783378362, 0.362683783378362, 0.313706645877887,
//               0.222381034453374, 0.101228536290377}};
//    case 16:
//    case 17:
//      return {{{-0.968160239507626},
//               {-0.836031107326636},
//               {-0.613371432700590},
//               {-0.324253423403809},
//               {0.000000000000000},
//               {0.324253423403809},
//               {0.613371432700590},
//               {0.836031107326636},
//               {0.968160239507626}},
//              {0.081274388361575, 0.180648160694857, 0.260610696402936,
//               0.312347077040003, 0.330239355001260, 0.312347077040003,
//               0.260610696402936, 0.180648160694857, 0.081274388361575}};
//    case 18:
//    case 19:
//      return {{{-0.973906528517172},
//               {-0.865063366688985},
//               {-0.679409568299024},
//               {-0.433395394129247},
//               {-0.148874338981631},
//               {0.148874338981631},
//               {0.433395394129247},
//               {0.679409568299024},
//               {0.865063366688985},
//               {0.973906528517172}},
//              {0.066671344308688, 0.149451349150581, 0.219086362515982,
//               0.269266719309996, 0.295524224714753, 0.295524224714753,
//               0.269266719309996, 0.219086362515982, 0.149451349150581,
//               0.066671344308688}};
//    case 20:
//    case 21:
//      return {{{-0.978228658146057},
//               {-0.887062599768095},
//               {-0.730152005574049},
//               {-0.519096129206812},
//               {-0.269543155952345},
//               {0.000000000000000},
//               {0.269543155952345},
//               {0.519096129206812},
//               {0.730152005574049},
//               {0.887062599768095},
//               {0.978228658146057}},
//              {0.055668567116174, 0.125580369464904, 0.186290210927734,
//               0.233193764591990, 0.262804544510247, 0.272925086777901,
//               0.262804544510247, 0.233193764591990, 0.186290210927734,
//               0.125580369464904, 0.055668567116174}};
//    default:
//      EXCEPTION("not supported integrand order");
//      return {};
//  }
//};
//
// bool Reference_Line::is_simplex(void) const { return false; };
//
// Vector_Function<Polynomial> Reference_Line::make_normal_vector_function(
//    const Vector_Function<Polynomial>& mapping_function) const {
//  constexpr auto r = 0;
//  const auto T_r = mapping_function.get_differentiate(r);
//
//  constexpr auto dimension = 2;
//  ms::sym::Polynomials normal_vector_function(2);
//  normal_vector_function[0] = -1 * T_r.at(1);
//  normal_vector_function[1] = T_r.at(0);
//
//  return normal_vector_function;
//}
//
// Euclidean_Vector Reference_Line::random_point(void) const {
//  std::random_device rd;
//  std::mt19937 gen(rd());
//  std::uniform_real_distribution<double> distribution(-1.0, 1.0);
//
//  return {distribution(gen)};
//}
//
// bool Reference_Line::is_line(void) const { return true; }
//
// std::vector<std::vector<ushort>>
// Reference_Line::set_of_face_vertex_index_sequences(void) const {
//  // 0 �������� 1
//  const std::vector<ushort> face0_node_index = {0};
//  const std::vector<ushort> face1_node_index = {1};
//  return {face0_node_index, face1_node_index};
//};
//
// std::vector<std::vector<ushort>>
// Reference_Line::set_of_face_node_index_sequences(void) const {
//  // 0 �������� 1
//  const std::vector<ushort> face0_node_index = {0};
//  const std::vector<ushort> face1_node_index = {1};
//  return {face0_node_index, face1_node_index};
//};
//
// std::vector<std::shared_ptr<const Reference_Geometry>>
// Reference_Line::sub_simplex_reference_geometries(void) const {
//  EXCEPTION("This routine is not supported for Line element");
//  return {};
//}
//
// std::vector<std::vector<ushort>>
// Reference_Line::set_of_sub_simplex_vertex_index_sequences(void) const {
//  EXCEPTION("This routine is not supported for Line element");
//  return {};
//}
//

// std::vector<Euclidean_Vector> Reference_Line::make_mapping_points(void)
// const
// {
//  // 0 �������� 1
//  switch (this->order_) {
//    case 1:
//      return {{-1, 0}, {1, 0}};
//    case 2:
//      return {{-1, 0}, {1, 0}, {0, 0}};
//    default:
//      EXCEPTION("unsuported figure order");
//      return {};
//  }
//}
//

// Vector_Function<Polynomial>
// Reference_Line::make_mapping_monomial_vector_function(void) const {
//  const auto num_monomial = this->order_ + 1;
//  ms::sym::Polynomials mapping_monomial_vector(num_monomial);
//
//  Polynomial r("x0");
//  for (ushort a = 0, index = 0; a <= this->order_; ++a)
//    mapping_monomial_vector[index++] = (r ^ a);
//
//  return mapping_monomial_vector;  // 1 r r^2 ...
//}
//
// Reference_Triangle::Reference_Triangle(const ushort order) {
//  this->order_ = order;
//  this->initialize();
//}
//
// Euclidean_Vector Reference_Triangle::center_point(void) const {
//  return {-1.0 / 3.0, -1.0 / 3.0};
//};
//
// std::vector<std::shared_ptr<const Reference_Geometry>>
// Reference_Triangle::face_reference_geometries(void) const {
//  constexpr auto num_face = 3;
//  std::vector<std::shared_ptr<const Reference_Geometry>>
//      face_reference_geometries(num_face);
//
//  for (ushort i = 0; i < num_face; ++i) {
//    face_reference_geometries[i] =
//    Reference_Geometry_Container::get_shared_ptr(
//        Figure::line, this->order_);
//  }
//
//  return face_reference_geometries;
//}
//
// Figure Reference_Triangle::figure(void) const { return Figure::triangle; }
//
// ushort Reference_Triangle::num_vertices(void) const { return 3; }
//
// ushort Reference_Triangle::num_post_nodes(const ushort post_order) const {
//  const auto n = post_order;
//  return static_cast<ushort>((n + 2) * (n + 3) * 0.5);
//}
//
// ushort Reference_Triangle::num_post_elements(const ushort post_order) const
// {
//  const auto n = post_order;
//  return (n + 1) * (n + 1);
//}
//
// Quadrature_Rule Reference_Triangle::make_quadrature_rule(
//    const ushort integrand_order) const {
//  switch (integrand_order) {
//    case 0:
//      return {{{-0.5, 0}}, {2}};
//    case 1:
//    case 2:
//      return {{{-0.6666666666666667, -0.577350269189626},
//               {-0.9106836025229592, 0.577350269189626},
//               {0.2440169358562927, -0.577350269189626},
//               {-0.6666666666666667, 0.577350269189626}},
//              {0.788675134594813, 0.211324865405187, 0.788675134594813,
//               0.211324865405187}};
//    case 3:
//    case 4:
//      return {{{-0.7999999999999996, -0.774596669241483},
//               {-0.8872983346207415, 0},
//               {-0.9745966692414834, 0.774596669241483},
//               {-0.1127016653792585, -0.774596669241483},
//               {-0.5, 0},
//               {-0.8872983346207415, 0.774596669241483},
//               {0.5745966692414826, -0.774596669241483},
//               {-0.1127016653792585, 0},
//               {-0.7999999999999996, 0.774596669241483}},
//              {0.2738575106854125, 0.2469135802469129, 0.03478446462322775,
//               0.4381720170966613, 0.3950617283950618, 0.05565514339716456,
//               0.2738575106854125, 0.2469135802469129,
//               0.03478446462322775}};
//    case 5:
//    case 6:
//      return {{{-0.8707778735729041, -0.861136311594052},
//               {-0.9069626449468777, -0.339981043584856},
//               {-0.9541736666471743, 0.339981043584856},
//               {-0.9903584380211479, 0.861136311594052},
//               {-0.3858073769376817, -0.861136311594052},
//               {-0.5577935549985239, -0.339981043584856},
//               {-0.7821874885863321, 0.339981043584856},
//               {-0.9541736666471743, 0.861136311594052},
//               {0.2469436885317338, -0.861136311594052},
//               {-0.1022254014166201, -0.339981043584856},
//               {-0.5577935549985238, 0.339981043584856},
//               {-0.9069626449468777, 0.861136311594052},
//               {0.7319141851669562, -0.861136311594052},
//               {0.2469436885317338, -0.339981043584856},
//               {-0.3858073769376817, 0.339981043584856},
//               {-0.8707778735729041, 0.861136311594052}},
//              {0.1126015323077027, 0.1519885905918008, 0.07486326126005106,
//               0.008401460977899435, 0.211101109416918, 0.2849424819989601,
//               0.140350821011734, 0.01575074243493392, 0.211101109416918,
//               0.2849424819989601, 0.140350821011734, 0.01575074243493392,
//               0.1126015323077027, 0.1519885905918008, 0.07486326126005106,
//               0.008401460977899435}};
//    case 7:
//    case 8:
//      return {{{-0.9105809565927103, -0.906179845938664},
//               {-0.9278302861536236, -0.538469310105683},
//               {-0.9530899229693319, 0},
//               {-0.9783495597850402, 0.538469310105683},
//               {-0.9955988893459535, 0.906179845938664},
//               {-0.5601197503206428, -0.906179845938664},
//               {-0.644974598962845, -0.538469310105683},
//               {-0.7692346550528415, 0},
//               {-0.893494711142838, 0.538469310105683},
//               {-0.9783495597850402, 0.906179845938664},
//               {-0.04691007703066802, -0.906179845938664},
//               {-0.2307653449471585, -0.538469310105683},
//               {-0.5, 0},
//               {-0.7692346550528415, 0.538469310105683},
//               {-0.9530899229693319, 0.906179845938664},
//               {0.4662995962593067, -0.906179845938664},
//               {0.1834439090685281, -0.538469310105683},
//               {-0.2307653449471585, 0},
//               {-0.644974598962845, 0.538469310105683},
//               {-0.9278302861536237, 0.906179845938664},
//               {0.8167608025313744, -0.906179845938664},
//               {0.4662995962593067, -0.538469310105683},
//               {-0.04691007703066802, 0},
//               {-0.5601197503206428, 0.538469310105683},
//               {-0.9105809565927103, 0.906179845938664}},
//              {0.05350108223322567,  0.08723120988299211,
//              0.06739253619376044,
//               0.02616879011700774,  0.002633266629202917,
//               0.1080803972647221, 0.1762204318958822,   0.1361432662753753,
//               0.05286497232810846, 0.005319602735277746,
//               0.1284622942592381,   0.2094522369422109, 0.1618172839506173,
//               0.06283429560853965, 0.006322778128282769,
//               0.1080803972647221,   0.1762204318958822, 0.1361432662753753,
//               0.05286497232810846,  0.005319602735277746,
//               0.05350108223322567, 0.08723120988299211,
//               0.06739253619376044, 0.02616879011700774,
//               0.002633266629202917}};
//    case 9:
//    case 10:
//      return {
//          {{-0.9347496974591312, -0.9324695142031521},
//           {-0.9439088615608247, -0.661209386466264},
//           {-0.9581777223232526, -0.238619186083197},
//           {-0.9742917918798994, 0.238619186083197},
//           {-0.9885606526423274, 0.661209386466264},
//           {-0.9977198167440209, 0.9324695142031521},
//           {-0.6726487338239366, -0.9324695142031521},
//           {-0.7185989263755467, -0.661209386466264},
//           {-0.7901837230061085, -0.238619186083197},
//           {-0.8710256634601555, 0.238619186083197},
//           {-0.9426104600907174, 0.661209386466264},
//           {-0.9885606526423274, 0.9324695142031521},
//           {-0.2643273942032976, -0.9324695142031521},
//           {-0.3675935226230415, -0.661209386466264},
//           {-0.5284695579835037, -0.238619186083197},
//           {-0.7101496280996933, 0.238619186083197},
//           {-0.8710256634601555, 0.661209386466264},
//           {-0.9742917918798994, 0.9324695142031521},
//           {0.1967969084064496, -0.9324695142031521},
//           {0.0288029090893055, -0.661209386466264},
//           {-0.2329112559332993, -0.238619186083197},
//           {-0.5284695579835037, 0.238619186083197},
//           {-0.7901837230061085, 0.661209386466264},
//           {-0.9581777223232526, 0.9324695142031521},
//           {0.6051182480270888, -0.9324695142031521},
//           {0.3798083128418107, -0.661209386466264},
//           {0.0288029090893055, -0.238619186083197},
//           {-0.3675935226230415, 0.238619186083197},
//           {-0.7185989263755467, 0.661209386466264},
//           {-0.9439088615608248, 0.9324695142031521},
//           {0.8672192116622832, -0.9324695142031521},
//           {0.6051182480270888, -0.661209386466264},
//           {0.1967969084064496, -0.238619186083197},
//           {-0.2643273942032976, 0.238619186083197},
//           {-0.6726487338239368, 0.661209386466264},
//           {-0.9347496974591312, 0.9324695142031521}},
//          {0.02836100152117781, 0.0513374279511389,  0.049647026182223,
//           0.03051809113558391, 0.01046986542124473, 0.0009910801678028134,
//           0.05972035509877095, 0.1081022976149208,  0.1045427831942518,
//           0.06426258389333622, 0.02204661497324695, 0.002086938273612684,
//           0.0774583226595905,  0.1402105301458922,  0.1355937790222319,
//           0.08334967114506463, 0.02859483694169573, 0.002706794658216405,
//           0.0774583226595905,  0.1402105301458922,  0.1355937790222319,
//           0.08334967114506463, 0.02859483694169573, 0.002706794658216405,
//           0.05972035509877095, 0.1081022976149208,  0.1045427831942518,
//           0.06426258389333622, 0.02204661497324695, 0.002086938273612684,
//           0.02836100152117781, 0.0513374279511389,  0.049647026182223,
//           0.03051809113558391, 0.01046986542124473,
//           0.0009910801678028134}};
//    case 11:
//    case 12:
//      return {
//          {{-0.9504029146358142, -0.949107912342758},
//           {-0.9556849211223275, -0.741531185599394},
//           {-0.9642268026617964, -0.405845151377397},
//           {-0.974553956171379, 0},
//           {-0.9848811096809615, 0.405845151377397},
//           {-0.9934229912204304, 0.741531185599394},
//           {-0.9987049977069438, 0.949107912342758},
//           {-0.7481081943789636, -0.949107912342758},
//           {-0.7749342496082214, -0.741531185599394},
//           {-0.8183164352463219, -0.405845151377397},
//           {-0.870765592799697, 0},
//           {-0.9232147503530721, 0.405845151377397},
//           {-0.9665969359911726, 0.741531185599394},
//           {-0.9934229912204304, 0.949107912342758},
//           {-0.4209640416964355, -0.949107912342758},
//           {-0.4826304010243249, -0.741531185599394},
//           {-0.5823551434482712, -0.405845151377397},
//           {-0.7029225756886985, 0},
//           {-0.8234900079291259, 0.405845151377397},
//           {-0.9232147503530722, 0.741531185599394},
//           {-0.9848811096809615, 0.949107912342758},
//           {-0.02544604382862098, -0.949107912342758},
//           {-0.129234407200303, -0.741531185599394},
//           {-0.2970774243113015, -0.405845151377397},
//           {-0.5, 0},
//           {-0.7029225756886985, 0.405845151377397},
//           {-0.870765592799697, 0.741531185599394},
//           {-0.974553956171379, 0.949107912342758},
//           {0.3700719540391935, -0.949107912342758},
//           {0.2241615866237189, -0.741531185599394},
//           {-0.01179970517433182, -0.405845151377397},
//           {-0.2970774243113015, 0},
//           {-0.5823551434482711, 0.405845151377397},
//           {-0.8183164352463218, 0.741531185599394},
//           {-0.9642268026617964, 0.949107912342758},
//           {0.6972161067217216, -0.949107912342758},
//           {0.5164654352076156, -0.741531185599394},
//           {0.2241615866237189, -0.405845151377397},
//           {-0.129234407200303, 0},
//           {-0.4826304010243249, 0.405845151377397},
//           {-0.7749342496082215, 0.741531185599394},
//           {-0.9556849211223276, 0.949107912342758},
//           {0.8995108269785723, -0.949107912342758},
//           {0.6972161067217215, -0.741531185599394},
//           {0.3700719540391935, -0.405845151377397},
//           {-0.02544604382862098, 0},
//           {-0.4209640416964355, 0.405845151377397},
//           {-0.7481081943789636, 0.741531185599394},
//           {-0.9504029146358142, 0.949107912342758}},
//          {0.01633971902233046,   0.03153707751100931, 0.03475337161903316,
//           0.02705971537896383,   0.01468787955288011, 0.004680565643230263,
//           0.0004266374414229519, 0.03529604741916743, 0.06812443847836634,
//           0.07507207749196593,   0.05845271854796313, 0.03172784626693878,
//           0.01011066754980336,   0.0009215957350721348,
//           0.04818316692765089, 0.09299769892288511,   0.1024820257759689,
//           0.07979468810555948, 0.04331216169277281,   0.0138022248360196,
//           0.001258084244262363, 0.05274230535088141,   0.1017972322343419,
//           0.1121789753788724, 0.08734493960849631,   0.04741040083224651,
//           0.01510820486158434, 0.001377125407046246,  0.04818316692765089,
//           0.09299769892288511, 0.1024820257759689,    0.07979468810555948,
//           0.04331216169277281, 0.0138022248360196,    0.001258084244262363,
//           0.03529604741916743, 0.06812443847836634,   0.07507207749196593,
//           0.05845271854796313, 0.03172784626693878,   0.01011066754980336,
//           0.0009215957350721348, 0.01633971902233046, 0.03153707751100931,
//           0.03475337161903316, 0.02705971537896383,   0.01468787955288011,
//           0.004680565643230263, 0.0004266374414229519}};
//    case 13:
//    case 14:
//      return {
//          {{-0.9610783042460291, -0.960289856497536},
//           {-0.9643270581779191, -0.796666477413627},
//           {-0.9697104445422814, -0.525532409916329},
//           {-0.9765028202603553, -0.18343464249565},
//           {-0.9837870362371807, 0.18343464249565},
//           {-0.9905794119552546, 0.525532409916329},
//           {-0.9959627983196169, 0.796666477413627},
//           {-0.9992115522515068, 0.960289856497536},
//           {-0.8007036790940101, -0.960289856497536},
//           {-0.8173387381173185, -0.796666477413627},
//           {-0.844904060636017, -0.525532409916329},
//           {-0.8796840326953073, -0.18343464249565},
//           {-0.9169824447183197, 0.18343464249565},
//           {-0.95176241677761, 0.525532409916329},
//           {-0.9793277392963085, 0.796666477413627},
//           {-0.9959627983196169, 0.960289856497536},
//           {-0.5349529979610744, -0.960289856497536},
//           {-0.573769993138719, -0.796666477413627},
//           {-0.6380921569362322, -0.525532409916329},
//           {-0.719249308576779, -0.18343464249565},
//           {-0.8062831013395499, 0.18343464249565},
//           {-0.8874402529800968, 0.525532409916329},
//           {-0.95176241677761, 0.796666477413627},
//           {-0.9905794119552546, 0.960289856497536},
//           {-0.1996476062584693, -0.960289856497536},
//           {-0.2664521977773303, -0.796666477413627},
//           {-0.3771515411561002, -0.525532409916329},
//           {-0.5168241340337535, -0.18343464249565},
//           {-0.6666105084618966, 0.18343464249565},
//           {-0.8062831013395499, 0.525532409916329},
//           {-0.9169824447183198, 0.796666477413627},
//           {-0.9837870362371808, 0.960289856497536},
//           {0.1599374627560052, -0.960289856497536},
//           {0.06311867519095721, -0.796666477413627},
//           {-0.0973160489275709, -0.525532409916329},
//           {-0.2997412234705965, -0.18343464249565},
//           {-0.5168241340337535, 0.18343464249565},
//           {-0.719249308576779, 0.525532409916329},
//           {-0.8796840326953073, 0.796666477413627},
//           {-0.9765028202603552, 0.960289856497536},
//           {0.4952428544586104, -0.960289856497536},
//           {0.370436470552346, -0.796666477413627},
//           {0.1636245668525612, -0.525532409916329},
//           {-0.0973160489275709, -0.18343464249565},
//           {-0.3771515411561001, 0.18343464249565},
//           {-0.6380921569362322, 0.525532409916329},
//           {-0.844904060636017, 0.796666477413627},
//           {-0.9697104445422814, 0.960289856497536},
//           {0.760993535591546, -0.960289856497536},
//           {0.6140052155309454, -0.796666477413627},
//           {0.370436470552346, -0.525532409916329},
//           {0.06311867519095721, -0.18343464249565},
//           {-0.2664521977773303, 0.18343464249565},
//           {-0.573769993138719, 0.525532409916329},
//           {-0.8173387381173185, 0.796666477413627},
//           {-0.9643270581779191, 0.960289856497536},
//           {0.9213681607435651, -0.960289856497536},
//           {0.760993535591546, -0.796666477413627},
//           {0.4952428544586104, -0.525532409916329},
//           {0.1599374627560052, -0.18343464249565},
//           {-0.1996476062584693, 0.18343464249565},
//           {-0.5349529979610744, 0.525532409916329},
//           {-0.8007036790940101, 0.796666477413627},
//           {-0.9610783042460291, 0.960289856497536}},
//          {0.01004375733945304,   0.02022265498028209, 0.02422245286926617,
//           0.02172427927521026,   0.01498966925243749, 0.007533611717515962,
//           0.002288651636172856,  0.00020345922003913, 0.02206434300837125,
//           0.0444255651490272,    0.05321240752324866, 0.04772436582622505,
//           0.0329296291009185,    0.01655000090197412, 0.005027759335525518,
//           0.0004469636080836971, 0.03112554564587481, 0.06266989030061787,
//           0.07506524072180071,   0.06732341526693157, 0.04645289793099652,
//           0.02334661894615327,   0.00709251812460491,
//           0.0006305189409073175, 0.03598499044536026, 0.07245416447754377,
//           0.08678472663211513, 0.07783421639230392,   0.05370531033333868,
//           0.02699158656581295, 0.008199830449599781, 0.0007289580822874853,
//           0.03598499044536026, 0.07245416447754377,   0.08678472663211513,
//           0.07783421639230392, 0.05370531033333868,   0.02699158656581295,
//           0.008199830449599781, 0.0007289580822874853, 0.03112554564587481,
//           0.06266989030061787, 0.07506524072180071,   0.06732341526693157,
//           0.04645289793099652, 0.02334661894615327,   0.00709251812460491,
//           0.0006305189409073175, 0.02206434300837125,   0.0444255651490272,
//           0.05321240752324866, 0.04772436582622505,   0.0329296291009185,
//           0.01655000090197412, 0.005027759335525518, 0.0004469636080836971,
//           0.01004375733945304, 0.02022265498028209,   0.02422245286926617,
//           0.02172427927521026, 0.01498966925243749,   0.007533611717515962,
//           0.002288651636172856, 0.00020345922003913}};
//    case 15:
//    case 16:
//      return {
//          {{-0.9686671246817318, -0.968160239507626},
//           {-0.9707706046430857, -0.836031107326636},
//           {-0.9743153199987874, -0.61337143270059},
//           {-0.9789180440838081, -0.324253423403809},
//           {-0.9840801197538129, 0},
//           {-0.9892421954238177, 0.324253423403809},
//           {-0.9938449195088385, 0.61337143270059},
//           {-0.9973896348645401, 0.836031107326636},
//           {-0.9994931148258941, 0.968160239507626},
//           {-0.8386414724620959, -0.968160239507626},
//           {-0.8494740062089006, -0.836031107326636},
//           {-0.8677286363546227, -0.61337143270059},
//           {-0.891431816272783, -0.324253423403809},
//           {-0.918015553663318, 0},
//           {-0.944599291053853, 0.324253423403809},
//           {-0.9683024709720133, 0.61337143270059},
//           {-0.9865571011177354, 0.836031107326636},
//           {-0.9973896348645401, 0.968160239507626},
//           {-0.6195265131917514, -0.968160239507626},
//           {-0.6450689617285768, -0.836031107326636},
//           {-0.6881122572265872, -0.61337143270059},
//           {-0.7440028980840231, -0.324253423403809},
//           {-0.806685716350295, 0},
//           {-0.8693685346165668, 0.324253423403809},
//           {-0.9252591754740027, 0.61337143270059},
//           {-0.9683024709720132, 0.836031107326636},
//           {-0.9938449195088385, 0.968160239507626},
//           {-0.3350112279799912, -0.968160239507626},
//           {-0.379654132349956, -0.836031107326636},
//           {-0.4548848887872422, -0.61337143270059},
//           {-0.552570141294545, -0.324253423403809},
//           {-0.6621267117019045, 0},
//           {-0.7716832821092641, 0.324253423403809},
//           {-0.8693685346165668, 0.61337143270059},
//           {-0.944599291053853, 0.836031107326636},
//           {-0.9892421954238177, 0.968160239507626},
//           {-0.01591988024618701, -0.968160239507626},
//           {-0.081984446336682, -0.836031107326636},
//           {-0.193314283649705, -0.61337143270059},
//           {-0.3378732882980955, -0.324253423403809},
//           {-0.5, 0},
//           {-0.6621267117019045, 0.324253423403809},
//           {-0.806685716350295, 0.61337143270059},
//           {-0.918015553663318, 0.836031107326636},
//           {-0.9840801197538129, 0.968160239507626},
//           {0.3031714674876171, -0.968160239507626},
//           {0.215685239676592, -0.836031107326636},
//           {0.06825632148783223, -0.61337143270059},
//           {-0.123176435301646, -0.324253423403809},
//           {-0.3378732882980955, 0},
//           {-0.552570141294545, 0.324253423403809},
//           {-0.7440028980840232, 0.61337143270059},
//           {-0.891431816272783, 0.836031107326636},
//           {-0.9789180440838081, 0.968160239507626},
//           {0.5876867526993774, -0.968160239507626},
//           {0.4811000690552128, -0.836031107326636},
//           {0.3014836899271773, -0.61337143270059},
//           {0.06825632148783223, -0.324253423403809},
//           {-0.193314283649705, 0},
//           {-0.4548848887872422, 0.324253423403809},
//           {-0.6881122572265872, 0.61337143270059},
//           {-0.8677286363546228, 0.836031107326636},
//           {-0.9743153199987875, 0.968160239507626},
//           {0.8068017119697218, -0.968160239507626},
//           {0.6855051135355366, -0.836031107326636},
//           {0.4811000690552127, -0.61337143270059},
//           {0.215685239676592, -0.324253423403809},
//           {-0.081984446336682, 0},
//           {-0.379654132349956, 0.324253423403809},
//           {-0.6450689617285768, 0.61337143270059},
//           {-0.8494740062089006, 0.836031107326636},
//           {-0.9707706046430858, 0.968160239507626},
//           {0.9368273641893579, -0.968160239507626},
//           {0.8068017119697217, -0.836031107326636},
//           {0.5876867526993774, -0.61337143270059},
//           {0.3031714674876172, -0.324253423403809},
//           {-0.01591988024618701, 0},
//           {-0.3350112279799912, 0.324253423403809},
//           {-0.6195265131917516, 0.61337143270059},
//           {-0.8386414724620959, 0.836031107326636},
//           {-0.9686671246817318, 0.968160239507626}},
//          {0.006500367017424581, 0.01347836749000478,  0.01708638995104882,
//           0.01680862795979199,  0.01342000079532422,  0.008577189683159995,
//           0.004094584999583912, 0.001203701279113231,
//           0.0001051591861235364, 0.01444833199254737,  0.02995829738399937,
//           0.03797783016022493, 0.0373604500255599,   0.02982856603501678,
//           0.01906447494013144, 0.009101012802371234, 0.002675460578435511,
//           0.0002337367765706412, 0.02084377636592117,  0.04321911008813613,
//           0.05478842811273874, 0.05389776935251935,  0.04303195414326739,
//           0.02750321991429733, 0.01312950696688454,  0.003859732874460045,
//           0.00033719858471156, 0.02498167846612465,  0.05179895873279031,
//           0.06566501533832468, 0.06459754318835401,  0.05157464862910972,
//           0.03296315334707955, 0.01573597392849199,  0.004625966232901031,
//           0.000404139176827337, 0.02641271197951784,  0.05476617512723753,
//           0.0694265255080294, 0.06829790500794712,  0.05452901579582412,
//           0.03485139225027233, 0.01663738277850538,  0.004890956942796017,
//           0.0004272896111305921, 0.02498167846612465,  0.05179895873279031,
//           0.06566501533832468, 0.06459754318835401,  0.05157464862910972,
//           0.03296315334707955, 0.01573597392849199,  0.004625966232901031,
//           0.000404139176827337, 0.02084377636592117,  0.04321911008813613,
//           0.05478842811273874, 0.05389776935251935,  0.04303195414326739,
//           0.02750321991429733, 0.01312950696688454,  0.003859732874460045,
//           0.00033719858471156, 0.01444833199254737,  0.02995829738399937,
//           0.03797783016022493, 0.0373604500255599,   0.02982856603501678,
//           0.01906447494013144, 0.009101012802371234, 0.002675460578435511,
//           0.0002337367765706412, 0.006500367017424581, 0.01347836749000478,
//           0.01708638995104882, 0.01680862795979199,  0.01342000079532422,
//           0.008577189683159995, 0.004094584999583912, 0.001203701279113231,
//           0.0001051591861235364}};
//    case 17:
//    case 18:
//      return {
//          {{-0.9742469631441846, -0.973906528517172},
//           {-0.9756670111138168, -0.865063366688985},
//           {-0.9780891871608004, -0.679409568299024},
//           {-0.9812988690798357, -0.433395394129247},
//           {-0.985010940099215, -0.148874338981631},
//           {-0.988895588417957, 0.148874338981631},
//           {-0.9926076594373363, 0.433395394129247},
//           {-0.9958173413563716, 0.679409568299024},
//           {-0.9982395174033551, 0.865063366688985},
//           {-0.9996595653729874, 0.973906528517172},
//           {-0.8668238492856299, -0.973906528517172},
//           {-0.8741673141936407, -0.865063366688985},
//           {-0.8866930634517123, -0.679409568299024},
//           {-0.903291225656342, -0.433395394129247},
//           {-0.9224873823002004, -0.148874338981631},
//           {-0.9425759843887846, 0.148874338981631},
//           {-0.9617721410326431, 0.433395394129247},
//           {-0.9783703032372728, 0.679409568299024},
//           {-0.9908960524953444, 0.865063366688985},
//           {-0.9982395174033551, 0.973906528517172},
//           {-0.6835922269426524, -0.973906528517172},
//           {-0.7010392650617513, -0.865063366688985},
//           {-0.730798680748133, -0.679409568299024},
//           {-0.770233575898957, -0.433395394129247},
//           {-0.8158409398478528, -0.148874338981631},
//           {-0.8635686284511712, 0.148874338981631},
//           {-0.909175992400067, 0.433395394129247},
//           {-0.948610887550891, 0.679409568299024},
//           {-0.9783703032372727, 0.865063366688985},
//           {-0.9958173413563716, 0.973906528517172},
//           {-0.4407877346919107, -0.973906528517172},
//           {-0.471623253096604, -0.865063366688985},
//           {-0.52421940172918, -0.679409568299024},
//           {-0.5939157838262227, -0.433395394129247},
//           {-0.6745212539831456, -0.148874338981631},
//           {-0.7588741401461014, 0.148874338981631},
//           {-0.8394796103030243, 0.433395394129247},
//           {-0.909175992400067, 0.679409568299024},
//           {-0.961772141032643, 0.865063366688985},
//           {-0.9926076594373363, 0.973906528517172},
//           {-0.1599787505636739, -0.973906528517172},
//           {-0.2062983545928464, -0.865063366688985},
//           {-0.2853057105304597, -0.679409568299024},
//           {-0.3900001988355295, -0.433395394129247},
//           {-0.5110817844036087, -0.148874338981631},
//           {-0.6377925545780222, 0.148874338981631},
//           {-0.7588741401461014, 0.433395394129247},
//           {-0.8635686284511712, 0.679409568299024},
//           {-0.9425759843887845, 0.865063366688985},
//           {-0.988895588417957, 0.973906528517172},
//           {0.1338852790808459, -0.973906528517172},
//           {0.07136172128183138, -0.865063366688985},
//           {-0.03528472117051629, -0.679409568299024},
//           {-0.1766044070352235, -0.433395394129247},
//           {-0.3400438766147603, -0.148874338981631},
//           {-0.5110817844036089, 0.148874338981631},
//           {-0.6745212539831456, 0.433395394129247},
//           {-0.8158409398478528, 0.679409568299024},
//           {-0.9224873823002004, 0.865063366688985},
//           {-0.985010940099215, 0.973906528517172},
//           {0.4146942632090826, -0.973906528517172},
//           {0.336686619785589, -0.865063366688985},
//           {0.203628970028204, -0.679409568299024},
//           {0.02731117795546967, -0.433395394129247},
//           {-0.1766044070352235, -0.148874338981631},
//           {-0.3900001988355296, 0.148874338981631},
//           {-0.5939157838262227, 0.433395394129247},
//           {-0.770233575898957, 0.679409568299024},
//           {-0.903291225656342, 0.865063366688985},
//           {-0.9812988690798357, 0.973906528517172},
//           {0.6574987554598244, -0.973906528517172},
//           {0.5661026317507363, -0.865063366688985},
//           {0.410208249047157, -0.679409568299024},
//           {0.203628970028204, -0.433395394129247},
//           {-0.03528472117051629, -0.148874338981631},
//           {-0.2853057105304597, 0.148874338981631},
//           {-0.52421940172918, 0.433395394129247},
//           {-0.730798680748133, 0.679409568299024},
//           {-0.8866930634517123, 0.865063366688985},
//           {-0.9780891871608004, 0.973906528517172},
//           {0.8407303778028019, -0.973906528517172},
//           {0.7392306808826257, -0.865063366688985},
//           {0.5661026317507363, -0.679409568299024},
//           {0.3366866197855889, -0.433395394129247},
//           {0.07136172128183144, -0.148874338981631},
//           {-0.2062983545928465, 0.148874338981631},
//           {-0.471623253096604, 0.433395394129247},
//           {-0.7010392650617512, 0.679409568299024},
//           {-0.8741673141936406, 0.865063366688985},
//           {-0.9756670111138168, 0.973906528517172},
//           {0.9481534916613564, -0.973906528517172},
//           {0.8407303778028018, -0.865063366688985},
//           {0.6574987554598245, -0.679409568299024},
//           {0.4146942632090827, -0.433395394129247},
//           {0.133885279080846, -0.148874338981631},
//           {-0.159978750563674, 0.148874338981631},
//           {-0.4407877346919107, 0.433395394129247},
//           {-0.6835922269426525, 0.679409568299024},
//           {-0.8668238492856298, 0.865063366688985},
//           {-0.9742469631441845, 0.973906528517172}},
//          {0.004387074522396848,  0.00929185979426592, 0.01226538498559636,
//           0.01286642521300538,   0.01131813402104741, 0.008384863316467971,
//           0.005085948940982216,  0.00234139732304471,
//           0.0006722625623504124, 5.799362953077529e-05,
//           0.009834123085334445,  0.02082875329379134, 0.02749424588563135,
//           0.0288415454460565,    0.02537087585156437, 0.01879561823873495,
//           0.01140072903617321,   0.005248506572875442,
//           0.001506952469137529,  0.0001299992712818888,
//           0.01441621147982787, 0.03053365406746336,   0.04030485074533407,
//           0.04227990792340959, 0.03719212262556008,   0.02755320480255082,
//           0.01671275815682937, 0.007693983495150223,  0.002209098391043433,
//           0.0001905708288132014, 0.01771815427246953, 0.03752719596452479,
//           0.04953642393731129, 0.05196385557058451,   0.0457107449708518,
//           0.03386409349471978, 0.02054071055738368,   0.009456242142927666,
//           0.002715078517704924, 0.0002342198815180672, 0.01944593753793903,
//           0.0411866550814514, 0.05436696119271134,   0.05703110347256456,
//           0.05016822169208674, 0.03716634570116908,   0.02254373499300701,
//           0.01037836623539956, 0.002979839008847915, 0.0002570597995763472,
//           0.01944593753793903, 0.0411866550814514,    0.05436696119271134,
//           0.05703110347256456, 0.05016822169208674,   0.03716634570116908,
//           0.02254373499300701, 0.01037836623539956,   0.002979839008847915,
//           0.0002570597995763472, 0.01771815427246953, 0.03752719596452479,
//           0.04953642393731129, 0.05196385557058451,   0.0457107449708518,
//           0.03386409349471978, 0.02054071055738368,   0.009456242142927666,
//           0.002715078517704924, 0.0002342198815180672, 0.01441621147982787,
//           0.03053365406746336, 0.04030485074533407,   0.04227990792340959,
//           0.03719212262556008, 0.02755320480255082,   0.01671275815682937,
//           0.007693983495150223, 0.002209098391043433,
//           0.0001905708288132014, 0.009834123085334445, 0.02082875329379134,
//           0.02749424588563135,   0.0288415454460565, 0.02537087585156437,
//           0.01879561823873495,   0.01140072903617321, 0.005248506572875442,
//           0.001506952469137529, 0.0001299992712818888,
//           0.004387074522396848,  0.00929185979426592, 0.01226538498559636,
//           0.01286642521300538,   0.01131813402104741, 0.008384863316467971,
//           0.005085948940982216,  0.00234139732304471,
//           0.0006722625623504124, 5.799362953077529e-05}};
//    case 19:
//    case 20:
//      return {
//          {{-0.9784656538091175, -0.978228658146057},
//           {-0.9794580575203291, -0.887062599768095},
//           {-0.9811661346136811, -0.730152005574049},
//           {-0.9834636194310185, -0.519096129206812},
//           {-0.9861801709767138, -0.269543155952345},
//           {-0.9891143290730284, 0},
//           {-0.992048487169343, 0.269543155952345},
//           {-0.9947650387150384, 0.519096129206812},
//           {-0.9970625235323758, 0.730152005574049},
//           {-0.9987706006257278, 0.887062599768095},
//           {-0.9997630043369393, 0.978228658146057},
//           {-0.8882919991423672, -0.978228658146057},
//           {-0.8934400279536657, -0.887062599768095},
//           {-0.9023005652422252, -0.730152005574049},
//           {-0.9142186162325163, -0.519096129206812},
//           {-0.9283105482422671, -0.269543155952345},
//           {-0.9435312998840475, 0},
//           {-0.9587520515258279, 0.269543155952345},
//           {-0.9728439835355787, 0.519096129206812},
//           {-0.9847620345258697, 0.730152005574049},
//           {-0.9936225718144293, 0.887062599768095},
//           {-0.9987706006257278, 0.978228658146057},
//           {-0.7330894820416732, -0.978228658146057},
//           {-0.7453899710481793, -0.887062599768095},
//           {-0.7665609756219032, -0.730152005574049},
//           {-0.7950374780966583, -0.519096129206812},
//           {-0.8287081627645337, -0.269543155952345},
//           {-0.8650760027870246, 0},
//           {-0.9014438428095154, 0.269543155952345},
//           {-0.9351145274773909, 0.519096129206812},
//           {-0.963591029952146, 0.730152005574049},
//           {-0.9847620345258699, 0.887062599768095},
//           {-0.9970625235323759, 0.978228658146057},
//           {-0.5243310904917735, -0.978228658146057},
//           {-0.5462521456712334, -0.887062599768095},
//           {-0.5839816017294213, -0.730152005574049},
//           {-0.6347303956787477, -0.519096129206812},
//           {-0.6947358910817587, -0.269543155952345},
//           {-0.7595480646034061, 0},
//           {-0.8243602381250534, 0.269543155952345},
//           {-0.8843657335280645, 0.519096129206812},
//           {-0.9351145274773909, 0.730152005574049},
//           {-0.9728439835355788, 0.887062599768095},
//           {-0.9947650387150386, 0.978228658146057},
//           {-0.2774946687830019, -0.978228658146057},
//           {-0.3107911044265171, -0.887062599768095},
//           {-0.3680993131428296, -0.730152005574049},
//           {-0.4451829178272916, -0.519096129206812},
//           {-0.5363267564603751, -0.269543155952345},
//           {-0.6347715779761725, 0},
//           {-0.7332163994919698, 0.269543155952345},
//           {-0.8243602381250532, 0.519096129206812},
//           {-0.9014438428095153, 0.730152005574049},
//           {-0.9587520515258279, 0.887062599768095},
//           {-0.992048487169343, 0.978228658146057},
//           {-0.01088567092697151, -0.978228658146057},
//           {-0.05646870011595251, -0.887062599768095},
//           {-0.1349239972129755, -0.730152005574049},
//           {-0.240451935396594, -0.519096129206812},
//           {-0.3652284220238275, -0.269543155952345},
//           {-0.5, 0},
//           {-0.6347715779761725, 0.269543155952345},
//           {-0.7595480646034061, 0.519096129206812},
//           {-0.8650760027870246, 0.730152005574049},
//           {-0.9435312998840475, 0.887062599768095},
//           {-0.9891143290730284, 0.978228658146057},
//           {0.2557233269290589, -0.978228658146057},
//           {0.1978537041946121, -0.887062599768095},
//           {0.0982513187168787, -0.730152005574049},
//           {-0.03572095296589628, -0.519096129206812},
//           {-0.1941300875872799, -0.269543155952345},
//           {-0.3652284220238275, 0},
//           {-0.5363267564603751, 0.269543155952345},
//           {-0.6947358910817587, 0.519096129206812},
//           {-0.8287081627645336, 0.730152005574049},
//           {-0.9283105482422671, 0.887062599768095},
//           {-0.986180170976714, 0.978228658146057},
//           {0.5025597486378306, -0.978228658146057},
//           {0.4333147454393283, -0.887062599768095},
//           {0.3141336073034702, -0.730152005574049},
//           {0.1538265248855596, -0.519096129206812},
//           {-0.03572095296589628, -0.269543155952345},
//           {-0.240451935396594, 0},
//           {-0.4451829178272917, 0.269543155952345},
//           {-0.6347303956787476, 0.519096129206812},
//           {-0.7950374780966583, 0.730152005574049},
//           {-0.9142186162325163, 0.887062599768095},
//           {-0.9834636194310185, 0.978228658146057},
//           {0.7113181401877302, -0.978228658146057},
//           {0.6324525708162743, -0.887062599768095},
//           {0.496712981195952, -0.730152005574049},
//           {0.3141336073034702, -0.519096129206812},
//           {0.0982513187168787, -0.269543155952345},
//           {-0.1349239972129755, 0},
//           {-0.3680993131428297, 0.269543155952345},
//           {-0.5839816017294213, 0.519096129206812},
//           {-0.766560975621903, 0.730152005574049},
//           {-0.9023005652422251, 0.887062599768095},
//           {-0.9811661346136811, 0.978228658146057},
//           {0.8665206572884241, -0.978228658146057},
//           {0.7805026277217607, -0.887062599768095},
//           {0.6324525708162743, -0.730152005574049},
//           {0.4333147454393284, -0.519096129206812},
//           {0.1978537041946121, -0.269543155952345},
//           {-0.05646870011595251, 0},
//           {-0.3107911044265171, 0.269543155952345},
//           {-0.5462521456712334, 0.519096129206812},
//           {-0.7453899710481793, 0.730152005574049},
//           {-0.8934400279536657, 0.887062599768095},
//           {-0.9794580575203291, 0.978228658146057},
//           {0.9566943119551745, -0.978228658146057},
//           {0.8665206572884241, -0.887062599768095},
//           {0.7113181401877302, -0.730152005574049},
//           {0.5025597486378304, -0.519096129206812},
//           {0.2557233269290589, -0.269543155952345},
//           {-0.01088567092697151, 0},
//           {-0.2774946687830019, 0.269543155952345},
//           {-0.5243310904917735, 0.519096129206812},
//           {-0.7330894820416731, 0.730152005574049},
//           {-0.8882919991423672, 0.887062599768095},
//           {-0.9784656538091177, 0.978228658146057}},
//          {0.003065254786336921,  0.006596113363469354,
//          0.008971278567846241,
//           0.009860120851096311,  0.009286677986218874,
//           0.007596674255491598, 0.005343274438285346, 0.003121441884166165,
//           0.001399230542270533,
//           0.0003947658625615832, 3.37345784310486e-05,
//           0.006914778815286163, 0.01487989355803276,   0.02023792843044774,
//           0.02224303019808677, 0.02094942465785366,   0.0171370166169049,
//           0.01205366713879896, 0.007041528916287183,  0.003156465085552,
//           0.0008905356369090305, 7.610041074477408e-05,
//           0.01025761916059888,   0.02207334252415016, 0.03002163452865245,
//           0.03299607100161919,   0.03107709234297891, 0.02542163599166264,
//           0.01788082168660207,   0.01044564459125499, 0.004682408158847183,
//           0.001321050991849573, 0.0001128899495178914, 0.01284024971520857,
//           0.02763089812771649,   0.03758038567929434, 0.04130371625698049,
//           0.03890158328739785,   0.03182221421866715, 0.02238279779882984,
//           0.01307561558760397,   0.005861329913579828,
//           0.001653660986657467,  0.0001413130200539035,
//           0.01447069557673382, 0.03113945010308819,   0.0423523165735007,
//           0.04654843304446183, 0.04384127892295795,   0.03586297655804295,
//           0.02522494969228044, 0.01473594804176587,   0.006605597456080277,
//           0.001863641693564429, 0.000159256847770402,  0.01502795871881384,
//           0.0323386231293656, 0.04398329449594855,   0.04834100244236725,
//           0.04552959644134281, 0.0372440514963624,    0.02619635667474308,
//           0.01530342599496706, 0.006859977487376735,  0.001935410104444195,
//           0.0001653897921693557, 0.01447069557673382, 0.03113945010308819,
//           0.0423523165735007, 0.04654843304446183,   0.04384127892295795,
//           0.03586297655804295, 0.02522494969228044,   0.01473594804176587,
//           0.006605597456080277, 0.001863641693564429, 0.000159256847770402,
//           0.01284024971520857, 0.02763089812771649,   0.03758038567929434,
//           0.04130371625698049, 0.03890158328739785,   0.03182221421866715,
//           0.02238279779882984, 0.01307561558760397,   0.005861329913579828,
//           0.001653660986657467, 0.0001413130200539035, 0.01025761916059888,
//           0.02207334252415016, 0.03002163452865245,   0.03299607100161919,
//           0.03107709234297891, 0.02542163599166264,   0.01788082168660207,
//           0.01044564459125499, 0.004682408158847183,  0.001321050991849573,
//           0.0001128899495178914, 0.006914778815286163, 0.01487989355803276,
//           0.02023792843044774, 0.02224303019808677,   0.02094942465785366,
//           0.0171370166169049, 0.01205366713879896,   0.007041528916287183,
//           0.003156465085552, 0.0008905356369090305, 7.610041074477408e-05,
//           0.003065254786336921, 0.006596113363469354, 0.008971278567846241,
//           0.009860120851096311, 0.009286677986218874, 0.007596674255491598,
//           0.005343274438285346, 0.003121441884166165, 0.001399230542270533,
//           0.0003947658625615832, 3.37345784310486e-05}};
//    default:
//      EXCEPTION("not supported integrand order");
//      return {};
//  }
//};
//
// bool Reference_Triangle::is_simplex(void) const { return true; };
//
// Vector_Function<Polynomial>
// Reference_Triangle::make_normal_vector_function(
//    const Vector_Function<Polynomial>& mapping_function) const {
//  constexpr ushort r = 0;
//  constexpr ushort s = 1;
//  const auto T_r = mapping_function.get_differentiate(r);
//  const auto T_s = mapping_function.get_differentiate(s);
//
//  return T_r.cross_product(T_s);
//}
//
// Euclidean_Vector Reference_Triangle::random_point(void) const {
//  std::random_device rd;
//  std::mt19937 gen(rd());
//
//  std::uniform_real_distribution<double> x_distribution(-1.0, 1.0);
//  const auto x_coordinate = x_distribution(gen);
//
//  std::uniform_real_distribution<double> y_distribution(-1.0,
//  -x_coordinate); const auto y_coordinate = y_distribution(gen);
//
//  return {x_coordinate, y_coordinate};
//}
//
// bool Reference_Triangle::is_line(void) const { return false; };
//
// std::vector<std::vector<ushort>>
// Reference_Triangle::set_of_face_vertex_index_sequences(void) const {
//  //      2
//  //  2  / \  1
//  //	  /   \
//	//   0����������1
//  //      0
//  const std::vector<ushort> face0_node_index = {0, 1};
//  const std::vector<ushort> face1_node_index = {1, 2};
//  const std::vector<ushort> face2_node_index = {2, 0};
//  return {face0_node_index, face1_node_index, face2_node_index};
//};
//
// std::vector<std::vector<ushort>>
// Reference_Triangle::set_of_face_node_index_sequences(void) const {
//  //      2
//  //  2  / \  1
//  //	  /   \
//	//   0����������1
//  //      0
//  constexpr ushort num_face = 3;
//  std::vector<std::vector<ushort>> set_of_face_node_index_orders(num_face);
//  set_of_face_node_index_orders[0] = {0, 1};
//  set_of_face_node_index_orders[1] = {1, 2};
//  set_of_face_node_index_orders[2] = {2, 0};
//
//  if (this->order_ > 1) {
//    const auto num_additional_point = this->order_ - 1;
//
//    ushort index = num_face;
//    for (ushort iface = 0; iface < num_face; ++iface)
//      for (ushort ipoint = 0; ipoint < num_additional_point; ++ipoint)
//        set_of_face_node_index_orders[iface].push_back(index++);
//  }
//
//  return set_of_face_node_index_orders;
//};
//
// std::vector<std::shared_ptr<const Reference_Geometry>>
// Reference_Triangle::sub_simplex_reference_geometries(void) const {
//  EXCEPTION("This routine is not supported for Triangle element");
//  return {};
//}
//
// std::vector<std::vector<ushort>>
// Reference_Triangle::set_of_sub_simplex_vertex_index_sequences(void) const {
//  EXCEPTION("This routine is not supported for Triangle element");
//  return {};
//}
//
// Irrational_Function Reference_Triangle::scale_function(
//    const Vector_Function<Polynomial>& mapping_function) const {
//  constexpr ushort r = 0;
//  constexpr ushort s = 1;
//  const auto mf_r = mapping_function.get_differentiate(r);
//  const auto mf_s = mapping_function.get_differentiate(s);
//  const auto cross_product = mf_r.cross_product(mf_s);
//  return cross_product.L2_norm();
//};
//
// ushort Reference_Triangle::scale_function_order(void) const {
//  REQUIRE(this->order_ == 1, "high order mesh is not supported yet");
//  return 0;
//}
//
// std::vector<Euclidean_Vector> Reference_Triangle::make_mapping_points(
//    void) const {
//  //	  2
//  //    �� \ 
//	//	  ��   \
//	//    0����������1
//
//  switch (this->order_) {
//    case 1:
//      return {{-1, -1}, {1, -1}, {-1, 1}};
//    default:
//      EXCEPTION("unsuported figure order");
//      return {};
//  }
//}
//
// std::vector<Euclidean_Vector> Reference_Triangle::make_post_points(
//    const ushort post_order) const {
//  const auto n = post_order;
//  const auto num_post_point = static_cast<ushort>((n + 2) * (n + 3) * 0.5);
//
//  std::vector<Euclidean_Vector> post_points;
//  post_points.reserve(num_post_point);
//
//  const auto num_division = n + 1;
//  const auto delta = 2.0 / num_division;
//
//  const auto X0_start_coord = -1.0;
//  const auto X1_start_coord = -1.0;
//
//  for (ushort j = 0; j <= n + 1; ++j) {
//    for (ushort i = 0; i <= n + 1 - j; ++i) {
//      const auto X0_coord = X0_start_coord + delta * i;
//      const auto X1_coord = X1_start_coord + delta * j;
//      post_points.push_back({X0_coord, X1_coord});
//    }
//  }
//
//  return post_points;
//}
//
// std::vector<std::vector<uint>> Reference_Triangle::make_connectivities(
//    const ushort post_order) const {
//  const auto n = post_order;
//  const auto num_connecitivity = (n + 1) * (n + 1);
//
//  std::vector<std::vector<uint>> connectivities;
//  connectivities.reserve(num_connecitivity);
//
//  for (ushort j = 0; j < n + 1; j++) {
//    for (ushort i = 0; i < n + 1 - j; i++) {
//      //   I2��������I2+1
//      //   ��     ��
//      //   I1��������I1+1
//
//      const uint I1 = (n + 3) * j - j * (j + 1) / 2 + i;
//      const uint I2 = (n + 3) * (j + 1) - (j + 1) * (j + 2) / 2 + i;
//
//      if (i == n - j)
//        connectivities.push_back({I1, I1 + 1, I2});
//      else {
//        auto quadrilateral_connectivities =
//            this->quadrilateral_connectivities({I1, I1 + 1, I2, I2 + 1});
//
//        for (auto& connectivity : quadrilateral_connectivities)
//          connectivities.push_back(std::move(connectivity));
//      }
//    }
//  }
//
//  return connectivities;
//}
//
// Vector_Function<Polynomial>
// Reference_Triangle::make_mapping_monomial_vector_function(void) const {
//  const auto num_monomial =
//      static_cast<size_t>((this->order_ + 2) * (this->order_ + 1) * 0.5);
//  ms::sym::Polynomials mapping_monomial_vector(num_monomial);
//
//  Polynomial r("x0");
//  Polynomial s("x1");
//  for (ushort a = 0, index = 0; a <= this->order_; ++a)
//    for (ushort b = 0; b <= a; ++b)
//      mapping_monomial_vector[index++] = (r ^ (a - b)) * (s ^ b);
//
//  return mapping_monomial_vector;  // 1 r s r^2 rs s^2 ...
//}
//
// Reference_Quadrilateral::Reference_Quadrilateral(const ushort order) {
//  this->order_ = order;
//  this->initialize();
//}
//
// Euclidean_Vector Reference_Quadrilateral::center_point(void) const {
//  return {0.0, 0.0};
//};
//
// std::vector<std::shared_ptr<const Reference_Geometry>>
// Reference_Quadrilateral::face_reference_geometries(void) const {
//  constexpr auto num_face = 4;
//  std::vector<std::shared_ptr<const Reference_Geometry>>
//      face_reference_geometries(num_face);
//
//  for (ushort i = 0; i < num_face; ++i) {
//    face_reference_geometries[i] =
//    Reference_Geometry_Container::get_shared_ptr(
//        Figure::line, this->order_);
//  }
//
//  return face_reference_geometries;
//}
//
// Figure Reference_Quadrilateral::figure(void) const {
//  return Figure::quadrilateral;
//}
//
// ushort Reference_Quadrilateral::num_vertices(void) const { return 4; }
//
// ushort Reference_Quadrilateral::num_post_nodes(const ushort post_order)
// const
// {
//  const auto n = post_order;
//  return (n + 2) * (n + 2);
//}
//
// ushort Reference_Quadrilateral::num_post_elements(
//    const ushort post_order) const {
//  const auto n = post_order;
//  return 2 * (n + 1) * (n + 1);
//}
//
// Vector_Function<Polynomial>
// Reference_Quadrilateral::make_normal_vector_function(
//    const Vector_Function<Polynomial>& mapping_function) const {
//  constexpr ushort r = 0;
//  constexpr ushort s = 1;
//  const auto T_r = mapping_function.get_differentiate(r);
//  const auto T_s = mapping_function.get_differentiate(s);
//
//  return T_r.cross_product(T_s);
//}
//
// Euclidean_Vector Reference_Quadrilateral::random_point(void) const {
//  std::random_device rd;
//  std::mt19937 gen(rd());
//
//  std::uniform_real_distribution<double> x_distribution(-1.0, 1.0);
//  const auto x_coordinate = x_distribution(gen);
//
//  std::uniform_real_distribution<double> y_distribution(-1.0, 1.0);
//  const auto y_coordinate = y_distribution(gen);
//
//  return {x_coordinate, y_coordinate};
//}
//
// Quadrature_Rule Reference_Quadrilateral::make_quadrature_rule(
//    const ushort integrand_order) const {
//  switch (integrand_order) {
//    case 0:
//    case 1:
//      return {{{0, 0}}, {4}};
//    case 2:
//    case 3:
//      return {{{-0.577350269189626, -0.577350269189626},
//               {-0.577350269189626, 0.577350269189626},
//               {0.577350269189626, -0.577350269189626},
//               {0.577350269189626, 0.577350269189626}},
//              {1, 1, 1, 1}};
//    case 4:
//    case 5:
//      return {{{-0.774596669241483, -0.774596669241483},
//               {-0.774596669241483, 0},
//               {-0.774596669241483, 0.774596669241483},
//               {0, -0.774596669241483},
//               {0, 0},
//               {0, 0.774596669241483},
//               {0.774596669241483, -0.774596669241483},
//               {0.774596669241483, 0},
//               {0.774596669241483, 0.774596669241483}},
//              {0.3086419753086403, 0.4938271604938259, 0.3086419753086403,
//               0.4938271604938259, 0.7901234567901235, 0.4938271604938259,
//               0.3086419753086403, 0.4938271604938259, 0.3086419753086403}};
//    case 6:
//    case 7:
//      return {{{-0.861136311594052, -0.861136311594052},
//               {-0.861136311594052, -0.339981043584856},
//               {-0.861136311594052, 0.339981043584856},
//               {-0.861136311594052, 0.861136311594052},
//               {-0.339981043584856, -0.861136311594052},
//               {-0.339981043584856, -0.339981043584856},
//               {-0.339981043584856, 0.339981043584856},
//               {-0.339981043584856, 0.861136311594052},
//               {0.339981043584856, -0.861136311594052},
//               {0.339981043584856, -0.339981043584856},
//               {0.339981043584856, 0.339981043584856},
//               {0.339981043584856, 0.861136311594052},
//               {0.861136311594052, -0.861136311594052},
//               {0.861136311594052, -0.339981043584856},
//               {0.861136311594052, 0.339981043584856},
//               {0.861136311594052, 0.861136311594052}},
//              {0.1210029932856021, 0.2268518518518519, 0.2268518518518519,
//               0.1210029932856021, 0.2268518518518519, 0.4252933030106941,
//               0.4252933030106941, 0.2268518518518519, 0.2268518518518519,
//               0.4252933030106941, 0.4252933030106941, 0.2268518518518519,
//               0.1210029932856021, 0.2268518518518519, 0.2268518518518519,
//               0.1210029932856021}};
//    case 8:
//    case 9:
//      return {{{-0.906179845938664, -0.906179845938664},
//               {-0.906179845938664, -0.538469310105683},
//               {-0.906179845938664, 0},
//               {-0.906179845938664, 0.538469310105683},
//               {-0.906179845938664, 0.906179845938664},
//               {-0.538469310105683, -0.906179845938664},
//               {-0.538469310105683, -0.538469310105683},
//               {-0.538469310105683, 0},
//               {-0.538469310105683, 0.538469310105683},
//               {-0.538469310105683, 0.906179845938664},
//               {0, -0.906179845938664},
//               {0, -0.538469310105683},
//               {0, 0},
//               {0, 0.538469310105683},
//               {0, 0.906179845938664},
//               {0.538469310105683, -0.906179845938664},
//               {0.538469310105683, -0.538469310105683},
//               {0.538469310105683, 0},
//               {0.538469310105683, 0.538469310105683},
//               {0.538469310105683, 0.906179845938664},
//               {0.906179845938664, -0.906179845938664},
//               {0.906179845938664, -0.538469310105683},
//               {0.906179845938664, 0},
//               {0.906179845938664, 0.538469310105683},
//               {0.906179845938664, 0.906179845938664}},
//              {0.05613434886242859, 0.1133999999999998,  0.1347850723875209,
//               0.1133999999999998,  0.05613434886242859, 0.1133999999999998,
//               0.2290854042239907,  0.2722865325507505,  0.2290854042239907,
//               0.1133999999999998,  0.1347850723875209,  0.2722865325507505,
//               0.3236345679012347,  0.2722865325507505,  0.1347850723875209,
//               0.1133999999999998,  0.2290854042239907,  0.2722865325507505,
//               0.2290854042239907,  0.1133999999999998, 0.05613434886242859,
//               0.1133999999999998,  0.1347850723875209,  0.1133999999999998,
//               0.05613434886242859}};
//    case 10:
//    case 11:
//      return {{{-0.9324695142031521, -0.9324695142031521},
//               {-0.9324695142031521, -0.661209386466264},
//               {-0.9324695142031521, -0.238619186083197},
//               {-0.9324695142031521, 0.238619186083197},
//               {-0.9324695142031521, 0.661209386466264},
//               {-0.9324695142031521, 0.9324695142031521},
//               {-0.661209386466264, -0.9324695142031521},
//               {-0.661209386466264, -0.661209386466264},
//               {-0.661209386466264, -0.238619186083197},
//               {-0.661209386466264, 0.238619186083197},
//               {-0.661209386466264, 0.661209386466264},
//               {-0.661209386466264, 0.9324695142031521},
//               {-0.238619186083197, -0.9324695142031521},
//               {-0.238619186083197, -0.661209386466264},
//               {-0.238619186083197, -0.238619186083197},
//               {-0.238619186083197, 0.238619186083197},
//               {-0.238619186083197, 0.661209386466264},
//               {-0.238619186083197, 0.9324695142031521},
//               {0.238619186083197, -0.9324695142031521},
//               {0.238619186083197, -0.661209386466264},
//               {0.238619186083197, -0.238619186083197},
//               {0.238619186083197, 0.238619186083197},
//               {0.238619186083197, 0.661209386466264},
//               {0.238619186083197, 0.9324695142031521},
//               {0.661209386466264, -0.9324695142031521},
//               {0.661209386466264, -0.661209386466264},
//               {0.661209386466264, -0.238619186083197},
//               {0.661209386466264, 0.238619186083197},
//               {0.661209386466264, 0.661209386466264},
//               {0.661209386466264, 0.9324695142031521},
//               {0.9324695142031521, -0.9324695142031521},
//               {0.9324695142031521, -0.661209386466264},
//               {0.9324695142031521, -0.238619186083197},
//               {0.9324695142031521, 0.238619186083197},
//               {0.9324695142031521, 0.661209386466264},
//               {0.9324695142031521, 0.9324695142031521}},
//              {0.02935208168898062, 0.06180729337238363,
//              0.08016511731780691,
//               0.08016511731780691, 0.06180729337238363,
//               0.02935208168898062, 0.06180729337238363, 0.1301489125881677,
//               0.168805367087588, 0.168805367087588,   0.1301489125881677,
//               0.06180729337238363, 0.08016511731780691, 0.168805367087588,
//               0.2189434501672965, 0.2189434501672965,  0.168805367087588,
//               0.08016511731780691, 0.08016511731780691, 0.168805367087588,
//               0.2189434501672965, 0.2189434501672965,  0.168805367087588,
//               0.08016511731780691, 0.06180729337238363, 0.1301489125881677,
//               0.168805367087588, 0.168805367087588,   0.1301489125881677,
//               0.06180729337238363, 0.02935208168898062,
//               0.06180729337238363, 0.08016511731780691,
//               0.08016511731780691, 0.06180729337238363,
//               0.02935208168898062}};
//    case 12:
//    case 13:
//      return {{{-0.949107912342758, -0.949107912342758},
//               {-0.949107912342758, -0.741531185599394},
//               {-0.949107912342758, -0.405845151377397},
//               {-0.949107912342758, 0},
//               {-0.949107912342758, 0.405845151377397},
//               {-0.949107912342758, 0.741531185599394},
//               {-0.949107912342758, 0.949107912342758},
//               {-0.741531185599394, -0.949107912342758},
//               {-0.741531185599394, -0.741531185599394},
//               {-0.741531185599394, -0.405845151377397},
//               {-0.741531185599394, 0},
//               {-0.741531185599394, 0.405845151377397},
//               {-0.741531185599394, 0.741531185599394},
//               {-0.741531185599394, 0.949107912342758},
//               {-0.405845151377397, -0.949107912342758},
//               {-0.405845151377397, -0.741531185599394},
//               {-0.405845151377397, -0.405845151377397},
//               {-0.405845151377397, 0},
//               {-0.405845151377397, 0.405845151377397},
//               {-0.405845151377397, 0.741531185599394},
//               {-0.405845151377397, 0.949107912342758},
//               {0, -0.949107912342758},
//               {0, -0.741531185599394},
//               {0, -0.405845151377397},
//               {0, 0},
//               {0, 0.405845151377397},
//               {0, 0.741531185599394},
//               {0, 0.949107912342758},
//               {0.405845151377397, -0.949107912342758},
//               {0.405845151377397, -0.741531185599394},
//               {0.405845151377397, -0.405845151377397},
//               {0.405845151377397, 0},
//               {0.405845151377397, 0.405845151377397},
//               {0.405845151377397, 0.741531185599394},
//               {0.405845151377397, 0.949107912342758},
//               {0.741531185599394, -0.949107912342758},
//               {0.741531185599394, -0.741531185599394},
//               {0.741531185599394, -0.405845151377397},
//               {0.741531185599394, 0},
//               {0.741531185599394, 0.405845151377397},
//               {0.741531185599394, 0.741531185599394},
//               {0.741531185599394, 0.949107912342758},
//               {0.949107912342758, -0.949107912342758},
//               {0.949107912342758, -0.741531185599394},
//               {0.949107912342758, -0.405845151377397},
//               {0.949107912342758, 0},
//               {0.949107912342758, 0.405845151377397},
//               {0.949107912342758, 0.741531185599394},
//               {0.949107912342758, 0.949107912342758}},
//              {0.01676635646375341, 0.03621764315423957,
//              0.04944125117191326,
//               0.05411943075792766, 0.04944125117191326,
//               0.03621764315423957, 0.01676635646375341,
//               0.03621764315423957, 0.0782351060281697, 0.1067999237589047,
//               0.1169054370959263,  0.1067999237589047, 0.0782351060281697,
//               0.03621764315423957, 0.04944125117191326, 0.1067999237589047,
//               0.1457941874687417,  0.159589376211119, 0.1457941874687417,
//               0.1067999237589047,  0.04944125117191326,
//               0.05411943075792766, 0.1169054370959263,  0.159589376211119,
//               0.1746898792169926,  0.159589376211119,   0.1169054370959263,
//               0.05411943075792766, 0.04944125117191326, 0.1067999237589047,
//               0.1457941874687417,  0.159589376211119,   0.1457941874687417,
//               0.1067999237589047,  0.04944125117191326,
//               0.03621764315423957, 0.0782351060281697,  0.1067999237589047,
//               0.1169054370959263, 0.1067999237589047,  0.0782351060281697,
//               0.03621764315423957, 0.01676635646375341,
//               0.03621764315423957, 0.04944125117191326,
//               0.05411943075792766, 0.04944125117191326,
//               0.03621764315423957, 0.01676635646375341}};
//    case 14:
//    case 15:
//      return {{{-0.960289856497536, -0.960289856497536},
//               {-0.960289856497536, -0.796666477413627},
//               {-0.960289856497536, -0.525532409916329},
//               {-0.960289856497536, -0.18343464249565},
//               {-0.960289856497536, 0.18343464249565},
//               {-0.960289856497536, 0.525532409916329},
//               {-0.960289856497536, 0.796666477413627},
//               {-0.960289856497536, 0.960289856497536},
//               {-0.796666477413627, -0.960289856497536},
//               {-0.796666477413627, -0.796666477413627},
//               {-0.796666477413627, -0.525532409916329},
//               {-0.796666477413627, -0.18343464249565},
//               {-0.796666477413627, 0.18343464249565},
//               {-0.796666477413627, 0.525532409916329},
//               {-0.796666477413627, 0.796666477413627},
//               {-0.796666477413627, 0.960289856497536},
//               {-0.525532409916329, -0.960289856497536},
//               {-0.525532409916329, -0.796666477413627},
//               {-0.525532409916329, -0.525532409916329},
//               {-0.525532409916329, -0.18343464249565},
//               {-0.525532409916329, 0.18343464249565},
//               {-0.525532409916329, 0.525532409916329},
//               {-0.525532409916329, 0.796666477413627},
//               {-0.525532409916329, 0.960289856497536},
//               {-0.18343464249565, -0.960289856497536},
//               {-0.18343464249565, -0.796666477413627},
//               {-0.18343464249565, -0.525532409916329},
//               {-0.18343464249565, -0.18343464249565},
//               {-0.18343464249565, 0.18343464249565},
//               {-0.18343464249565, 0.525532409916329},
//               {-0.18343464249565, 0.796666477413627},
//               {-0.18343464249565, 0.960289856497536},
//               {0.18343464249565, -0.960289856497536},
//               {0.18343464249565, -0.796666477413627},
//               {0.18343464249565, -0.525532409916329},
//               {0.18343464249565, -0.18343464249565},
//               {0.18343464249565, 0.18343464249565},
//               {0.18343464249565, 0.525532409916329},
//               {0.18343464249565, 0.796666477413627},
//               {0.18343464249565, 0.960289856497536},
//               {0.525532409916329, -0.960289856497536},
//               {0.525532409916329, -0.796666477413627},
//               {0.525532409916329, -0.525532409916329},
//               {0.525532409916329, -0.18343464249565},
//               {0.525532409916329, 0.18343464249565},
//               {0.525532409916329, 0.525532409916329},
//               {0.525532409916329, 0.796666477413627},
//               {0.525532409916329, 0.960289856497536},
//               {0.796666477413627, -0.960289856497536},
//               {0.796666477413627, -0.796666477413627},
//               {0.796666477413627, -0.525532409916329},
//               {0.796666477413627, -0.18343464249565},
//               {0.796666477413627, 0.18343464249565},
//               {0.796666477413627, 0.525532409916329},
//               {0.796666477413627, 0.796666477413627},
//               {0.796666477413627, 0.960289856497536},
//               {0.960289856497536, -0.960289856497536},
//               {0.960289856497536, -0.796666477413627},
//               {0.960289856497536, -0.525532409916329},
//               {0.960289856497536, -0.18343464249565},
//               {0.960289856497536, 0.18343464249565},
//               {0.960289856497536, 0.525532409916329},
//               {0.960289856497536, 0.796666477413627},
//               {0.960289856497536, 0.960289856497536}},
//              {0.01024721655949217, 0.02251130661645495,
//              0.03175606458678213,
//               0.03671394852764775, 0.03671394852764775,
//               0.03175606458678213, 0.02251130661645495,
//               0.01024721655949217, 0.02251130661645495,
//               0.04945332448455272, 0.06976240842522279,
//               0.08065399492714355, 0.08065399492714355,
//               0.06976240842522279, 0.04945332448455272,
//               0.02251130661645495, 0.03175606458678213,
//               0.06976240842522279, 0.09841185966795399, 0.1137763131979281,
//               0.1137763131979281, 0.09841185966795399, 0.06976240842522279,
//               0.03175606458678213, 0.03671394852764775,
//               0.08065399492714355, 0.1137763131979281, 0.1315395267256426,
//               0.1315395267256426,  0.1137763131979281, 0.08065399492714355,
//               0.03671394852764775, 0.03671394852764775,
//               0.08065399492714355, 0.1137763131979281,  0.1315395267256426,
//               0.1315395267256426,  0.1137763131979281, 0.08065399492714355,
//               0.03671394852764775, 0.03175606458678213,
//               0.06976240842522279, 0.09841185966795399, 0.1137763131979281,
//               0.1137763131979281, 0.09841185966795399, 0.06976240842522279,
//               0.03175606458678213, 0.02251130661645495,
//               0.04945332448455272, 0.06976240842522279,
//               0.08065399492714355, 0.08065399492714355,
//               0.06976240842522279, 0.04945332448455272,
//               0.02251130661645495, 0.01024721655949217,
//               0.02251130661645495, 0.03175606458678213,
//               0.03671394852764775, 0.03671394852764775,
//               0.03175606458678213, 0.02251130661645495,
//               0.01024721655949217}};
//    case 16:
//    case 17:
//      return {
//          {{-0.968160239507626, -0.968160239507626},
//           {-0.968160239507626, -0.836031107326636},
//           {-0.968160239507626, -0.61337143270059},
//           {-0.968160239507626, -0.324253423403809},
//           {-0.968160239507626, 0},
//           {-0.968160239507626, 0.324253423403809},
//           {-0.968160239507626, 0.61337143270059},
//           {-0.968160239507626, 0.836031107326636},
//           {-0.968160239507626, 0.968160239507626},
//           {-0.836031107326636, -0.968160239507626},
//           {-0.836031107326636, -0.836031107326636},
//           {-0.836031107326636, -0.61337143270059},
//           {-0.836031107326636, -0.324253423403809},
//           {-0.836031107326636, 0},
//           {-0.836031107326636, 0.324253423403809},
//           {-0.836031107326636, 0.61337143270059},
//           {-0.836031107326636, 0.836031107326636},
//           {-0.836031107326636, 0.968160239507626},
//           {-0.61337143270059, -0.968160239507626},
//           {-0.61337143270059, -0.836031107326636},
//           {-0.61337143270059, -0.61337143270059},
//           {-0.61337143270059, -0.324253423403809},
//           {-0.61337143270059, 0},
//           {-0.61337143270059, 0.324253423403809},
//           {-0.61337143270059, 0.61337143270059},
//           {-0.61337143270059, 0.836031107326636},
//           {-0.61337143270059, 0.968160239507626},
//           {-0.324253423403809, -0.968160239507626},
//           {-0.324253423403809, -0.836031107326636},
//           {-0.324253423403809, -0.61337143270059},
//           {-0.324253423403809, -0.324253423403809},
//           {-0.324253423403809, 0},
//           {-0.324253423403809, 0.324253423403809},
//           {-0.324253423403809, 0.61337143270059},
//           {-0.324253423403809, 0.836031107326636},
//           {-0.324253423403809, 0.968160239507626},
//           {0, -0.968160239507626},
//           {0, -0.836031107326636},
//           {0, -0.61337143270059},
//           {0, -0.324253423403809},
//           {0, 0},
//           {0, 0.324253423403809},
//           {0, 0.61337143270059},
//           {0, 0.836031107326636},
//           {0, 0.968160239507626},
//           {0.324253423403809, -0.968160239507626},
//           {0.324253423403809, -0.836031107326636},
//           {0.324253423403809, -0.61337143270059},
//           {0.324253423403809, -0.324253423403809},
//           {0.324253423403809, 0},
//           {0.324253423403809, 0.324253423403809},
//           {0.324253423403809, 0.61337143270059},
//           {0.324253423403809, 0.836031107326636},
//           {0.324253423403809, 0.968160239507626},
//           {0.61337143270059, -0.968160239507626},
//           {0.61337143270059, -0.836031107326636},
//           {0.61337143270059, -0.61337143270059},
//           {0.61337143270059, -0.324253423403809},
//           {0.61337143270059, 0},
//           {0.61337143270059, 0.324253423403809},
//           {0.61337143270059, 0.61337143270059},
//           {0.61337143270059, 0.836031107326636},
//           {0.61337143270059, 0.968160239507626},
//           {0.836031107326636, -0.968160239507626},
//           {0.836031107326636, -0.836031107326636},
//           {0.836031107326636, -0.61337143270059},
//           {0.836031107326636, -0.324253423403809},
//           {0.836031107326636, 0},
//           {0.836031107326636, 0.324253423403809},
//           {0.836031107326636, 0.61337143270059},
//           {0.836031107326636, 0.836031107326636},
//           {0.836031107326636, 0.968160239507626},
//           {0.968160239507626, -0.968160239507626},
//           {0.968160239507626, -0.836031107326636},
//           {0.968160239507626, -0.61337143270059},
//           {0.968160239507626, -0.324253423403809},
//           {0.968160239507626, 0},
//           {0.968160239507626, 0.324253423403809},
//           {0.968160239507626, 0.61337143270059},
//           {0.968160239507626, 0.836031107326636},
//           {0.968160239507626, 0.968160239507626}},
//          {0.006605526203548117, 0.01468206876911802, 0.02118097495063273,
//           0.02538581764295198,  0.02684000159064844, 0.02538581764295198,
//           0.02118097495063273,  0.01468206876911802, 0.006605526203548117,
//           0.01468206876911802,  0.03263375796243488, 0.04707884296259617,
//           0.05642492496569134,  0.05965713207003355, 0.05642492496569134,
//           0.04707884296259617,  0.03263375796243488, 0.01468206876911802,
//           0.02118097495063273,  0.04707884296259617, 0.06791793507962328,
//           0.08140098926681667,  0.08606390828653478, 0.08140098926681667,
//           0.06791793507962328,  0.04707884296259617, 0.02118097495063273,
//           0.02538581764295198,  0.05642492496569134, 0.08140098926681667,
//           0.09756069653543355,  0.1031492972582194,  0.09756069653543355,
//           0.08140098926681667,  0.05642492496569134, 0.02538581764295198,
//           0.02684000159064844,  0.05965713207003355, 0.08606390828653478,
//           0.1031492972582194,   0.1090580315916482,  0.1031492972582194,
//           0.08606390828653478,  0.05965713207003355, 0.02684000159064844,
//           0.02538581764295198,  0.05642492496569134, 0.08140098926681667,
//           0.09756069653543355,  0.1031492972582194,  0.09756069653543355,
//           0.08140098926681667,  0.05642492496569134, 0.02538581764295198,
//           0.02118097495063273,  0.04707884296259617, 0.06791793507962328,
//           0.08140098926681667,  0.08606390828653478, 0.08140098926681667,
//           0.06791793507962328,  0.04707884296259617, 0.02118097495063273,
//           0.01468206876911802,  0.03263375796243488, 0.04707884296259617,
//           0.05642492496569134,  0.05965713207003355, 0.05642492496569134,
//           0.04707884296259617,  0.03263375796243488, 0.01468206876911802,
//           0.006605526203548117, 0.01468206876911802, 0.02118097495063273,
//           0.02538581764295198,  0.02684000159064844, 0.02538581764295198,
//           0.02118097495063273,  0.01468206876911802,
//           0.006605526203548117}};
//    case 18:
//    case 19:
//      return {{{-0.973906528517172, -0.973906528517172},
//               {-0.973906528517172, -0.865063366688985},
//               {-0.973906528517172, -0.679409568299024},
//               {-0.973906528517172, -0.433395394129247},
//               {-0.973906528517172, -0.148874338981631},
//               {-0.973906528517172, 0.148874338981631},
//               {-0.973906528517172, 0.433395394129247},
//               {-0.973906528517172, 0.679409568299024},
//               {-0.973906528517172, 0.865063366688985},
//               {-0.973906528517172, 0.973906528517172},
//               {-0.865063366688985, -0.973906528517172},
//               {-0.865063366688985, -0.865063366688985},
//               {-0.865063366688985, -0.679409568299024},
//               {-0.865063366688985, -0.433395394129247},
//               {-0.865063366688985, -0.148874338981631},
//               {-0.865063366688985, 0.148874338981631},
//               {-0.865063366688985, 0.433395394129247},
//               {-0.865063366688985, 0.679409568299024},
//               {-0.865063366688985, 0.865063366688985},
//               {-0.865063366688985, 0.973906528517172},
//               {-0.679409568299024, -0.973906528517172},
//               {-0.679409568299024, -0.865063366688985},
//               {-0.679409568299024, -0.679409568299024},
//               {-0.679409568299024, -0.433395394129247},
//               {-0.679409568299024, -0.148874338981631},
//               {-0.679409568299024, 0.148874338981631},
//               {-0.679409568299024, 0.433395394129247},
//               {-0.679409568299024, 0.679409568299024},
//               {-0.679409568299024, 0.865063366688985},
//               {-0.679409568299024, 0.973906528517172},
//               {-0.433395394129247, -0.973906528517172},
//               {-0.433395394129247, -0.865063366688985},
//               {-0.433395394129247, -0.679409568299024},
//               {-0.433395394129247, -0.433395394129247},
//               {-0.433395394129247, -0.148874338981631},
//               {-0.433395394129247, 0.148874338981631},
//               {-0.433395394129247, 0.433395394129247},
//               {-0.433395394129247, 0.679409568299024},
//               {-0.433395394129247, 0.865063366688985},
//               {-0.433395394129247, 0.973906528517172},
//               {-0.148874338981631, -0.973906528517172},
//               {-0.148874338981631, -0.865063366688985},
//               {-0.148874338981631, -0.679409568299024},
//               {-0.148874338981631, -0.433395394129247},
//               {-0.148874338981631, -0.148874338981631},
//               {-0.148874338981631, 0.148874338981631},
//               {-0.148874338981631, 0.433395394129247},
//               {-0.148874338981631, 0.679409568299024},
//               {-0.148874338981631, 0.865063366688985},
//               {-0.148874338981631, 0.973906528517172},
//               {0.148874338981631, -0.973906528517172},
//               {0.148874338981631, -0.865063366688985},
//               {0.148874338981631, -0.679409568299024},
//               {0.148874338981631, -0.433395394129247},
//               {0.148874338981631, -0.148874338981631},
//               {0.148874338981631, 0.148874338981631},
//               {0.148874338981631, 0.433395394129247},
//               {0.148874338981631, 0.679409568299024},
//               {0.148874338981631, 0.865063366688985},
//               {0.148874338981631, 0.973906528517172},
//               {0.433395394129247, -0.973906528517172},
//               {0.433395394129247, -0.865063366688985},
//               {0.433395394129247, -0.679409568299024},
//               {0.433395394129247, -0.433395394129247},
//               {0.433395394129247, -0.148874338981631},
//               {0.433395394129247, 0.148874338981631},
//               {0.433395394129247, 0.433395394129247},
//               {0.433395394129247, 0.679409568299024},
//               {0.433395394129247, 0.865063366688985},
//               {0.433395394129247, 0.973906528517172},
//               {0.679409568299024, -0.973906528517172},
//               {0.679409568299024, -0.865063366688985},
//               {0.679409568299024, -0.679409568299024},
//               {0.679409568299024, -0.433395394129247},
//               {0.679409568299024, -0.148874338981631},
//               {0.679409568299024, 0.148874338981631},
//               {0.679409568299024, 0.433395394129247},
//               {0.679409568299024, 0.679409568299024},
//               {0.679409568299024, 0.865063366688985},
//               {0.679409568299024, 0.973906528517172},
//               {0.865063366688985, -0.973906528517172},
//               {0.865063366688985, -0.865063366688985},
//               {0.865063366688985, -0.679409568299024},
//               {0.865063366688985, -0.433395394129247},
//               {0.865063366688985, -0.148874338981631},
//               {0.865063366688985, 0.148874338981631},
//               {0.865063366688985, 0.433395394129247},
//               {0.865063366688985, 0.679409568299024},
//               {0.865063366688985, 0.865063366688985},
//               {0.865063366688985, 0.973906528517172},
//               {0.973906528517172, -0.973906528517172},
//               {0.973906528517172, -0.865063366688985},
//               {0.973906528517172, -0.679409568299024},
//               {0.973906528517172, -0.433395394129247},
//               {0.973906528517172, -0.148874338981631},
//               {0.973906528517172, 0.148874338981631},
//               {0.973906528517172, 0.433395394129247},
//               {0.973906528517172, 0.679409568299024},
//               {0.973906528517172, 0.865063366688985},
//               {0.973906528517172, 0.973906528517172}},
//              {0.004445068151927624, 0.009964122356616333,
//              0.01460678230864107,
//               0.01795237415398759,  0.01970299733751538,
//               0.01970299733751538, 0.01795237415398759,
//               0.01460678230864107, 0.009964122356616333,
//               0.004445068151927624, 0.009964122356616333,
//               0.02233570576292887, 0.03274275245850679,
//               0.04024227448222971,  0.04416649409029931,
//               0.04416649409029931, 0.04024227448222971,
//               0.03274275245850679, 0.02233570576292887,
//               0.009964122356616333, 0.01460678230864107,
//               0.03274275245850679, 0.04799883424048429,
//               0.05899266608023896, 0.0647453274281109, 0.0647453274281109,
//               0.05899266608023896, 0.04799883424048429,
//               0.03274275245850679,  0.01460678230864107,
//               0.01795237415398759, 0.04024227448222971,
//               0.05899266608023896, 0.07250456612796818,
//               0.07957483846557158,  0.07957483846557158,
//               0.07250456612796818, 0.05899266608023896,
//               0.04024227448222971, 0.01795237415398759,
//               0.01970299733751538,  0.04416649409029931,
//               0.0647453274281109, 0.07957483846557158, 0.08733456739325582,
//               0.08733456739325582, 0.07957483846557158, 0.0647453274281109,
//               0.04416649409029931, 0.01970299733751538,
//               0.01970299733751538, 0.04416649409029931, 0.0647453274281109,
//               0.07957483846557158, 0.08733456739325582,
//               0.08733456739325582,  0.07957483846557158,
//               0.0647453274281109, 0.04416649409029931, 0.01970299733751538,
//               0.01795237415398759, 0.04024227448222971,
//               0.05899266608023896, 0.07250456612796818,
//               0.07957483846557158,  0.07957483846557158,
//               0.07250456612796818, 0.05899266608023896,
//               0.04024227448222971, 0.01795237415398759,
//               0.01460678230864107,  0.03274275245850679,
//               0.04799883424048429, 0.05899266608023896, 0.0647453274281109,
//               0.0647453274281109, 0.05899266608023896, 0.04799883424048429,
//               0.03274275245850679, 0.01460678230864107,
//               0.009964122356616333, 0.02233570576292887,
//               0.03274275245850679, 0.04024227448222971,
//               0.04416649409029931,  0.04416649409029931,
//               0.04024227448222971, 0.03274275245850679,
//               0.02233570576292887, 0.009964122356616333,
//               0.004445068151927624, 0.009964122356616333,
//               0.01460678230864107, 0.01795237415398759,
//               0.01970299733751538,  0.01970299733751538,
//               0.01795237415398759, 0.01460678230864107,
//               0.009964122356616333, 0.004445068151927624}};
//    case 20:
//    case 21:
//      return {{{-0.978228658146057, -0.978228658146057},
//               {-0.978228658146057, -0.887062599768095},
//               {-0.978228658146057, -0.730152005574049},
//               {-0.978228658146057, -0.519096129206812},
//               {-0.978228658146057, -0.269543155952345},
//               {-0.978228658146057, 0},
//               {-0.978228658146057, 0.269543155952345},
//               {-0.978228658146057, 0.519096129206812},
//               {-0.978228658146057, 0.730152005574049},
//               {-0.978228658146057, 0.887062599768095},
//               {-0.978228658146057, 0.978228658146057},
//               {-0.887062599768095, -0.978228658146057},
//               {-0.887062599768095, -0.887062599768095},
//               {-0.887062599768095, -0.730152005574049},
//               {-0.887062599768095, -0.519096129206812},
//               {-0.887062599768095, -0.269543155952345},
//               {-0.887062599768095, 0},
//               {-0.887062599768095, 0.269543155952345},
//               {-0.887062599768095, 0.519096129206812},
//               {-0.887062599768095, 0.730152005574049},
//               {-0.887062599768095, 0.887062599768095},
//               {-0.887062599768095, 0.978228658146057},
//               {-0.730152005574049, -0.978228658146057},
//               {-0.730152005574049, -0.887062599768095},
//               {-0.730152005574049, -0.730152005574049},
//               {-0.730152005574049, -0.519096129206812},
//               {-0.730152005574049, -0.269543155952345},
//               {-0.730152005574049, 0},
//               {-0.730152005574049, 0.269543155952345},
//               {-0.730152005574049, 0.519096129206812},
//               {-0.730152005574049, 0.730152005574049},
//               {-0.730152005574049, 0.887062599768095},
//               {-0.730152005574049, 0.978228658146057},
//               {-0.519096129206812, -0.978228658146057},
//               {-0.519096129206812, -0.887062599768095},
//               {-0.519096129206812, -0.730152005574049},
//               {-0.519096129206812, -0.519096129206812},
//               {-0.519096129206812, -0.269543155952345},
//               {-0.519096129206812, 0},
//               {-0.519096129206812, 0.269543155952345},
//               {-0.519096129206812, 0.519096129206812},
//               {-0.519096129206812, 0.730152005574049},
//               {-0.519096129206812, 0.887062599768095},
//               {-0.519096129206812, 0.978228658146057},
//               {-0.269543155952345, -0.978228658146057},
//               {-0.269543155952345, -0.887062599768095},
//               {-0.269543155952345, -0.730152005574049},
//               {-0.269543155952345, -0.519096129206812},
//               {-0.269543155952345, -0.269543155952345},
//               {-0.269543155952345, 0},
//               {-0.269543155952345, 0.269543155952345},
//               {-0.269543155952345, 0.519096129206812},
//               {-0.269543155952345, 0.730152005574049},
//               {-0.269543155952345, 0.887062599768095},
//               {-0.269543155952345, 0.978228658146057},
//               {0, -0.978228658146057},
//               {0, -0.887062599768095},
//               {0, -0.730152005574049},
//               {0, -0.519096129206812},
//               {0, -0.269543155952345},
//               {0, 0},
//               {0, 0.269543155952345},
//               {0, 0.519096129206812},
//               {0, 0.730152005574049},
//               {0, 0.887062599768095},
//               {0, 0.978228658146057},
//               {0.269543155952345, -0.978228658146057},
//               {0.269543155952345, -0.887062599768095},
//               {0.269543155952345, -0.730152005574049},
//               {0.269543155952345, -0.519096129206812},
//               {0.269543155952345, -0.269543155952345},
//               {0.269543155952345, 0},
//               {0.269543155952345, 0.269543155952345},
//               {0.269543155952345, 0.519096129206812},
//               {0.269543155952345, 0.730152005574049},
//               {0.269543155952345, 0.887062599768095},
//               {0.269543155952345, 0.978228658146057},
//               {0.519096129206812, -0.978228658146057},
//               {0.519096129206812, -0.887062599768095},
//               {0.519096129206812, -0.730152005574049},
//               {0.519096129206812, -0.519096129206812},
//               {0.519096129206812, -0.269543155952345},
//               {0.519096129206812, 0},
//               {0.519096129206812, 0.269543155952345},
//               {0.519096129206812, 0.519096129206812},
//               {0.519096129206812, 0.730152005574049},
//               {0.519096129206812, 0.887062599768095},
//               {0.519096129206812, 0.978228658146057},
//               {0.730152005574049, -0.978228658146057},
//               {0.730152005574049, -0.887062599768095},
//               {0.730152005574049, -0.730152005574049},
//               {0.730152005574049, -0.519096129206812},
//               {0.730152005574049, -0.269543155952345},
//               {0.730152005574049, 0},
//               {0.730152005574049, 0.269543155952345},
//               {0.730152005574049, 0.519096129206812},
//               {0.730152005574049, 0.730152005574049},
//               {0.730152005574049, 0.887062599768095},
//               {0.730152005574049, 0.978228658146057},
//               {0.887062599768095, -0.978228658146057},
//               {0.887062599768095, -0.887062599768095},
//               {0.887062599768095, -0.730152005574049},
//               {0.887062599768095, -0.519096129206812},
//               {0.887062599768095, -0.269543155952345},
//               {0.887062599768095, 0},
//               {0.887062599768095, 0.269543155952345},
//               {0.887062599768095, 0.519096129206812},
//               {0.887062599768095, 0.730152005574049},
//               {0.887062599768095, 0.887062599768095},
//               {0.887062599768095, 0.978228658146057},
//               {0.978228658146057, -0.978228658146057},
//               {0.978228658146057, -0.887062599768095},
//               {0.978228658146057, -0.730152005574049},
//               {0.978228658146057, -0.519096129206812},
//               {0.978228658146057, -0.269543155952345},
//               {0.978228658146057, 0},
//               {0.978228658146057, 0.269543155952345},
//               {0.978228658146057, 0.519096129206812},
//               {0.978228658146057, 0.730152005574049},
//               {0.978228658146057, 0.887062599768095},
//               {0.978228658146057, 0.978228658146057}},
//              {0.00309898936476797,  0.006990879226030937,
//              0.01037050911011677,
//               0.01298156273526248,  0.01462995242450422,
//               0.0151933485109832, 0.01462995242450422, 0.01298156273526248,
//               0.01037050911011677, 0.006990879226030937,
//               0.00309898936476797, 0.006990879226030937,
//               0.01577042919494179, 0.02339439351599974,
//               0.02928455911437395, 0.03300309179665262,
//               0.03427403323380979, 0.03300309179665262,
//               0.02928455911437395,  0.02339439351599974,
//               0.01577042919494179, 0.006990879226030937,
//               0.01037050911011677, 0.02339439351599974,
//               0.03470404268749963,  0.04344171559287417,
//               0.04895791402958097, 0.05084327198332528,
//               0.04895791402958097, 0.04344171559287417,
//               0.03470404268749963,  0.02339439351599974,
//               0.01037050911011677, 0.01298156273526248,
//               0.02928455911437395, 0.04344171559287417,
//               0.05437933184458445,  0.0612843810862277, 0.0636444284373343,
//               0.0612843810862277,   0.05437933184458445,
//               0.04344171559287417, 0.02928455911437395,
//               0.01298156273526248, 0.01462995242450422,
//               0.03300309179665262,  0.04895791402958097,
//               0.0612843810862277, 0.0690662286152384,   0.0717259531160859,
//               0.0690662286152384, 0.0612843810862277, 0.04895791402958097,
//               0.03300309179665262, 0.01462995242450422, 0.0151933485109832,
//               0.03427403323380979, 0.05084327198332528, 0.0636444284373343,
//               0.0717259531160859, 0.07448810299272479,  0.0717259531160859,
//               0.0636444284373343, 0.05084327198332528, 0.03427403323380979,
//               0.0151933485109832, 0.01462995242450422, 0.03300309179665262,
//               0.04895791402958097, 0.0612843810862277, 0.0690662286152384,
//               0.0717259531160859, 0.0690662286152384,   0.0612843810862277,
//               0.04895791402958097, 0.03300309179665262,
//               0.01462995242450422, 0.01298156273526248,
//               0.02928455911437395,  0.04344171559287417,
//               0.05437933184458445, 0.0612843810862277, 0.0636444284373343,
//               0.0612843810862277, 0.05437933184458445, 0.04344171559287417,
//               0.02928455911437395, 0.01298156273526248,
//               0.01037050911011677, 0.02339439351599974,
//               0.03470404268749963,  0.04344171559287417,
//               0.04895791402958097, 0.05084327198332528,
//               0.04895791402958097, 0.04344171559287417,
//               0.03470404268749963,  0.02339439351599974,
//               0.01037050911011677, 0.006990879226030937,
//               0.01577042919494179, 0.02339439351599974,
//               0.02928455911437395,  0.03300309179665262,
//               0.03427403323380979, 0.03300309179665262,
//               0.02928455911437395, 0.02339439351599974,
//               0.01577042919494179, 0.006990879226030937,
//               0.00309898936476797, 0.006990879226030937,
//               0.01037050911011677, 0.01298156273526248,
//               0.01462995242450422,  0.0151933485109832,
//               0.01462995242450422, 0.01298156273526248,
//               0.01037050911011677, 0.006990879226030937,
//               0.00309898936476797}};
//    default:
//      EXCEPTION("not supported order");
//      return {};
//  }
//}
// bool Reference_Quadrilateral::is_simplex(void) const { return false; };
//
// bool Reference_Quadrilateral::is_line(void) const { return false; }
//
// std::vector<std::vector<ushort>>
// Reference_Quadrilateral::set_of_face_vertex_index_sequences(void) const {
//  //      2
//  //   3����������2
//  // 3  ��     ��   1
//  //   0����������1
//  //      0
//  std::vector<ushort> face0_node_index = {0, 1};
//  std::vector<ushort> face1_node_index = {1, 2};
//  std::vector<ushort> face2_node_index = {2, 3};
//  std::vector<ushort> face3_node_index = {3, 0};
//  return {face0_node_index, face1_node_index, face2_node_index,
//          face3_node_index};
//};
//
// std::vector<std::vector<ushort>>
// Reference_Quadrilateral::set_of_face_node_index_sequences(void) const {
//  //      2
//  //   3����������2
//  // 3  ��     ��   1
//  //   0����������1
//  //      0
//  constexpr ushort num_face = 4;
//  std::vector<std::vector<ushort>> set_of_face_node_index_orders(num_face);
//  set_of_face_node_index_orders[0] = {0, 1};
//  set_of_face_node_index_orders[1] = {1, 2};
//  set_of_face_node_index_orders[2] = {2, 3};
//  set_of_face_node_index_orders[3] = {3, 0};
//
//  if (this->order_ > 1) {
//    const ushort num_additional_point = this->order_ - 1;
//
//    ushort index = num_face;
//    for (ushort iface = 0; iface < num_face; ++iface)
//      for (ushort ipoint = 0; ipoint < num_additional_point; ++ipoint)
//        set_of_face_node_index_orders[iface].push_back(index++);
//  }
//
//  return set_of_face_node_index_orders;
//};
//
// std::vector<std::shared_ptr<const Reference_Geometry>>
// Reference_Quadrilateral::sub_simplex_reference_geometries(void) const {
//  constexpr auto num_sub_simplex = 4;
//
//  std::vector<std::shared_ptr<const Reference_Geometry>>
//      sub_simplex_reference_geometries(num_sub_simplex);
//
//  for (ushort i = 0; i < num_sub_simplex; ++i) {
//    sub_simplex_reference_geometries[i] =
//        Reference_Geometry_Container::get_shared_ptr(Figure::triangle,
//                                                     this->order_);
//  }
//
//  return sub_simplex_reference_geometries;
//}
//
// std::vector<std::vector<ushort>>
// Reference_Quadrilateral::set_of_sub_simplex_vertex_index_sequences(void)
// const {
//  //   3����������2
//  //   ��     ��
//  //   0����������1
//
//  constexpr auto num_sub_simplex = 4;
//
//  const std::vector<ushort> simplex1 = {0, 1, 3};
//  const std::vector<ushort> simplex2 = {1, 2, 0};
//  const std::vector<ushort> simplex3 = {2, 3, 1};
//  const std::vector<ushort> simplex4 = {3, 0, 2};
//
//  return {simplex1, simplex2, simplex3, simplex4};
//}
//
// Irrational_Function Reference_Quadrilateral::scale_function(
//    const Vector_Function<Polynomial>& mapping_function) const {
//  constexpr ushort r = 0;
//  constexpr ushort s = 1;
//  const auto mf_r = mapping_function.get_differentiate(r);
//  const auto mf_s = mapping_function.get_differentiate(s);
//  const auto cross_product = mf_r.cross_product(mf_s);
//  return cross_product.L2_norm();
//};
//
// ushort Reference_Quadrilateral::scale_function_order(void) const {
//  REQUIRE(this->order_ == 1, "high order mesh is not supported yet");
//  return 1;
//}
//
// std::vector<Euclidean_Vector> Reference_Quadrilateral::make_mapping_points(
//    void) const {
//  //   3����������2
//  //   ��     ��
//  //   0����������1
//  switch (this->order_) {
//    case 1:
//      return {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
//    default:
//      EXCEPTION("unsuported figure order");
//      return {};
//  }
//}
//
// std::vector<Euclidean_Vector> Reference_Quadrilateral::make_post_points(
//    const ushort post_order) const {
//  const auto n = post_order;
//  const auto num_post_nodes = (n + 2) * (n + 2);
//
//  std::vector<Euclidean_Vector> post_points;
//  post_points.reserve(num_post_nodes);
//
//  const auto num_division = n + 1;
//  const auto delta = 2.0 / num_division;
//
//  const auto X0_start_coord = -1.0;
//  const auto X1_start_coord = -1.0;
//
//  for (ushort j = 0; j <= n + 1; ++j) {
//    for (ushort i = 0; i <= n + 1; ++i) {
//      const double X0_coord = X0_start_coord + delta * i;
//      const double X1_coord = X1_start_coord + delta * j;
//      post_points.push_back({X0_coord, X1_coord});
//    }
//  }
//
//  return post_points;
//}
//
// std::vector<std::vector<uint>>
// Reference_Quadrilateral::make_connectivities(
//    const ushort post_order) const {
//  const auto n = post_order;
//  const auto num_connectivity = 2 * (n + 1) * (n + 1);
//
//  std::vector<std::vector<uint>> connectivities;
//  connectivities.reserve(num_connectivity);
//
//  for (ushort j = 0; j <= n; j++) {
//    for (ushort i = 0; i <= n; i++) {
//      //   I2��������I2+1
//      //   ��     ��
//      //   I1��������I1+1
//
//      const uint I1 = (n + 2) * j + i;
//      const uint I2 = (n + 2) * (j + 1) + i;
//
//      auto quadrilateral_connectivityies =
//          this->quadrilateral_connectivities({I1, I1 + 1, I2, I2 + 1});
//
//      for (auto& connectivity : quadrilateral_connectivityies)
//        connectivities.push_back(std::move(connectivity));
//    }
//  }
//
//  return connectivities;
//}
//
// Vector_Function<Polynomial>
// Reference_Quadrilateral::make_mapping_monomial_vector_function(void) const
// {
//  Polynomial r("x0");
//  Polynomial s("x1");
//
//  const auto n = this->order_;
//  const auto num_monomial = (n + 1) * (n + 1);
//  ms::sym::Polynomials mapping_monomial_vector(num_monomial);
//
//  for (ushort a = 0, index = 0; a <= n; ++a) {
//    for (ushort b = 0; b <= a; ++b) {
//      mapping_monomial_vector[index++] = (r ^ a) * (s ^ b);
//    }
//
//    if (a == 0) {
//      continue;
//    }
//
//    for (int c = static_cast<int>(a - 1); 0 <= c; --c) {
//      mapping_monomial_vector[index++] = (r ^ c) * (s ^ a);
//    }
//  }
//
//  return mapping_monomial_vector;  // 1 r rs s r^2 r^2s r^2s^2 rs^2 s^2...
//}
//
// const std::shared_ptr<const Reference_Geometry>&
// Reference_Geometry_Container::get_shared_ptr(const Figure figure,
//                                             const ushort order) {
//  switch (figure) {
//    case Figure::line: {
//      if (!Reference_Geometry_Container::order_to_reference_line_.contains(
//              order)) {
//        Reference_Geometry_Container::store_line(order);
//      }
//
//      return
//      Reference_Geometry_Container::order_to_reference_line_.at(order);
//    }
//    case Figure::triangle: {
//      if
//      (!Reference_Geometry_Container::order_to_reference_triangle_.contains(
//              order)) {
//        Reference_Geometry_Container::store_triangle(order);
//      }
//
//      return Reference_Geometry_Container::order_to_reference_triangle_.at(
//          order);
//    }
//    case Figure::quadrilateral: {
//      if (!Reference_Geometry_Container::order_to_reference_quadrilateral_
//               .contains(order)) {
//        Reference_Geometry_Container::store_quadrilateral(order);
//      }
//
//      return
//      Reference_Geometry_Container::order_to_reference_quadrilateral_.at(
//          order);
//    }
//    default: {
//      EXCEPTION("not supproted figure");
//      return Reference_Geometry_Container::error_;
//    }
//  }
//}
//
// void Reference_Geometry_Container::store_line(const ushort order) {
//  Reference_Geometry_Container::order_to_reference_line_.emplace(
//      order, std::make_shared<Reference_Line>(order));
//}
//
// void Reference_Geometry_Container::store_triangle(const ushort order) {
//  Reference_Geometry_Container::order_to_reference_triangle_.emplace(
//      order, std::make_shared<Reference_Triangle>(order));
//}
//
// void Reference_Geometry_Container::store_quadrilateral(const ushort order)
// {
//  Reference_Geometry_Container::order_to_reference_quadrilateral_.emplace(
//      order, std::make_shared<Reference_Quadrilateral>(order));
//}