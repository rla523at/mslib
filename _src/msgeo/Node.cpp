#include "msgeo/Node.h"

#include "msexception/Exception.h"
#include "msmath/BLAS.h"
#include <algorithm>
#include <iomanip>
#include <sstream>

namespace ms::geo
{

  double Node_View::operator[]( const int index ) const
  {
    return this->at( index );
  }

  bool Node_View::operator==( const Node_View& other ) const
  {
    if ( this->_dimension != other._dimension ) return false;

    constexpr auto epsilon = 1.0e-10;

    for ( int i = 0; i < this->_dimension; ++i )
    {
      const auto abs_diff = std::abs( ( *this )[i] - other[i] );
      if ( epsilon < abs_diff ) return false;
    }

    return true;
  }

  ms::math::Vector<0> Node_View::operator-( const Node_View& other ) const
  {
    ms::math::Vector<0> result( this->_dimension );
    for ( int i = 0; i < this->_dimension; ++i )
    {
      result[i] = ( *this )[i] - other[i];
    }

    return result;
  }

  double Node_View::at( const int index ) const
  {
    REQUIRE( index < this->_dimension, "index can not exceed dimension" );
    return this->_coordinates_view[index * this->_inc];
  }

  int Node_View::dimension( void ) const
  {
    return this->_dimension;
  }

  double Node_View::distance( const Node_View other ) const
  {
    double result = 0.0;

    for ( int i = 0; i < this->_dimension; ++i )
    {
      const auto diff  = ( *this )[i] - other[i];
      result          += diff * diff;
    }

    return std::sqrt( result );
  }

  void Node_View::other_to_this_vector( const Node_View& other, ms::math::Vector_Wrap vector_wrap ) const
  {
    REQUIRE( vector_wrap.dimension() == this->_dimension, "vector wrap dimension should be equal to dimension" );

    for ( int i = 0; i < this->_dimension; ++i )
    {
      vector_wrap[i] = ( *this )[i] - other[i];
    }
  }

  int Node_View::size( void ) const
  {
    return this->_dimension;
  }

  std::string Node_View::to_string( void ) const
  {
    std::ostringstream oss;
    oss << std::setprecision( 16 ) << std::showpoint << std::left;
    for ( int i = 0; i < this->_dimension; ++i )
    {
      oss << std::setw( 25 ) << this->at( i );
    }
    return oss.str();
  }

  ms::math::Vector_View Node_View::to_vector_view( void ) const
  {
    const auto result = ms::math::Vector_View( this->_coordinates_view, this->_inc );
    return result;
  }

  const double* Node_View::data( void ) const
  {
    return this->_coordinates_view.data();
  }

} // namespace ms::geo

/*










*/

namespace ms::geo
{

  double& Node_Wrap::operator[]( const int index )
  {
    return this->at( index );
  }

  double& Node_Wrap::at( const int index )
  {
    REQUIRE( index < this->_dimension, "index can not exceed dimension" );
    return this->_coordinates_wrap[index * this->_inc];
  }

  ms::math::Vector_Wrap Node_Wrap::to_vector_wrap( void )
  {
    auto result = ms::math::Vector_Wrap( this->_coordinates_wrap, this->_inc );
    return result;
  }

} // namespace ms::geo

/*










*/

namespace ms::geo
{

  Node::Node( const int dim )
  {
    REQUIRE( 0 < dim, "dimension should be natural number" );
    REQUIRE( dim <= 3, "dimension can't not excced 3" );

    this->_dimension = dim;
    this->_inc       = 1;

    this->reallocate_coordinates();
  }

  Node::Node( std::initializer_list<double> coordinates )
  {
    this->_dimension = static_cast<int>( coordinates.size() );
    this->_inc       = 1;
    REQUIRE( 0 < this->_dimension, "dimension should be natural number" );
    REQUIRE( this->_dimension <= 3, "dimension can't not excced 3" );

    ms::math::blas::copy( this->_coordinates.data(), std::data( coordinates ), this->_dimension );
    this->reallocate_coordinates();
  }

  Node::Node( const Node& other )
    : _coordinates( other._coordinates )
  {
    this->_dimension = other._dimension;
    this->_inc       = other._inc;
    this->reallocate_coordinates();
  }

  Node::Node( Node&& other ) noexcept
    : _coordinates( other._coordinates )
  {
    this->_dimension = other._dimension;
    this->_inc       = other._inc;
    this->reallocate_coordinates();

    other._coordinates_view = std::span<const double>();
    other._coordinates_wrap = std::span<double>();
  }

  void Node::reallocate_coordinates( void )
  {
    this->_coordinates_view = this->_coordinates;
    this->_coordinates_wrap = this->_coordinates;
  }

} // namespace ms::geo

namespace ms::geo
{

  bool Node_Compare::operator()( Node_View lhs, Node_View rhs ) const
  {
    const auto dim = lhs.dimension();
    REQUIRE( dim == rhs.dimension(), "dimension should be matched" );

    constexpr auto epsilon = 1.0e-10;

    for ( int i = 0; i < dim; ++i )
    {
      if ( std::abs( rhs[i] - lhs[i] ) <= epsilon )
      {
        continue;
      }

      return lhs[i] < rhs[i];
    }

    return false;
  }

} // namespace ms::geo

/*```










*/

namespace ms::geo
{

  Node_View Nodes_View::operator[]( const int index ) const
  {
    return this->at( index );
  }

  bool Nodes_View::operator==( const Nodes_View other ) const
  {
    if ( this->_num_nodes != other._num_nodes ) return false;
    if ( this->_dimension != other._dimension ) return false;

    for ( int i = 0; i < this->_num_nodes; ++i )
    {
      if ( this->at( i ) != other.at( i ) ) return false;
    }

    return true;
  }

  Node_View Nodes_View::at( const int index ) const
  {
    // convention : index start with 0
    REQUIRE( index < this->_num_nodes, "index can not exceed number of nodes" );

    const auto start_index     = index * this->_dimension;
    const auto num             = this->_dimension;
    const auto sub_coordinates = this->_coordinates_view.subspan( start_index, num );
    const auto result          = Node_View( sub_coordinates );
    return result;
  }

  ms::math::Vector_View Nodes_View::axis_coordinates_vector_view( const int axis_index ) const
  {
    REQUIRE( axis_index < this->_dimension, "axis index can not exceed dimension" );

    const auto start_index          = axis_index;
    const auto sub_coordinates_view = this->_coordinates_view.subspan( start_index );

    const auto result = ms::math::Vector_View( sub_coordinates_view, this->_num_nodes );
    return result;
  }

  void Nodes_View::copy_coordinates( double* dest ) const
  {
    const auto data_ptr = this->_coordinates_view.data();
    std::copy( data_ptr, data_ptr + this->num_coordinates(), dest );
  }

  int Nodes_View::dimension( void ) const
  {
    return this->_dimension;
  }

  int Nodes_View::num_coordinates( void ) const
  {
    return this->_dimension * this->_num_nodes;
  }

  int Nodes_View::num_nodes( void ) const
  {
    return this->_num_nodes;
  }

  std::vector<Node_View> Nodes_View::nodes_at_indexes( const std::vector<int>& numbers ) const
  {
    const auto num_nodes = numbers.size();
    REQUIRE( num_nodes < this->_num_nodes, "number of numbers can not excced number of nodes" );

    std::vector<Node_View> result;
    result.reserve( num_nodes );

    for ( const auto index : numbers )
    {
      result.push_back( ( *this )[index] );
    }

    return result;
  }

  void Nodes_View::nodes_at_indexes( std::vector<Node_View>& result, const std::vector<int>& numbers ) const
  {
    const auto num_nodes = numbers.size();
    REQUIRE( num_nodes < this->_num_nodes, "number of numbers can not excced number of nodes" );

    result.clear();
    result.reserve( num_nodes );

    for ( const auto index : numbers )
    {
      result.push_back( ( *this )[index] );
    }
  }

}; // namespace ms::geo

/*










*/

namespace ms::geo
{

  Node_Wrap Nodes_Wrapper::operator[]( const int index )
  {
    return this->at( index );
  }

  Node_Wrap Nodes_Wrapper::at( const int index )
  {
    // convention : index start with 0
    REQUIRE( index < this->_num_nodes, "index can not exceed number of nodes" );

    const auto start_index     = index * this->_dimension;
    const auto num             = this->_dimension;
    const auto sub_coordinates = this->_coordinates_wrap.subspan( start_index, num );
    const auto result          = Node_Wrap( sub_coordinates );
    return result;
  }

} // namespace ms::geo

/*










*/

namespace ms::geo
{

  Nodes::Nodes( const int num_nodes, const int dimension )
  {
    this->_num_nodes = num_nodes;
    this->_dimension = dimension;
    REQUIRE( 0 < this->_num_nodes, "number of nodes should be natural number" );
    REQUIRE( 0 < this->_dimension, "dimension should be natural number" );

    this->_coordinates.resize( this->_num_nodes * this->_dimension );
    this->_coordinates_view = this->_coordinates;
    this->_coordinates_wrap = this->_coordinates;
  }

  Nodes::Nodes( const int num_nodes, const int dimension, std::vector<double>&& coordinates )
  {
    this->_num_nodes = num_nodes;
    this->_dimension = dimension;
    REQUIRE( 0 < this->_num_nodes, "number of nodes should be natural number" );
    REQUIRE( 0 < this->_dimension, "dimension should be natural number" );

    this->_coordinates = std::move( coordinates );
    REQUIRE( this->_coordinates.size() == this->_num_nodes * this->_dimension, "Given coordinates size should be matched with given size" );

    this->_coordinates_view = this->_coordinates;
    this->_coordinates_wrap = this->_coordinates;
  }

  Nodes::Nodes( const Nodes& other )
    : _coordinates( other._coordinates )
  {
    this->_num_nodes = other._num_nodes;
    this->_dimension = other._dimension;

    this->_coordinates_view = this->_coordinates;
    this->_coordinates_wrap = this->_coordinates;
  }

  Nodes::Nodes( Nodes&& other ) noexcept
    : _coordinates( std::move( other._coordinates ) )
  {
    this->_num_nodes = other._num_nodes;
    this->_dimension = other._dimension;

    this->_coordinates_view = this->_coordinates;
    this->_coordinates_wrap = this->_coordinates;

    other._coordinates_view = std::span<const double>();
    other._coordinates_wrap = std::span<double>();
  }

  void Nodes::operator=( const Nodes& other )
  {
    this->_coordinates = other._coordinates;
    this->_num_nodes   = other._num_nodes;
    this->_dimension   = other._dimension;

    this->_coordinates_view = this->_coordinates;
    this->_coordinates_wrap = this->_coordinates;
  }

  void Nodes::operator=( Nodes&& other )
  {
    this->_num_nodes = other._num_nodes;
    this->_dimension = other._dimension;

    this->_coordinates      = std::move( other._coordinates );
    this->_coordinates_view = this->_coordinates;
    this->_coordinates_wrap = this->_coordinates;

    other._coordinates_view = std::span<const double>();
    other._coordinates_wrap = std::span<double>();
  }

  void Nodes::add_node( const Node_View node )
  {
    const auto dim = node.dimension();

    if ( this->is_null_node() )
    {
      this->_coordinates.resize( dim );
      this->_dimension = dim;

      for ( int i = 0; i < dim; ++i )
      {
        this->_coordinates[i] = node[i];
      }
    }
    else
    {
      REQUIRE( this->_dimension == dim, "dimension should be matched for adding" );

      const auto c1 = this->_dimension * this->_num_nodes;
      this->_coordinates.resize( c1 + dim );

      for ( int i = 0; i < dim; ++i )
      {
        this->_coordinates[c1 + i] = node[i];
      }
    }

    this->_coordinates_view = this->_coordinates;
    this->_coordinates_wrap = this->_coordinates;
    this->_num_nodes++;
  }

  void Nodes::add_nodes( Nodes&& other )
  {
    if ( this->is_null_node() )
    {
      this->_coordinates = std::move( other._coordinates );
      this->_dimension   = other._dimension;
    }
    else
    {
      REQUIRE( this->_dimension == other._dimension, "dimension should be matched for adding" );
      this->_coordinates.insert( this->_coordinates.end(), other._coordinates.begin(), other._coordinates.end() );
    }

    this->_coordinates_view  = this->_coordinates;
    this->_coordinates_wrap  = this->_coordinates;
    this->_num_nodes        += other._num_nodes;

    other._coordinates.clear();
    other._coordinates_view = std::span<const double>();
    other._coordinates_wrap = std::span<double>();
  }

  bool Nodes::is_null_node( void ) const
  {
    return this->_coordinates.empty();
  }

} // namespace ms::geo

namespace ms::geo
{

  void A_to_B_vector( ms::math::Vector_Wrap vector_wrap, const Node_View A, const Node_View B )
  {
    const auto dimension = A.dimension();
    REQUIRE( dimension == B.dimension(), "dimension should be matched" );
    REQUIRE( vector_wrap.dimension() == dimension, "vector wrap dimension should be equal to dimension" );

    for ( int i = 0; i < dimension; ++i )
    {
      vector_wrap[i] = B[i] - A[i];
    }
  }

} // namespace ms::geo
