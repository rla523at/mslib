#include "msconfig/Config_Data.h"

namespace ms::config
{

  void Data::set_string( std::string_view key, std::string_view str )
  {
    // key shold be upper case!
    const auto KEY = ms::string::upper_case( key );

    if ( this->_key_to_str.contains( KEY.data() ) )
    {
      EXCEPTION( std::string( key ) + "is duplicated" );
    }

    this->_key_to_str.insert( { KEY, std::string( str ) } );
  }

  void Data::set_number( std::string_view key, std::string_view number )
  {
    // key shold be upper case!
    const auto KEY = ms::string::upper_case( key );

    if ( this->_key_to_number.contains( KEY.data() ) )
    {
      EXCEPTION( std::string( key ) + "is duplicated" );
    }

    this->_key_to_number.insert( { KEY, std::string( number ) } );
  }

  void Data::set_bool( std::string_view key, const bool boolean )
  {
    // key shold be upper case!
    const auto KEY = ms::string::upper_case( key );

    if ( this->_key_to_bool.contains( KEY.data() ) )
    {
      EXCEPTION( std::string( key ) + "is duplicated" );
    }

    this->_key_to_bool.insert( { KEY, boolean } );
  }

  void Data::set_sub_data( std::string_view key, Data&& sub_data )
  {
    // key shold be upper case!
    const auto KEY = ms::string::upper_case( key );

    if ( this->_key_to_data.contains( KEY.data() ) )
    {
      EXCEPTION( std::string( key ) + "is duplicated" );
    }

    this->_key_to_data.insert( { KEY, std::move( sub_data ) } );
  }

} // namespace ms::config
