#include "stdafx.h"

#include "file_reader.h"

#include "exception.h"
#include "file_system.h"
#include <thread>

namespace mslib::file_reader
{
  namespace
  {
    uint64 file_data_size( _In_ const std::string_view file_path );
    void   read_file( _In_ const std::string_view file_path, _In_ const std::streamoff read_offset, _In_ const uint64 read_size, _Inout_ utf8* data_ptr );

  } // namespace

  std::string read_file_parallel( _In_ const std::string_view file_path, _In_ const uint64 num_threads )
  {
    REQUIRE( num_threads != 0, "쓰레드의 개수가 0 일 수는 없습니다.", num_threads );
    REQUIRE( mslib::filesystem::is_exist_file( file_path ), "{} 는 존재하지 않은 파일입니다. 경로를 확인해주세요.", file_path );

    const uint64 data_size = file_data_size( file_path );

    std::string result;
    result.resize( data_size );

    std::vector<std::thread> workers;
    workers.reserve( num_threads );

    std::streamoff read_offset = 0;
    const uint64   read_size   = ( data_size + num_threads - 1 ) / num_threads;
    utf8*          data_ptr    = result.data();
    for ( uint64 i = 0; i < num_threads; ++i )
    {
      workers.push_back( std::thread( read_file, file_path, read_offset, read_size, data_ptr ) );
      read_offset += read_size;
      data_ptr    += read_size;
    }

    for ( std::thread& worker : workers )
      worker.join();

    return result;
  }

  std::vector<std::streampos> find_parallel_reading_start_positions_by_new_line( const std::string_view file_path, const uint64 num_threads )
  {
    REQUIRE( num_threads != 0, "쓰레드의 개수가 0 일 수는 없습니다.", num_threads );
    REQUIRE( mslib::filesystem::is_exist_file( file_path ), "{} 는 존재하지 않은 파일입니다. 경로를 확인해주세요.", file_path );

    std::ifstream file( file_path.data(), std::ios::binary );
    REQUIRE( file.is_open(), "{} 를 open 하는데 실패하였습니다.", file_path );

    file.seekg( 0, std::ios_base::end );
    const std::streampos end_position = file.tellg();

    std::vector<std::streampos> start_positions;

    const uint64 num_maximum_positions = num_threads;
    start_positions.reserve( num_maximum_positions );

    const std::streampos first_start_position = 0;
    start_positions.push_back( first_start_position );

    const std::streampos offset_unit = end_position / num_threads; // streamoff (long long) 으로 형변환되어 / 연산이 수행되고 다시 streamoff 에서 streampos로 형변환된다.

    std::string temp;

    for ( uint64 i = 1; i < num_maximum_positions; ++i )
    {
      const std::streamoff position_candidate = i * offset_unit;

      if ( position_candidate <= start_positions.back() )
        continue;

      file.seekg( position_candidate );
      std::getline( file, temp );

      // 공백문자를 제외하고 파일의 끝인지를 판단할 수 있어야 한다.
      // 동작이 명확하게 정의가 되어 있지 않다.
      // 병렬로 읽을 때 줄을 나눠서 읽어야 되는 이유가 뭐지?
      // 하고 싶은게 뭐지? --> Matrix 형태로 적힌 파일을 읽어서 각각의 행을 병렬로 읽어서 처리하고 싶다.

      if ( file.eof() )
      {
        // file.clear();
        // start_positions.push_back( file.tellg() );
        break;
      }

      start_positions.push_back( file.tellg() );
    }

    return start_positions;
  }

  std::streampos find_line_start_position( std::ifstream& file, const std::streampos start_position )
  {
    REQUIRE( file.is_open(), "파일이 열려있지 않습니다." );
    REQUIRE( start_position >= 0, "주어진 시작 위치({:d})는 0 이상이어야 합니다.", long long( start_position ) );
    REQUIRE( file.fail() == false, "파일의 상태가 fail 상태입니다." );

    if ( file.rdstate() == std::ios_base::eofbit )
      file.clear();

    const std::streampos initial_position = file.tellg();

    std::streampos line_start_position = start_position;
    while ( line_start_position > 0 )
    {
      const std::streampos next_position = line_start_position - std::streampos( 1 );

      file.seekg( next_position );

      if ( file.get() == '\n' )
        break;

      line_start_position = next_position;
    }

    file.seekg( initial_position );
    return line_start_position;
  }

  namespace
  {
    uint64 file_data_size( _In_ const std::string_view file_path )
    {
      std::ifstream file( file_path.data(), std::ios::binary );
      REQUIRE( file.is_open(), "{} 를 open 하는데 실패하였습니다.", file_path );

      file.seekg( 0, std::ios_base::end );
      const std::streampos end_position = file.tellg();
      const uint64         data_size    = static_cast<uint64>( end_position );

      return data_size;
    }

    void read_file( _In_ const std::string_view file_path, _In_ const std::streamoff read_offset, _In_ const uint64 read_size, _Inout_ utf8* data_ptr )
    {
      std::ifstream file( file_path.data(), std::ios::binary );
      REQUIRE( file.is_open(), "{} 를 open 하는데 실패하였습니다.", file_path );

      // input validation requirement 추가하기
      file.seekg( read_offset );
      file.read( data_ptr, read_size );
    }
  } // namespace

} // namespace mslib::file_reader