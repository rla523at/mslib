#pragma once

namespace mslib::file_reader
{
  std::string read_file_parallel( _In_ const std::string_view file_path, _In_ const uint64 num_threads );

  /*
  num threads 만큼의 쓰레드가 병렬로 읽기 위한 행의 시작 위치를 찾아서 반환한다.    
    * 행의 수 보다 thread 가 많을 경우 반환되는 시작 위치의 수가 thread 수보다 적을 수 있다.
  */
  std::vector<std::streampos> find_parallel_reading_start_positions_by_new_line( const std::string_view file_path, const uint64 num_threads );

  /*
  start position 위치를 포함하는 행의 시작 위치를 찾아서 반환한다.
  */
  std::streampos find_line_start_position( std::ifstream& file, const std::streampos start_position );
} // namespace mslib::file_reader