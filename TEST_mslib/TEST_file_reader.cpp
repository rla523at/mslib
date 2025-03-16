#include "pch.h"

#include "file_reader.h"

TEST( file_reader, read_file_parallel1 )
{
  constexpr const char* file_path = "rsc/file_reader_test.txt";
  {
    std::ofstream file( file_path );
    file << "a b c d e f g h i";
  }

  for (size_t num_threads = 1; num_threads < 10; ++num_threads )
  {
    const std::string result = mslib::file_reader::read_file_parallel( file_path, num_threads );
    EXPECT_EQ( "a b c d e f g h i", result );
  }
}


//TEST( file_reader, find_parallel_reading_start_positions_by_new_line1 )
//{
//  constexpr const char* file_path = "rsc/test.txt";
//  {
//    std::ofstream file( file_path );
//    file << "a b c";
//    file << "d e f";
//    file << "g h i";
//  }
//
//  std::vector<std::streampos> reference;
//  reference.push_back( 0 );
//
//  const std::vector<std::streampos> start_positions = mslib::file_reader::find_parallel_reading_start_positions_by_new_line( file_path, 1 );
//  EXPECT_EQ( reference, start_positions );
//}

//TEST( file_reader, find_parallel_reading_start_positions_by_new_line2 )
//{
//  constexpr const char* file_path = "rsc/test.txt";
//  {
//    std::ofstream file( file_path );
//    file << "a b c\n";
//    file << "d e f\n";
//    file << "g h i\n";
//  }
//
//  std::vector<std::streampos> reference;
//  reference.push_back( 0 );
//  {
//    std::ifstream file( file_path, std::ios::binary );
//    std::string   temp;
//    std::getline( file, temp );
//    std::getline( file, temp );
//
//    reference.push_back( file.tellg() );
//  }
//
//  const std::vector<std::streampos> start_positions = mslib::file_reader::find_parallel_reading_start_positions_by_new_line( file_path, 2 );
//  EXPECT_EQ( reference, start_positions );
//  {
//    std::ifstream file( file_path, std::ios::binary );
//    file.seekg( start_positions[0] );
//    EXPECT_EQ( file.get(), 'a' );
//
//    file.seekg( start_positions[1] );
//    EXPECT_EQ( file.get(), 'g' );
//  }
//}
//
//TEST( file_reader, find_parallel_reading_start_positions_by_new_line3 )
//{
//  constexpr const char* file_path = "rsc/test.txt";
//  {
//    std::ofstream file( file_path );
//    file << "a b c\n";
//    file << "d e f\n";
//    file << "g h i g h i g h i\n";
//  }
//
//  // ghi\n 이 있어서 EOF 라고 인식을 못한다. 그 다음줄도 있는걸로 보는듯하다..
//
//  std::vector<std::streampos> reference;
//  reference.push_back( 0 );
//
//  const std::vector<std::streampos> start_positions = mslib::file_reader::find_parallel_reading_start_positions_by_new_line( file_path, 2 );
//  EXPECT_EQ( reference, start_positions );
//
//  {
//    std::ifstream file( file_path, std::ios::binary );
//    file.seekg( start_positions[0] );
//    EXPECT_EQ( file.get(), 'a' );
//
//    //return eof 네 뭐지..
//    file.seekg( start_positions[1] );
//    EXPECT_EQ( file.get(), ' ' );
//  }
//}
//
//TEST( file_reader, find_parallel_reading_start_positions_by_new_line4 )
//{
//  constexpr const char* file_path = "rsc/test.txt";
//  {
//    std::ofstream file( file_path );
//    file << "a b c";
//    file << "d e f";
//    file << "g h i";
//  }
//
//  std::vector<std::streampos> reference;
//  reference.push_back( 0 );
//
//  const std::vector<std::streampos> start_positions = mslib::file_reader::find_parallel_reading_start_positions_by_new_line( file_path, 2 );
//  EXPECT_EQ( reference, start_positions );
//}
//
//TEST( file_reader, find_parallel_reading_start_positions_by_new_line5 )
//{
//  constexpr const char* file_path = "rsc/test.txt";
//  {
//    std::ofstream file( file_path );
//    file << "a b c\n";
//    file << "d e f\n";
//    file << "g h i g h i g h i";
//  }
//
//  std::vector<std::streampos> reference;
//  reference.push_back( 0 );
//
//  const std::vector<std::streampos> start_positions = mslib::file_reader::find_parallel_reading_start_positions_by_new_line( file_path, 2 );
//  EXPECT_EQ( reference, start_positions );
//}