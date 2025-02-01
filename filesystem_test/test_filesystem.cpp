#include "msfilesystem.h"
#include "gtest/gtest.h"

#include <fstream>

#ifdef _DEBUG

// TEST( msfilesystem, copy_file_exception1 )
//{
//   constexpr auto src = "Test/Copy/From/folder";
//   constexpr auto dst = "Test/Copy/To/";
//   EXPECT_ANY_THROW( ms::filesystem::copy_file( src, dst ) );
// }
#endif

TEST( msfilesystem, base_name_view1 )
{
  constexpr auto file_path = L"Test/move_file/A.txt";
  const auto     result    = ms::filesystem::base_name_view( file_path );
  const auto     ref       = L"A";

  EXPECT_EQ( result, ref );
}

TEST( msfilesystem, file_name_view1 )
{
  constexpr auto file_path = L"Test/move_file/A.txt";
  const auto     result    = ms::filesystem::file_name_view( file_path );
  const auto     ref       = L"A.txt";

  EXPECT_EQ( result, ref );
}

TEST( msfilesystem, folder_path_view1 )
{
  constexpr auto file_path = L"Test/move_file/A.txt";
  const auto     result    = ms::filesystem::folder_path_view( file_path );
  const auto     ref       = L"Test/move_file/";

  EXPECT_EQ( result, ref );
}

TEST( msfilesystem, move_file1 )
{
  constexpr auto file_path       = "Test/move_file/A.txt";
  constexpr auto new_folder_path = "Test/move_file/A/";
  ms::filesystem::move_file( file_path, new_folder_path );

  EXPECT_TRUE( ms::filesystem::is_exist_file( "Test/move_file/A/A.txt" ) );

  ms::filesystem::move_file( "Test/move_file/A/A.txt", "Test/move_file/" );
}

TEST( msfilesystem, copy_file1 )
{
  constexpr auto src = "Test/Copy/From/file1.txt";
  constexpr auto dst = "Test/Copy/To/file1.txt";
  ms::filesystem::copy_file( src, dst );

  const bool is_exist = ms::filesystem::is_exist_file( "Test/Copy/To/file1.txt" );
  EXPECT_TRUE( is_exist );

  if ( is_exist )
  {
    ms::filesystem::remove_file( "Test/Copy/To/file1.txt" );
  }
}

TEST( msfilesystem, copy_file2 )
{
  constexpr auto src = "Test/Copy/From/AlreadyExist.txt";
  constexpr auto dst = "Test/Copy/To/AlreadyExist.txt";

  EXPECT_NO_THROW( ms::filesystem::copy_file( src, dst, std::filesystem::copy_options::overwrite_existing ) );
}

TEST( msfilesystem, replace_file1 )
{
  constexpr auto src = "Test/Replace/A/file.txt";
  constexpr auto dst = "Test/Replace/B/B.txt";

  std::string temp;

  std::ifstream file( dst );
  std::getline( file, temp );
  EXPECT_EQ( temp, "B" );
  file.close();

  ms::filesystem::replace_file( src, dst );

  file.open( dst );
  std::getline( file, temp );
  EXPECT_EQ( temp, "A" );
  file.close();

  // roll back
  std::ofstream out_file( dst );
  out_file << "B";
}

TEST( msfilesystem, extract_file_name )
{
  constexpr auto src = "Test/Copy/From/file1.txt";

  constexpr auto ref = "file1.txt";
  EXPECT_EQ( ref, ms::filesystem::extract_file_name( src ) );
}
TEST( msfilesystem, extract_folder_path )
{
  constexpr auto src = "D:/Code/mslib/x64/Debug/filesystem_test.exe";

  constexpr auto ref = "D:/Code/mslib/x64/Debug/";
  EXPECT_EQ( ref, ms::filesystem::extract_folder_path( src ) );
}

TEST( msfilesystem, has_this_extension )
{
  EXPECT_TRUE( ms::filesystem::has_this_extension( "A.txt", ".txt" ) );
}

TEST( msfilesystem, is_folder_path1 )
{
  EXPECT_TRUE( ms::filesystem::is_folder_path( "A/" ) );
}
TEST( msfilesystem, is_folder_path2 )
{
  EXPECT_FALSE( ms::filesystem::is_folder_path( "Test" ) );
}

TEST( msfilesystem, folder_name1 )
{
  const auto     folder_path = ms::filesystem::extract_folder_path( "Test/C/a.txt" );
  constexpr auto ref         = "Test/C/";
  EXPECT_EQ( folder_path, ref );
}

TEST( msfilesystem, is_exist_file1 )
{
  EXPECT_TRUE( ms::filesystem::is_exist_file( "Test/C/a.txt" ) );
}
TEST( msfilesystem, is_exist_file2 )
{
  EXPECT_FALSE( ms::filesystem::is_exist_file( "Test/C/c.txt" ) );
}
TEST( msfilesystem, is_exist_file3 )
{
  // Folder paths are not case-sensitive
  EXPECT_TRUE( ms::filesystem::is_exist_file( "Test/C/A.txt" ) );
}

TEST( msfilesystem, is_exist_folder1 )
{
  EXPECT_TRUE( ms::filesystem::is_exist_folder( "Test/A/" ) );
}
TEST( msfilesystem, is_exist_folder2 )
{
  // Folder paths are not case-sensitive
  EXPECT_TRUE( ms::filesystem::is_exist_folder( "Test/a/" ) );
}

TEST( msfilesystem, is_exist_folder3 )
{
  EXPECT_FALSE( ms::filesystem::is_exist_folder( "Test/Z/" ) );
}
TEST( msfilesystem, remove_empty_folder1 )
{
  constexpr auto folder_path = "Test/C/";
  EXPECT_ANY_THROW( ms::filesystem::remove_empty_folder( folder_path ) );
}
TEST( msfilesystem, remove_folder1 )
{
  std::string reference_path = "Test/is_modefied_later/reference.txt";
  std::string target_path    = "Test/is_modefied_later/target.txt";

  EXPECT_TRUE( ms::filesystem::is_modified_later( reference_path, target_path ) );
}
TEST( msfilesystem, executable_path1 )
{
  std::string reference_path = "D:/Code/mslib/x64/Debug/filesystem_test.exe";
  std::string result         = ms::filesystem::excutable_file_path();

  EXPECT_EQ( reference_path, result );
}
TEST( msfilesystem, parent_path_view1 )
{
  constexpr auto src = "D:/Code/mslib/x64/Debug/filesystem_test.exe";

  constexpr auto ref = "D:/Code/mslib/";
  EXPECT_EQ( ref, ms::filesystem::parent_path_view( src, 2 ) );
}
TEST( msfilesystem, parent_path_view2 )
{
  constexpr auto src = "D:/Code/mslib/x64/Debug/filesystem_test.exe";

  constexpr auto ref = "";
  EXPECT_EQ( ref, ms::filesystem::parent_path_view( src, 5 ) );
}