#include "msfilesystem/filesystem.h"
#include "gtest/gtest.h"

#ifdef _DEBUG

TEST(msfilesystem, copy_file_exception1)
{
  constexpr auto src = "Test/Copy/From/folder";
	constexpr auto dst = "Test/Copy/To/";
	EXPECT_ANY_THROW(ms::filesystem::copy_file(src, dst));	
}
#endif

TEST(msfilesystem, copy_file1)
{
  constexpr auto src = "Test/Copy/From/file1.txt";
  constexpr auto dst = "Test/Copy/To/file1.txt";
  ms::filesystem::copy_file(src, dst);

  const bool is_exist = ms::filesystem::is_exist_file("Test/Copy/To/file1.txt");
  EXPECT_TRUE(is_exist);

  if (is_exist)
	{
		ms::filesystem::remove_file("Test/Copy/To/file1.txt");
	}
}

TEST(msfilesystem, copy_file2)
{
  constexpr auto src = "Test/Copy/From/AlreadyExist.txt";
  constexpr auto dst = "Test/Copy/To/AlreadyExist.txt";

  EXPECT_NO_THROW(ms::filesystem::copy_file(src, dst, std::filesystem::copy_options::overwrite_existing));
}

TEST(msfilesystem, extract_file_name)
{
  constexpr auto src = "Test/Copy/From/file1.txt";

  constexpr auto ref = "file1.txt";
  EXPECT_EQ(ref, ms::filesystem::extract_file_name(src));
}

TEST(msfilesystem, has_this_extension)
{
  EXPECT_TRUE(ms::filesystem::has_this_extension("A.txt", ".txt"));
}

TEST(msfilesystem, is_folder_path1)
{
  EXPECT_TRUE(ms::filesystem::is_folder_path("A/"));
}
TEST(msfilesystem, is_folder_path2)
{
  EXPECT_FALSE(ms::filesystem::is_folder_path("Test"));
}

TEST(msfilesystem, folder_name1)
{
  const auto     folder_path = ms::filesystem::extract_folder_path("Test/C/a.txt");
  constexpr auto ref         = "Test/C/";
  EXPECT_EQ(folder_path, ref);
}

TEST(msfilesystem, is_exist_file1)
{
  EXPECT_TRUE(ms::filesystem::is_exist_file("Test/C/a.txt"));
}
TEST(msfilesystem, is_exist_file2)
{
  EXPECT_FALSE(ms::filesystem::is_exist_file("Test/C/c.txt"));
}
TEST(msfilesystem, is_exist_file3)
{
  // Folder paths are not case-sensitive
  EXPECT_TRUE(ms::filesystem::is_exist_file("Test/C/A.txt"));
}

TEST(msfilesystem, is_exist_folder1)
{
  EXPECT_TRUE(ms::filesystem::is_exist_folder("Test/A/"));
}
TEST(msfilesystem, is_exist_folder2)
{
  // Folder paths are not case-sensitive
  EXPECT_TRUE(ms::filesystem::is_exist_folder("Test/a/"));
}

TEST(msfilesystem, is_exist_folder3)
{
  EXPECT_FALSE(ms::filesystem::is_exist_folder("Test/Z/"));
}
TEST(msfilesystem, remove_empty_folder1)
{
  constexpr auto folder_path = "Test/C/";
  EXPECT_ANY_THROW(ms::filesystem::remove_empty_folder(folder_path));
}
TEST(msfilesystem, remove_folder1)
{
  std::string folder_path = "Test/D/";
  ms::filesystem::make_folder(folder_path);
  ms::filesystem::make_folder(folder_path + "A/");
  ms::filesystem::make_folder(folder_path + "B/");

  ms::filesystem::remove_folder(folder_path);
  EXPECT_FALSE(ms::filesystem::is_exist_folder(folder_path));
}