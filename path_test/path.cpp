#include "mspath/path.h"
#include "gtest/gtest.h"

TEST(mspath, is_folder_path1)
{
  EXPECT_TRUE(ms::path::is_folder_path("A/"));
}
TEST(mspath, is_folder_path2)
{
  EXPECT_FALSE(ms::path::is_folder_path("Test"));
}

TEST(mspath, folder_name1)
{
  const auto     folder_path = ms::path::extract_folder_path("Test/C/a.txt");
  constexpr auto ref         = "Test/C/";
  EXPECT_EQ(folder_path, ref);
}

TEST(mspath, is_exist_file1)
{
  EXPECT_TRUE(ms::path::is_exist_file("Test/C/a.txt"));
}
TEST(mspath, is_exist_file2)
{
  EXPECT_FALSE(ms::path::is_exist_file("Test/C/c.txt"));
}
TEST(mspath, is_exist_file3)
{
  // Folder paths are not case-sensitive
  EXPECT_TRUE(ms::path::is_exist_file("Test/C/A.txt"));
}

TEST(mspath, is_exist_folder1)
{
  EXPECT_TRUE(ms::path::is_exist_folder("Test/A/"));
}
TEST(mspath, is_exist_folder2)
{
  // Folder paths are not case-sensitive
  EXPECT_TRUE(ms::path::is_exist_folder("Test/a/"));
}

TEST(mspath, is_exist_folder3)
{
  EXPECT_FALSE(ms::path::is_exist_folder("Test/Z/"));
}
TEST(mspath, remove_empty_folder1)
{
  constexpr auto folder_path = "Test/C/";
  EXPECT_ANY_THROW(ms::path::remove_empty_folder(folder_path));
}
TEST(mspath, remove_folder1)
{
  std::string folder_path = "Test/D/";
  ms::path::make_folder(folder_path);
  ms::path::make_folder(folder_path + "A/");
  ms::path::make_folder(folder_path + "B/");

  ms::path::remove_folder(folder_path);
  EXPECT_FALSE(ms::path::is_exist_folder(folder_path));
}