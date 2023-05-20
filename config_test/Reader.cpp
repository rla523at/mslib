#include "../_src/msconfig/Reader.h"
#include "../_src/msconfig/Data.h"
#include "gtest/gtest.h"

TEST(Reader, read_bool1)
{

  const auto config = ms::config::Reader::read("TEST/config.dat");

  const auto     result = config.get_data<bool>("is_config");
  constexpr auto ref    = true;
  EXPECT_EQ(result, ref);
}
TEST(Reader, read_bool2)
{
  const auto config = ms::config::Reader::read("TEST/config.dat");

  const auto     result = config.get_data<bool>("Log_on");
  constexpr auto ref    = false;
  EXPECT_EQ(result, ref);
}
TEST(Reader, read_Data_Set)
{
  const auto config = ms::config::Reader::read("TEST/config.dat");

  const auto problem_option   = config.get_data<ms::config::Data>("Problem_Option");
  const auto sine_wave_option = problem_option.get_data<ms::config::Data>("sine_wave_option");

  const auto     result = sine_wave_option.get_data<double>("x_wave_length");
  constexpr auto ref    = 1.0;
  EXPECT_EQ(result, ref);
}