#include "../_src/msconfig/Data.h"
#include "gtest/gtest.h"

TEST(DataSet, extract)
{
  ms::config::Data minseok;

  minseok.set<std::string>("name", "KimMinseok");
  minseok.set<std::string>("company", "MidasIT");
  minseok.set<int>("age", 29);
  minseok.set<double>("height", 171.1);

  const auto     data = minseok.extract_data<std::string>("name");
  constexpr auto ref  = "KimMinseok";
  EXPECT_EQ(data, ref);
}

TEST(DataSET, duplicated_key_test)
{
  ms::config::Data minseok;

  minseok.set<std::string>("name", "KimMinseok");
  EXPECT_ANY_THROW(minseok.set<std::string>("name", "John"));
}