#include "pch.h"

#include "time_util.h"

using namespace mslib::time_util;

TEST(time_util, 1)
{
  constexpr int N = 1000000;

  std::vector<double> value(N);
  for (int i = 0; i < N; ++i)
  {
    value[i] = static_cast<double>(i);
  }


  Time_Recorder time_recorder1;
  {
    double sum = 0;

    RECORD_THIS_SCOPE(time_recorder1);
    for (int i = 0; i < N; ++i)
    {
      sum += value[i];
    }
  }

  Time_Recorder time_recorder2;
  {
    double sum = 0;

    for (int i = 0; i < N; ++i)
    {
      RECORD_THIS_SCOPE(time_recorder2);
      sum += value[i];
    }
  }

  EXPECT_TRUE(true);
}
