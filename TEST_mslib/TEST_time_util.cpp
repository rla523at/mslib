#include "pch.h"

#include "time_util.h"

using namespace mslib::time_util;

TEST( time_util, 1 )
{
  constexpr int N = 1000000;

  std::vector<double> value( N );
  for ( int i = 0; i < N; ++i )
  {
    value[i] = static_cast<double>( i );
  }

  Time_Recorder time_recorder1;
  {
    double sum = 0;

    RECORD_THIS_SCOPE( time_recorder1 );
    for ( int i = 0; i < N; ++i )
    {
      sum += value[i];
    }
  }

  Time_Recorder time_recorder2;
  {
    double sum = 0;

    for ( int i = 0; i < N; ++i )
    {
      RECORD_THIS_SCOPE( time_recorder2 );
      sum += value[i];
    }
  }

  EXPECT_TRUE( time_recorder1._count == 1 );
  EXPECT_TRUE( time_recorder2._count == N );
}

TEST( time_util, 2 )
{
  constexpr int N = 100000;

  std::vector<double> value( N );
  for ( int i = 0; i < N; ++i )
  {
    value[i] = static_cast<double>( i );
  }

  double              sum  = 0;
  double              sum2 = 0;
  std::vector<double> temp;
  std::vector<double> temp2;

  Time_Recorder time_recorder1;
  {
    {
      RECORD_THIS_SCOPE( time_recorder1 );
      for ( int i = 0; i < N; ++i )
      {
         //RECORD_THIS_SCOPE( time_recorder1 );
        temp.push_back( value[i] );
      }
    }

    sum += temp[1234];
  }

  Time_Recorder time_recorder2;
  {
    temp2.reserve( N );

    {
      RECORD_THIS_SCOPE( time_recorder2 );
      for ( int i = 0; i < N; ++i )
      {
        //RECORD_THIS_SCOPE( time_recorder2 );
        temp2.push_back( value[i] );
      }
    }

    sum2 += temp2[1234];
  }
  std::cout << time_recorder1._total_time << "\n";
  std::cout << time_recorder2._total_time << "\n";

  EXPECT_EQ( sum, sum2 );

  // EXPECT_TRUE( time_recorder1._count == 1 );
  // EXPECT_TRUE( time_recorder2._count == N );
}