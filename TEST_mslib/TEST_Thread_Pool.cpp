#include "pch.h"

#include "Thread_Pool.h"

TEST( Thread_Pool, test1 )
{
  using namespace mslib;

  constexpr uint64 num_worker_thread = 3;
  Thread_Pool      thread_pool( num_worker_thread );

  constexpr uint64 num_job = 5;
  auto             job     = []( const uint64 i, uint64* value ) { std::this_thread::sleep_for(std::chrono::milliseconds(300 * i)); *value = i; };
  // auto             job     = []( const uint64 i, uint64& value ) { std::this_thread::sleep_for(std::chrono::milliseconds(300 * i)); printf( " %llu th job done \n", i ); };

  std::vector<uint64> result( num_job );
  for ( uint64 i = 0; i < num_job; ++i )
    thread_pool.addJob( job, i, result.data() + i );

  thread_pool.wait();

  std::vector<uint64> ref = { 0, 1, 2, 3, 4 };
  EXPECT_EQ( result, ref );
}