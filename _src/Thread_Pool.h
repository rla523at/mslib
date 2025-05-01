#pragma once

#include <condition_variable>
#include <mutex>
#include <queue>
#include <semaphore>
#include <thread>

namespace mslib
{
  inline constexpr uint64 MAX_NUM_WORKER_THREAD = 16;

  class Thread_Pool_Worker;

  class Thread_Pool
  {
  public:
    Thread_Pool( _In_ const uint64 num_worker_thread );
    ~Thread_Pool( void );

    template <typename Func, typename... Args>
    void addJob( Func&& func, Args&&... args )
    {
      {
        std::unique_lock<std::mutex> lock( _job_queue_mutex );
        _job_queue.emplace( std::bind( std::forward<Func>( func ), std::forward<Args>( args )... ) );
      }

      _semaphore.release();
      ++_num_pending_job;
    }

    void wait( void );

  protected:
    std::vector<std::unique_ptr<Thread_Pool_Worker>>  _worker_thread_arr;
    std::counting_semaphore<MAX_NUM_WORKER_THREAD>    _semaphore;
    std::queue<std::function<void()>>                 _job_queue;
    std::mutex                                        _job_queue_mutex;
    std::condition_variable                           _empty_job_state;
    std::atomic<uint64>                               _num_pending_job = 0;
    bool                                              _quit            = false;
  };
} // namespace ms
