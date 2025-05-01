#include "stdafx.h"

#include "Thread_Pool.h"
#include "exception.h"
#include <windows.h> // SetThreadDescription

//
// class Private_Thread_Pool
//
namespace mslib
{
  class Private_Thread_Pool : public Thread_Pool
  {
  public:
    using Thread_Pool::_empty_job_state;
    using Thread_Pool::_job_queue;
    using Thread_Pool::_job_queue_mutex;
    using Thread_Pool::_num_pending_job;
    using Thread_Pool::_quit;
    using Thread_Pool::_semaphore;
  };
} // namespace mslib

//
// class Thread_Pool_Worker
//
namespace mslib
{
  class Thread_Pool_Worker
  {
  public:
    Thread_Pool_Worker( const std::string_view name, Thread_Pool* thread_pool )
      : _name( name ), _thread_pool( thread_pool )
    {
      _thread = std::thread( [this]() { this->initalize(); } );
    }

    void join( void )
    {
      _thread.join();
    }

  private:
    void initalize( void ) const
    {
      set_thread_name();

      Private_Thread_Pool& thread_pool = static_cast<Private_Thread_Pool&>( *_thread_pool );

      while ( true )
      {
        thread_pool._semaphore.acquire();

        if ( thread_pool._quit )
          break;

        std::function<void()> job = getJob();

        job();

        check_empty_job_state();
      }
    }

    std::function<void()> getJob( void ) const
    {
      Private_Thread_Pool& thread_pool = static_cast<Private_Thread_Pool&>( *_thread_pool );

      std::unique_lock<std::mutex> lock( thread_pool._job_queue_mutex );
      REQUIRE( thread_pool._job_queue.empty() == false, "semaphore acquire 에 성공했을 때, job queue 는 Empty 일 수 없습니다." );

      std::function<void()> job = std::move( thread_pool._job_queue.front() );
      thread_pool._job_queue.pop();

      return job;
    }

    void check_empty_job_state( void ) const
    {
      Private_Thread_Pool& thread_pool = static_cast<Private_Thread_Pool&>( *_thread_pool );

      const uint64 num_pending_job = --thread_pool._num_pending_job;
      if ( num_pending_job == 0 )
        thread_pool._empty_job_state.notify_all();
    }

    void set_thread_name( void ) const
    {
      SetThreadDescription( GetCurrentThread(), std::wstring( _name.begin(), _name.end() ).c_str() );
    }

  protected:
    Thread_Pool* _thread_pool;
    std::string  _name;
    std::thread  _thread;
  };

} // namespace mslib

//
// class Thread_Pool
//
namespace mslib
{
  Thread_Pool::Thread_Pool( _In_ const uint64 num_worker_thread )
    : _semaphore( 0 )
  {
    REQUIRE( 0 < num_worker_thread && num_worker_thread <= MAX_NUM_WORKER_THREAD, "num_worker_thread({:d}) should be in range [1, MAX_NUM_WORKER_THREAD]", num_worker_thread );

    _worker_thread_arr.reserve( num_worker_thread );

    for ( uint64 i = 0; i < num_worker_thread; ++i )
      _worker_thread_arr.emplace_back( std::make_unique<Thread_Pool_Worker>( std::format( "{:d} thread pool worker", i ), this ) );
  }

  Thread_Pool::~Thread_Pool( void )
  {
    _quit = true;

    _semaphore.release( _worker_thread_arr.size() );
    for ( auto& worker : _worker_thread_arr )
      worker->join();
  }

  void Thread_Pool::wait( void )
  {
    std::mutex                   mutex;
    std::unique_lock<std::mutex> lock( mutex );
    _empty_job_state.wait( lock, [this] { return ( _num_pending_job == 0 ); } );
  }

} // namespace mslib
