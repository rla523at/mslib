#pragma once

namespace mslib::time_util
{
  using Clock_Type = std::chrono::steady_clock;
  using Time_Point = Clock_Type::time_point;
  using Nano_Second = std::chrono::nanoseconds;
  using Mill_Second = std::chrono::duration<double, std::milli>;

  constexpr Time_Point Invalid_Time_Point = Time_Point();

  inline Time_Point get_time_point( void )
  {
    return Clock_Type::now();
  }
  inline uint64 calculate_duration_ns( const Time_Point& tp1, const Time_Point& tp2 )
  {
    Nano_Second duration = tp2 - tp1;
    return static_cast<uint64>( duration.count() );
  }
  inline double cal_ellapsed_time_ms( const Time_Point& tp1, const Time_Point& tp2 )
  {
    Nano_Second duration    = tp2 - tp1;
    Mill_Second duration_ms = std::chrono::duration_cast<Mill_Second>( duration );
    return duration_ms.count();
  }

  struct Time_Recorder
  {
    inline void set_start_time_point(void)
    {
      _start_time_point = get_time_point();
    }
    inline void record(void)
    {
      const Time_Point end_time_point = get_time_point();

      const double ellapsed_time = cal_ellapsed_time_ms( _start_time_point, end_time_point );

      _count      += 1;
      _total_time += ellapsed_time;
      _min_time    = std::min( _min_time, ellapsed_time );
      _max_time    = std::max( _max_time, ellapsed_time );
      _avg_time    = _total_time / _count;  
    }

    uint64     _count            = 0;
    double     _total_time       = 0;
    double     _max_time         = 0;
    double     _min_time         = std::numeric_limits<double>::max();
    double     _avg_time         = 0.0;
    Time_Point _start_time_point = Invalid_Time_Point;
  };

  class Time_Recorder_Helper_Scope
  {
  public:
    Time_Recorder_Helper_Scope( Time_Recorder& time_recorder )
      : _time_recorder( time_recorder )
    {
      _time_recorder.set_start_time_point();
    }
    ~Time_Recorder_Helper_Scope()
    {
      _time_recorder.record();
    }

  private:
    Time_Recorder& _time_recorder;
  };

#define RECORD_THIS_SCOPE( time_recorder ) mslib::time_util::Time_Recorder_Helper_Scope time_recorder_helper_scope( time_recorder );

} // namespace mslib::time_util
