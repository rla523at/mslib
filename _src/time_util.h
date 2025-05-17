#pragma once

namespace mslib::time_util
{
using Time_Point = std::chrono::steady_clock::time_point;
using NanoSecond = std::chrono::nanoseconds;
using MillSecond = std::chrono::duration<double, std::milli>;

constexpr Time_Point Invalid_Time_Point = Time_Point();

struct Time_Recorder
{
  void set_start_time_point(void);
  void record(void);

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
  Time_Recorder_Helper_Scope(Time_Recorder& time_recorder)
    : _time_recorder(time_recorder)
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

#define RECORD_THIS_SCOPE(time_recorder) mslib::time_util::Time_Recorder_Helper_Scope time_recorder_helper_scope(time_recorder);

} // namespace mslib::time
