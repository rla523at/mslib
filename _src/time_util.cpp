#include "stdafx.h"

#include "time_util.h"

#include "exception.h"

namespace mslib::time_util
{

//
// Local Functions and Variables
//
namespace
{

Time_Point get_time_point(void)
{
  return std::chrono::steady_clock::now();
}
uint64 calculate_duration_ns(const Time_Point& tp1, const Time_Point& tp2)
{
  NanoSecond duration = tp2 - tp1;
  return static_cast<uint64>(duration.count());
}
double cal_ellapsed_time_ms(const Time_Point& tp1, const Time_Point& tp2)
{
  NanoSecond duration    = tp2 - tp1;
  MillSecond duration_ms = std::chrono::duration_cast<MillSecond>(duration);
  return duration_ms.count();
}

} // namespace

//
// TimeRecorder
//

void Time_Recorder::set_start_time_point(void)
{
  _start_time_point = get_time_point();
}
void Time_Recorder::record(void)
{
  REQUIRE(_start_time_point != Invalid_Time_Point, "start time point 가 설정되지 않았습니다.");

  const Time_Point end_time_point = get_time_point();

  const double ellapsed_time = cal_ellapsed_time_ms(_start_time_point, end_time_point);

  _count      += 1;
  _total_time += ellapsed_time;
  _min_time    = std::min(_min_time, ellapsed_time);
  _max_time    = std::max(_max_time, ellapsed_time);
  _avg_time    = _total_time / _count;
}

} // namespace mslib::time
