#include <time.h>
#include "Timer.h"

Timer::Timer() : start_sec_(0), start_ns_(0), total_(0.0) {}

void Timer::GetWallTime(int& sec, int& ns) {
  struct timespec wall_time;
  clock_gettime(CLOCK_REALTIME, &wall_time);
  sec = wall_time.tv_sec;
  ns = wall_time.tv_nsec;
};
