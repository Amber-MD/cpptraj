#ifndef INC_PROGRESSTIMER_H
#define INC_PROGRESSTIMER_H
#include "Timer.h"
/// Print estimate of time remaining based on current progress.
class ProgressTimer {
  public:
    ProgressTimer() : target_(0.0), interval_(0.0), total_it_(0) {}
    ProgressTimer(int t) : target_(5.0), interval_(5.0), total_it_(t) { progress_.Start(); }
    ProgressTimer(int t, double i) : target_(i), interval_(i), total_it_(t) { progress_.Start(); }
    void Remaining(int);
  private:
    Timer progress_;  ///< Keep track of elapsed time.
    double target_;   ///< Next target time to print progress in seconds.
    double interval_; ///< Interval between updates in seconds.
    int total_it_;    ///< Total number of iterations process should take.
};
#endif
