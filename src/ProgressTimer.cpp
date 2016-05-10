#include "ProgressTimer.h"
#include "CpptrajStdio.h"

void ProgressTimer::Remaining(int iterations) {
  double elapsed = progress_.Elapsed();
  if (elapsed > target_) {
    target_ += interval_;
    int it_remaining = total_it_ - iterations;
    double it_per_s = ((double)iterations) / elapsed;
    double time_remaining = ((double)it_remaining) / it_per_s;
    mprintf("\t%i iterations in %g s, %g s remaining.\n", iterations, elapsed, time_remaining);
  }
}
