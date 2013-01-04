#ifndef INC_ACTIONFRAMECOUNTER_H
#define INC_ACTIONFRAMECOUNTER_H
#include "ArgList.h"
/// Internal frame counter, for processing a subset of frames.
class ActionFrameCounter {
  public:
    ActionFrameCounter();
    static void Help();
    int InitFrameCounter(ArgList&);
    /// \return true if frame should not be processed.
    bool CheckFrameCounter(int frameNum) {
      if (frameNum != start_) return true;
      start_ += offset_;
      // Since frameNum will never be -1, this disables all future checks.
      if (start_ > stop_ && stop_ != -1) start_ = -1;
      return false;
    }
    void FrameCounterInfo();
  private:
    int start_;
    int stop_;
    int offset_;
};
#endif
