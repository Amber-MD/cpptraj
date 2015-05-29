#ifndef INC_TRAJFRAMECOUNTER_H
#define INC_TRAJFRAMECOUNTER_H
#include "ArgList.h"
/// Used to keep track of frame # during traj read.
class TrajFrameCounter {
  public:
    TrajFrameCounter();
    /// \return Current frame.
    int Current() const { return current_; }
    /// \return Total number of frames
    int TotalFrames() const { return total_frames_; }
    /// \return Start frame number.
    int Start() const { return start_; }
    /// Prepare counter for use.
    void Begin() { numFramesProcessed_ = 0; current_ = start_; }
    /// Set total number of frames.
    int SetTotalFrames(int);
    /// Get start/stop/offset from ArgList, check against total_frames_
    int CheckFrameArgs(ArgList&);
    /// Print start/stop/offset info to screen
    void PrintFrameInfo() const;
    /// Check if processing is complete.
    inline int CheckFinished() { return (current_ > stop_ && stop_ != -1); }
    /// Update current according to offset
    inline void UpdateCounters() { ++numFramesProcessed_; current_ += offset_; }
  private:
    int start_;  ///< Frame to begin processing.
    int stop_;   ///< Frame to end processing. -1 means until no more frames.
    int offset_; ///< Number of frames to skip between processed frames.
    int total_frames_; ///< Number of frames in the trajectory.
    int total_read_frames_; ///< # of frames that will be read based on start/stop/offset
    int current_; ///< The current frame number being read.
    int numFramesProcessed_; ///< Number of frames that have been read since beginning.
};
#endif
