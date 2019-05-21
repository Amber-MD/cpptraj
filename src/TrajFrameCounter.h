#ifndef INC_TRAJFRAMECOUNTER_H
#define INC_TRAJFRAMECOUNTER_H
// Forward declarations
class ArgList;
/// Used to keep track of frame # during traj read.
class TrajFrameCounter {
  public:
    TrajFrameCounter();
    /// \return Current frame.
    int Current() const { return current_; }
    /// \return Total number of frames
    int TotalFrames() const { return total_frames_; }
    /// \return Total number of frames that will be read based on start/stop/offset.
    int TotalReadFrames() const { return total_read_frames_; }
    /// \return Start frame number.
    int Start() const { return start_; }
    /// \return End frame number.
    int Stop() const { return stop_; }
    /// \return offset
    int Offset() const { return offset_; }
    /// \return Number of frames processed since Begin called.
    int NumFramesProcessed() const { return numFramesProcessed_; }
    /// Prepare counter for use.
    void Begin() { numFramesProcessed_ = 0; current_ = start_; }
    /// Set total_frames_ and get start/stop/offset from ArgList, check for validity.
    int CheckFrameArgs(int, ArgList&);
    /// Print start/stop/offset info to screen
    void PrintFrameInfo() const;
    /// Print brief info to screen
    void PrintInfoLine(const char*) const;
    /// Check if processing is complete.
    inline int CheckFinished() { return !(current_ < stop_ || stop_ == -1); }
    /// Update current according to offset
    inline void UpdateCounters() { ++numFramesProcessed_; current_ += offset_; }
    /// \return previous frame number (last # before UpdateCounters was called).
    inline int PreviousFrameNumber() const { return current_ - offset_; }
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
