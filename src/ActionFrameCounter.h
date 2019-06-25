#ifndef INC_ACTIONFRAMECOUNTER_H
#define INC_ACTIONFRAMECOUNTER_H
// Forward declarations
class ArgList;
/// Internal frame counter, for processing a subset of frames.
class ActionFrameCounter {
  public:
    ActionFrameCounter();
    static const char* HelpText;
    int InitFrameCounter(ArgList&);
    /// \return true if frame should not be processed.
    int CheckFrameCounter(int frameNum) const {
/*      if (start_ == -1) return true;
      // Need this in case frames were skipped
      while (start_ < frameNum)
        start_ += offset_;
      if (frameNum != start_) return true;
      start_ += offset_;
      // Since frameNum will never be -1, this disables all future checks.
      if (start_ > stop_ && stop_ != -1) start_ = -1;*/
      if (stop_ != -1 && frameNum > stop_) return 1;
      if (frameNum < start_) return 2;
      if (offset_ == 1) return 0;
      if ( ((frameNum - start_) % offset_) != 0 ) return 3;
      return 0;
    }
    void FrameCounterInfo() const;
    void FrameCounterBrief() const;
    bool DefaultSettings() const { return (start_ == 0 && stop_ == -1 && offset_ == 1); }
  private:
    int start_;
    int stop_;
    int offset_;
};
#endif
