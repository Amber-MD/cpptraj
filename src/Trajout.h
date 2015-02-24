#ifndef INC_TRAJOUT_H
#define INC_TRAJOUT_H
#include "TrajectoryFile.h"
#include "Range.h"
#include "ActionFrameCounter.h"
/// Output trajectory class.
// FIXME: InitTrajWrite should also take # frames to write?
class Trajout : public TrajectoryFile {
  public:
    Trajout();
    virtual ~Trajout() {}
    /// Close output trajectory.
    virtual void EndTraj() = 0;
    /// Write a single frame.
    virtual int WriteSingle(int, Frame const&) = 0;
    /// Write an array of frames.
    virtual int WriteEnsemble(int, FramePtrArray const&) = 0;
    /// Print information on trajectory to be written.
    virtual void PrintInfo(int) const = 0;
    /// Prepare trajectory for writing for the given topology.
    virtual int InitTrajWrite(std::string const&, ArgList const&, Topology*,
                              TrajectoryFile::TrajFormatType) = 0;
    /// Perform topology-related setup if given topology Pindex matches.
    virtual int SetupTrajWrite(Topology*) = 0;

    int NumFramesProcessed()          const { return numFramesProcessed_; }
    bool TrajIsOpen()                 const { return trajIsOpen_;         }
    bool TrajoutAppend()              const { return append_;             }
    std::string const& TrajoutTitle() const { return title_;              }
    void SetTrajIsOpen(bool o)              { trajIsOpen_ = o;            }
    inline int CheckFrameRange(int);
    /// Write single frame, performing set up if needed.
    inline int WriteFrame(int, Topology*, Frame const&);
  protected:
    /// Grab keywords common to all trajouts, set/determine format if necessary
    int CommonTrajoutSetup(std::string const&, ArgList&, Topology*, 
                           TrajectoryFile::TrajFormatType&);
    /// Set up anything topology-related.
    int FirstFrameSetup(std::string const&, TrajectoryIO*, Topology*);
    /// Check format for append
    int CheckAppendFormat(std::string const&, TrajFormatType&);
    /// For ensemble trajouts, get range of members to write.
    Range MembersToWrite(std::string const&,int) const;
    /// Print information for TrajectoryIO
    void CommonInfo(TrajectoryIO*) const;
  private:
    std::string title_;                ///< Output trajectory title.
    int numFramesProcessed_;
    bool trajIsOpen_;                  ///< If true trajectory has been opened.
    bool nobox_;                       ///< If true do not put box information in output traj
    bool append_;                      ///< If true, append to this file.
    bool hasRange_;                    ///< If true a frame range is defined.
    Range FrameRange_;                 ///< List of frame numbers to write.
    Range::const_iterator rangeframe_; ///< If frame range defined, this is next frame in range.
    ActionFrameCounter frameCount_;    ///< Hold start/stop/offset values
};
// ---------- INLINE FUNCTION DEFINITIONS --------------------------------------
/** \return 1 if set should not be written. */
int Trajout::CheckFrameRange(int set) {
  if (hasRange_) {
    // If no more frames in the framerange, skip.
    if (rangeframe_ == FrameRange_.end()) return 1;
    // If this frame is not the next in the range, skip.
    if ( *rangeframe_ != set ) return 1;
    // This frame is next in the range. Advance FrameRange iterator.
    ++rangeframe_;
  } else {
    if (frameCount_.CheckFrameCounter( set )) return 1;
  }
  // Frame will be processed. Increment frame count.
  ++numFramesProcessed_;
  return 0;
}

int Trajout::WriteFrame(int set, Topology* tparmIn, Frame const& f) {
  int err = SetupTrajWrite(tparmIn);
  return err + WriteSingle(set, f);
}
#endif
