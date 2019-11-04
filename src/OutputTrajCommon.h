#ifndef INC_OUTPUTTRAJCOMMON_H
#define INC_OUTPUTTRAJCOMMON_H
#include "TrajectoryFile.h"
#include "Range.h"
#include "ActionFrameCounter.h"
#include "CoordinateInfo.h"
// Forward declarations
class Topology;
/// Output trajectory/ensemble common functionality.
class OutputTrajCommon {
  public:
    OutputTrajCommon();
    /// \return trajout file name
    FileName const& Filename()                    const { return trajName_; }
    // TODO: return const&, modify all TrajectoryIO routines?
    /// \return Pointer to associated topology
    Topology* Parm()                              const { return trajParm_; }
    /// \return true if trajectory should be appended to.
    bool Append()                                 const { return append_; }
    /// \return Title
    std::string const& Title()                    const { return title_; }
    /// \return Coordinate Info
    CoordinateInfo const& CoordInfo()             const { return cInfo_; }
    /// \return Number of expected frames to write.
    int NframesToWrite()                          const { return NframesToWrite_; }
    /// \return Number of frames that have been written.
    int NframesWritten()                          const { return numFramesWritten_; }
    /// \return Write format, can be changed.
    TrajectoryFile::TrajFormatType& WriteFormat()       { return writeFormat_; }
    /// \return Write format.
    TrajectoryFile::TrajFormatType  WriteFormat() const { return writeFormat_; }
    /// Set append status
    void SetAppend(bool a) { append_ = a; }
    /// Process common arguments
    int CommonTrajoutSetup(FileName const&, ArgList&, TrajectoryFile::TrajFormatType);
    /// Check if file can be appended to with given format.
    static int CheckAppendFormat(FileName const&, TrajectoryFile::TrajFormatType&);
    /// Set up CoordinateInfo based on given topology etc and current options.
    int SetupCoordInfo(Topology*, int, CoordinateInfo const&);
    /// Print common info to STDOUT
    void CommonInfo() const;
    /// \return 1 if set should not be written.
    inline int CheckFrameRange(int);
    /// \return 'true' if output frame range has been set up
    inline bool HasRange() const { return hasRange_ || !frameCount_.DefaultSettings(); }
  private:
    FileName trajName_; // FIXME: Save this here?
    Topology* trajParm_;// FIXME: Save this here?
    CoordinateInfo cInfo_;             ///< Coordinate info for trajectory
    int NframesToWrite_;               ///< Expected number of frames to be written.
    // Track frame numbers
    Range FrameRange_;                 ///< List of frame numbers to write.
    Range::const_iterator rangeframe_; ///< If frame range defined, this is next frame in range.
    ActionFrameCounter frameCount_;    ///< Hold start/stop/offset values
    int numFramesWritten_;           ///< Number of frames that have been written so far.
    // Trajout arguments
    TrajectoryFile::TrajFormatType writeFormat_;
    std::string title_;                ///< Output traj title.
    bool noBox_;                       ///< If true do not put box information in output traj
    bool noVel_;
    bool noTemp_;
    bool noTime_;
    bool noFrc_;
    bool noReps_;
    bool append_;                      ///< If true, append to this file.
    bool hasRange_;                    ///< If true a frame range is defined.
};
// ----- INLINE ROUTINES -------------------------------------------------------
int OutputTrajCommon::CheckFrameRange(int set) {
  if (hasRange_) {
    // If no more frames in the framerange, skip.
    if (rangeframe_ == FrameRange_.end()) return 1;
    // If this frame is not the next in the range, skip.
    if ( *rangeframe_ != set ) return 1;
    // This frame is next in the range. Advance FrameRange iterator.
    ++rangeframe_;
  } else {
    //printf("DEBUG: CheckFrameRange(%i)= %i\n", set, frameCount_.CheckFrameCounter(set));
    if (frameCount_.CheckFrameCounter( set )) return 1;
  }
  // Frame will be processed. Increment frame count.
  ++numFramesWritten_;
  return 0;
}
#endif
