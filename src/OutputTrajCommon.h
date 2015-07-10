#ifndef INC_OUTPUTTRAJCOMMON_H
#define INC_OUTPUTTRAJCOMMON_Hi
#include "TrajectoryFile.h"
#include "Range.h"
#include "ActionFrameCounter.h"
/// Output trajectory/ensemble common functionality.
class OutputTrajCommon {
  public:
    OutputTrajCommon();
    FileName const& Filename() const { return trajName_; }
    // TODO: return const&, modify all TrajectoryIO routines?
    Topology* Parm()                  const { return trajParm_; }
    /// \return true if trajectory should be appended to.
    bool Append() const { return append_; }
    /// \return Write format, can be changed.
    TrajectoryFile::TrajFormatType& WriteFormat() { return writeFormat_; }
    /// \return Write format.
    TrajectoryFile::TrajFormatType const& WriteFormat() const { return writeFormat_; }
    /// Set append status
    void SetAppend(bool a) { append_ = a; }
    /// Process common arguments
    int CommonTrajoutSetup(std::string const&, ArgList&,
                           TrajectoryFile::TrajFormatType);
    /// Check if file can be appended to with given format.
    static int CheckAppendFormat(std::string const&, TrajectoryFile::TrajFormatType&);
  private:
    FileName trajName_;
    Topology* trajParm_;
    std::string title_;                ///< Output traj title.
    int numFramesProcessed_;           ///< Number of frames that have been written so far.
    // Trajout arguments
    TrajectoryFile::TrajFormatType writeFormat_;
    bool nobox_;                       ///< If true do not put box information in output traj
    bool append_;                      ///< If true, append to this file.
    bool hasRange_;                    ///< If true a frame range is defined.
    Range FrameRange_;                 ///< List of frame numbers to write.
    Range::const_iterator rangeframe_; ///< If frame range defined, this is next frame in range.
    ActionFrameCounter frameCount_;    ///< Hold start/stop/offset values
};
#endif 
