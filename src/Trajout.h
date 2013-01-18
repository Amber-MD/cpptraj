#ifndef INC_TRAJOUT_H
#define INC_TRAJOUT_H
#include "TrajectoryFile.h"
#include "Range.h"
#include "ActionFrameCounter.h"
/// Output trajectory class.
class Trajout : public TrajectoryFile {
  public:
    Trajout();
    ~Trajout();
    /*virtual ~Trajout() {}
    virtual int SetupTrajWrite(std::string const&, ArgList*, Topology*, 
                               TrajectoryFile::TrajFormatType) = 0;
    virtual int EndTraj() = 0;
    virtual int WriteFrame(int, Topology*, Frame&) = 0;
    virtual void PrintInfo(int) = 0;*/
    int SetupTrajWrite(std::string const&, ArgList*, Topology*, TrajectoryFile::TrajFormatType);
    int SetupTrajWriteWithArgs(std::string const&, const char*, Topology*,
                               TrajectoryFile::TrajFormatType);
    void EndTraj();
    int WriteFrame(int, Topology*, Frame&);
    void PrintInfo(int);
    bool TrajIsOpen()        { return trajIsOpen_;         }
    int NumFramesProcessed() { return numFramesProcessed_; }
  private:
    int numFramesProcessed_;
    TrajectoryIO* trajio_;
    bool trajIsOpen_;                  ///< If true trajectory has been opened.
    bool nobox_;                       ///< If true do not put box information in output traj
    bool append_;                      ///< If true, append to this file.
    bool hasRange_;                    ///< If true a frame range is defined.
    Range FrameRange_;                 ///< List of frame numbers to write.
    Range::const_iterator rangeframe_; ///< If frame range defined, this is next frame in range.
    ActionFrameCounter frameCount_;    ///< Hold start/stop/offset values
};
#endif
