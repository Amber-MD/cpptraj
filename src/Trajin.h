#ifndef INC_TRAJIN_H
#define INC_TRAJIN_H
#include "TrajectoryFile.h"
#include "ProgressBar.h"
class Trajin : public TrajectoryFile {
  public:
    Trajin();
    virtual ~Trajin() {}
    virtual int SetupTrajRead(std::string const&, ArgList *, Topology *) = 0;
    virtual int BeginTraj(bool) = 0;
    virtual int EndTraj() = 0;
    virtual int GetNextFrame(Frame&) = 0;
    virtual void PrintInfo(int) = 0;
    virtual bool HasVelocity() = 0;

    int StartStopOffset( TrajectoryIO*, ArgList* ); 
    int setupFrameInfo();

    int TotalFrames() { return total_frames_; }

    void SetTotalFrames(int nfIn) { total_frames_ = nfIn; }
  private:
    int start_;              ///< Frame to begin processing
    int stop_;               ///< Frame to end processing
    int offset_;             ///< Number of frames to skip between processed frames
    int total_frames_;       ///< The total number of frames in the traj
    int total_read_frames_;  ///< # frames that will actually be read based on start/stop/offset
    int currentFrame_;       ///< The current frame number being read
    int targetSet_;          ///< The next frame to process
    int frameskip_;          ///< The number of frames to skip over while reading
    ProgressBar progress_;   ///< Keep track of trajectory progress
    bool useProgress_;       ///< Indicate whether progress should be shown.
};
#endif
