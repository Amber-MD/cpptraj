#ifndef INC_TRAJIN_H
#define INC_TRAJIN_H
#include "TrajectoryFile.h"
#include "ProgressBar.h"
/// Class that all input trajectories will inherit.
class Trajin : public TrajectoryFile {
  public:
    Trajin();
    virtual ~Trajin() {}
    virtual int SetupTrajRead(std::string const&, ArgList*, Topology *) = 0;
    virtual int BeginTraj(bool) = 0;
    virtual void EndTraj() = 0;
    virtual int GetNextFrame(Frame&) = 0;
    virtual void PrintInfo(int) const = 0;
    virtual bool HasVelocity() const = 0;
    virtual int NreplicaDimension() const = 0;

    int SetupTrajIO( std::string const&, TrajectoryIO&, ArgList* );
    int CheckBoxInfo(const char*, Box&, Box const&) const; 
    int setupFrameInfo();
    void SingleFrame();
    void PrepareForRead(bool,bool);
    void PrintInfoLine() const;
    void PrintFrameInfo() const;
 
    int TotalFrames()        const { return total_frames_;       }
    int TotalReadFrames()    const { return total_read_frames_;  }
    int CurrentFrame()       const { return currentFrame_;       }
    int Start()              const { return start_;              }
    int NumFramesProcessed() const { return numFramesProcessed_; }
    bool IsEnsemble()        const { return isEnsemble_;         }
    void SetEnsemble(bool b)       { isEnsemble_ = b;            }

    inline void SetTotalFrames(int); 
    inline bool CheckFinished();
    inline bool ProcessFrame();
  private:
    int start_;              ///< Frame to begin processing
    int stop_;               ///< Frame to end processing
    int offset_;             ///< Number of frames to skip between processed frames
    int total_frames_;       ///< The total number of frames in the traj
    int total_read_frames_;  ///< # frames that will actually be read based on start/stop/offset
    int currentFrame_;       ///< The current frame number being read
    int targetSet_;          ///< The next frame to process
    int frameskip_;          ///< The number of frames to skip over while reading
    int numFramesProcessed_; ///< Number of frames that have been read.
    ProgressBar progress_;   ///< Keep track of trajectory progress
    bool useProgress_;       ///< Indicate whether progress should be shown.
    bool isEnsemble_;        ///< True if this will be processed as an ensemble.
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
// Trajin::SetTotalFrames()
void Trajin::SetTotalFrames(int nfIn) { 
  total_frames_ = nfIn;
  if (stop_ > total_frames_) stop_ = total_frames_;
}
// Trajin::CheckFinished()
bool Trajin::CheckFinished() {
  if (currentFrame_ > stop_ && stop_ != -1) return true;
  if (useProgress_)
    progress_.Update(numFramesProcessed_);
  return false;
}
// Trajin::ProcessFrame()
bool Trajin::ProcessFrame() {
  bool tgtFrameFound = false;
  if (currentFrame_ == targetSet_) {
    tgtFrameFound=true;
    targetSet_ += offset_;
  }
  ++numFramesProcessed_;
  currentFrame_ += frameskip_;
  return tgtFrameFound;
}
#endif
