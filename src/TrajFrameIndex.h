#ifndef INC_TRAJFRAMEINDEX_H
#define INC_TRAJFRAMEINDEX_H
#include <vector>
/// Class used to find absolute position inside a group of trajectories.
class TrajFrameIndex {
  public:
    TrajFrameIndex() : currentTrajNum_(-1), maxFrames_(0), trajHasChanged_(false) {}
    /// \return Internal index for current traj given a global index
    inline int FindIndex(int);
    inline void AddTraj(int, int, int); 
    inline int CurrentTrajNum()  const { return currentTrajNum_; }
    inline int MaxFrames()       const { return maxFrames_;      }
    inline bool TrajHasChanged() const { return trajHasChanged_; }
    inline size_t DataSize()     const;
  private:
    typedef std::vector<int> Iarray;
    Iarray TotalReadFrames_; ///< Total number of read frames in each trajectory.
    Iarray Starts_;          ///< Start argument for each trajectory.
    Iarray Offsets_;         ///< Offset argument for each trajectory.
    int currentTrajNum_;     ///< # of currently open input trajectory.
    int maxFrames_;          ///< Total # of readable frames in all trajectories.
    bool trajHasChanged_;    ///< True if new traj opened after last call to FindIndex()
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
// TrajFrameIndex::FindIndex()
int TrajFrameIndex::FindIndex(int idx) {
  // Determine which trajectory has the desired index.
  int globalOffset = 0; // Internal offset for converting global index to traj index.
  int currentMax = 0;   // Index after which we are in next trajectory.
  int desiredTrajNum = 0;
  for (; desiredTrajNum < (int)TotalReadFrames_.size(); ++desiredTrajNum) {
    currentMax += TotalReadFrames_[desiredTrajNum];
    if (idx < currentMax) break;
    globalOffset += TotalReadFrames_[desiredTrajNum];
  }
  if (desiredTrajNum == (int)TotalReadFrames_.size()) return -1;
  trajHasChanged_ = (desiredTrajNum != currentTrajNum_);
  currentTrajNum_ = desiredTrajNum;
  // Convert desired index into trajectory internal index.
  int internalIdx = ((idx - globalOffset) * Offsets_[currentTrajNum_]) +
                    Starts_[currentTrajNum_];
  return internalIdx;
}
// TrajFrameIndex::AddTraj()
void TrajFrameIndex::AddTraj(int total, int start, int offset) {
  TotalReadFrames_.push_back( total );
  maxFrames_ += total;
  Starts_.push_back( start );
  Offsets_.push_back( offset );
}
// TrajFrameIndex::DataSize()
size_t TrajFrameIndex::DataSize() const {
  return (TotalReadFrames_.size() * sizeof(int)) +
         (Starts_.size() * sizeof(int)) +
         (Offsets_.size() * sizeof(int)) +
         (2*sizeof(int) + sizeof(bool));
}
#endif
