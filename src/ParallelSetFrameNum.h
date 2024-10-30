#ifndef INC_PARALLELSETFRAMENUM_H
#define INC_PARALLELSETFRAMENUM_H
#ifdef MPI
# include "Parallel.h"
#endif
class DataSet;
namespace Cpptraj {
/// Used to pick the correct frame from a data set in parallel
/** In parallel, Actions that access data from DataSets need to handle
  * two cases:
  *   1) The data set exists already before trajectory processing.
  *   2) The data set is being created during trajectory processing.
  * In case 1, the data set should be accessed using ActionFrame::TrajoutNum().
  * In case 2, the data set should be accessed using the frameNum variable
  * since the set has not yet been synced.
  * The SetForParallel() routine will also ensure the data set is broadcast
  * to all processes if it needs to be.
  */
class ParallelSetFrameNum {
  public:
    /// CONSTRUCTOR
    ParallelSetFrameNum() : set_(0), exists_(false) {}
#   ifdef MPI
    void SetComm(Parallel::Comm const& commIn) { trajComm_ = commIn; }
#   endif
    /// Initialize with DataSet, record if it has any data yet.
    int SetForParallel(DataSet*);
    /// \return Correct frame number in parallel based on whether set exists or is being generated
    int ParallelFrameNum(int frameNum, int trajoutNum) const {
#     ifdef MPI
      // In parallel, number to use depends on whether the set is being generated or not.
      if (exists_)
        return trajoutNum;
      else
        return frameNum;
#     else
      // In serial, just return the frame number.
      return frameNum;
#     endif
    }
  private:
    DataSet* set_; ///< The DataSet in question
    bool exists_;  ///< True if set exists, false if it is being generated.
#   ifdef MPI
    Parallel::Comm trajComm_; ///< The communicator associated with the set
#   endif
};
} // END namespace Cpptraj
#endif
