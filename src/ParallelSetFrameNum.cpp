#include "ParallelSetFrameNum.h"
#include "CpptrajStdio.h"
#include "DataSet.h"
using namespace Cpptraj;

int ParallelSetFrameNum::SetForParallel(DataSet* setIn) {
  set_ = setIn;
# ifdef MPI
  // In parallel, number to use will depend on whether the set exists or not.
  // Also check if the set needs to be synced. If it does, then it really 
  // does not exist yet but will be generated.
  if (set_->Size() == 0 && set_->NeedsSync())
    exists_ = false;
  else
    exists_ = true;
  // In parallel if the set exists it must be the same on all processes.
  if (exists_) {
    if (trajComm_.IsNull()) {
      mprinterr("Internal Error: ParallelSetFrameNum::SetForParallel(): Null communicator.\n");
      return 1;
    }
    long int mysize = set_->Size();
    std::vector<long int> sizeOnRank( trajComm_.Size() );
    trajComm_.AllGather( &mysize, 1, MPI_LONG, &sizeOnRank[0] );
    bool set_needs_bcast = false;
    for (int rank = 1; rank < trajComm_.Size(); rank++) {
      if (sizeOnRank[rank] != sizeOnRank[0]) {
#       ifdef DEBUG_CPPTRAJ_PARALLELSETFRAMENUM
        mprintf("DEBUG: Size on rank (%li) does not match rank 0 (%li)\n", sizeOnRank[rank], sizeOnRank[0]);
#       endif
        set_needs_bcast = true;
        break;
      }
    }
    if (set_needs_bcast) {
      int err = set_->Bcast(trajComm_);
      if (trajComm_.CheckError(err)) {
        mprinterr("Internal Error: ParallelSetFrameNum::SetForParallel(): DataSet broadcast error.\n");
        return 1;
      }
    }
  }
# ifdef DEBUG_CPPTRAJ_PARALLELSETFRAMENUM
  if (!exists_)
    rprintf("DEBUG: Set '%s' does not yet exist (needs sync %i).\n", set_->legend(), (int)set_->NeedsSync());
  else
    rprintf("DEBUG: Set '%s' exists with size %zu (needs sync %i).\n", set_->legend(), set_->Size(), (int)set_->NeedsSync());
# endif
# endif /*MPI */
  return 0;
}
