#include "ParallelSetFrameNum.h"
#include "CpptrajStdio.h"
#include "DataSet.h"
using namespace Cpptraj;

void ParallelSetFrameNum::SetForParallel(DataSet const* setIn) {
  set_ = setIn;
# ifdef MPI
  // In parallel, number to use will depend on whether the set exists or not
  if (set_->Size() > 0)
    exists_ = true;
  else
    exists_ = false;
  if (!exists_)
    rprintf("DEBUG: Set '%s' does not yet exist.\n", set_->legend());
  else
    rprintf("DEBUG: Set '%s' exists with size %zu.\n", set_->legend(), set_->Size());
# endif /*MPI */
}
