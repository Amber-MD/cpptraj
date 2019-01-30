#include "DataSet_PairwiseCache_MEM.h"
#include "CpptrajStdio.h"

/** Print cached distances to stdout. */
void DataSet_PairwiseCache_MEM::PrintCached() const {
  for (Cframes::const_iterator it1 = FrameToMat().begin(); it1 != FrameToMat().end(); ++it1)
  {
    if (*it1 != -1) {
      for (Cframes::const_iterator it2 = it1 + 1; it2 != FrameToMat().end(); ++it2)
      {
        if (*it2 != -1)
          mprintf("\t%i %i %f\n", *it1+1, *it2+1, Mat_.element(*it1, *it2));
      }
    }
  }
}
