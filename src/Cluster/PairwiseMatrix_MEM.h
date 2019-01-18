#ifndef INC_CLUSTER_PAIRWISE_MATRIX_MEM_H
#define INC_CLUSTER_PAIRWISE_MATRIX_MEM_H
#include "PairwiseMatrix.h"
#include "../Matrix.h"
namespace Cpptraj {
namespace Cluster {

class PairwiseMatrix_MEM : public PairwiseMatrix {
  public:
    PairwiseMatrix_MEM(Metric* m) : PairwiseMatrix(MEM, m) {}
    // -------------------------------------------
    double GetFdist(int f1, int f2) const { return Mat_.element(frameToMat_[f1], frameToMat_[f2]); }
    double Frame_Distance(int, int) const;
    int CacheDistances(Cframes const&);
    // -------------------------------------------
  protected:
    void SetElement(int col, int row, double val) { Mat_.setElement(col, row, val); }
  private:
    Matrix<float> Mat_;  ///< Hold cached distances
    Cframes frameToMat_; ///< Hold indices into Mat_ for all cached frames, -1 for non-cached.
};

} /* END namespace Cluster */
} /* END namespace Cpptraj */
#endif
