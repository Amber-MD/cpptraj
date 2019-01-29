#ifndef INC_CLUSTER_PAIRWISE_MATRIX_NC_H
#define INC_CLUSTER_PAIRWISE_MATRIX_NC_H
#include "PairwiseMatrix.h"
#include "../NC_Cmatrix.h"
namespace Cpptraj {
namespace Cluster {

/// Cache pairwise distances on disk via NetCDF
class PairwiseMatrix_NC : public PairwiseMatrix {
  public:
    PairwiseMatrix_NC()          : PairwiseMatrix(NC)    {}
    PairwiseMatrix_NC(Metric* m) : PairwiseMatrix(NC, m) {}
    // -------------------------------------------
    inline double GetFdist(int f1, int f2) const;
    double Frame_Distance(int, int) const;
    int CacheDistances(Cframes const&);
    void PrintCached() const;
    // -------------------------------------------
    void SetFilename(FileName const& f) { fname_ = f; }
  protected:
    void SetElement(int x, int y, double val) { file_.WriteCmatrixElement(x, y, val); }
  private:
    NC_Cmatrix file_;
    FileName fname_; // TODO should this be part of NC_Cmatrix?
};

double PairwiseMatrix_NC::GetFdist(int x, int y) const {
  return file_.GetCmatrixElement( frameToMat_[x], frameToMat_[y] );
}


}
}
#endif
