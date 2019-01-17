#include "PairwiseMatrix_MEM.h"

double Cpptraj::Cluster::PairwiseMatrix_MEM::Frame_Distance(int f1, int f2) const {
  int idx1 = frameToMat_[f1];
  if (idx1 != -1) {
    int idx2 = frameToMat_[f2];
    if (idx2 != -1) {
      return Mat_.element(idx1, idx2);
    }
  }
  // If here, distance was not cached.
  return distMetric()->FrameDist(f1, f2);
}
