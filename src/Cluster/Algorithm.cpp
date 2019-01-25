#include "Algorithm.h"
#include "../CpptrajStdio.h"

double Cpptraj::Cluster::Algorithm::ClusterDistance(Node const& C1, Node const& C2,
                                                    PairwiseMatrix const& pmatrix,
                                                    bool includeSieved,
                                                    Cframes const& sievedOut) const
{
  if (C1.Cent() == 0 || C2.Cent() == 0) {
    mprinterr("Internal Error: One or both centroids are null in ClusterDistance().\n");
    return -1.0;
  }
  return pmatrix.MetricPtr()->CentroidDist( C1.Cent(), C2.Cent() );
}
