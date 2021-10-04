#include "Algorithm.h"
#include "MetricArray.h"
#include "Node.h"
#include "../CpptrajStdio.h"

/** \return Distance between given cluster centroids. */
double Cpptraj::Cluster::Algorithm::ClusterDistance(Node const& C1, Node const& C2,
                                                    MetricArray& pmatrix,
                                                    bool includeSieved,
                                                    Cframes const& sievedOut) const
{
  if (C1.Cent().empty() || C2.Cent().empty()) {
    mprinterr("Internal Error: One or both centroids are null in ClusterDistance().\n");
    return -1.0;
  }
  return pmatrix.CentroidDist( C1.Cent(), C2.Cent() );
}
