#ifndef INC_CLUSTER_ALGORITHM_NONE_H
#define INC_CLUSTER_ALGORITHM_NONE_H
#include "Algorithm.h"
namespace Cpptraj {
namespace Cluster {
/// No clustering. Intended for use when reading in previous clusters.
class Algorithm_None : public Algorithm {
  public:
    Algorithm_None();
    static void Help();
    int Setup(ArgList&);
    void Info() const;
    void Results(CpptrajFile&) const;
    int DoClustering(List&, Cframes const&, MetricArray&);
    void Timing(double) const {}
};
}
}
#endif
