#ifndef INC_CLUSTER_ALGORITHM_H
#define INC_CLUSTER_ALGORITHM_H
#include "PairwiseMatrix.h"
#include "List.h"
namespace Cpptraj {
namespace Cluster {

/// Abstract base class for implementing clustering algorithms.
class Algorithm {
  public:
    Algorithm() {}
    virtual ~Algorithm() {}
    /// Set up clustering algorithm
    virtual int Setup(ArgList&) = 0;
    /// Report details on algorithm setup.
    virtual void Info() const = 0;
    /// Report details on algorithm setup as a string.
    virtual void Results(std::string&) const = 0;
    /// Perform clustering on specified frames using given distance matrix.
    virtual int DoClustering(List&, Cframes const&, PairwiseMatrix const&) = 0;
    /// Report any timing data
    virtual void Timing(double) const = 0;
};

} /* END namespace Cluster */
} /* END namespace Cpptraj */
#endif
