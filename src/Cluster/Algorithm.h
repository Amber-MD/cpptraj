#ifndef INC_CLUSTER_ALGORITHM_H
#define INC_CLUSTER_ALGORITHM_H
#include "PairwiseMatrix.h"
#include "List.h"
#include "../CpptrajFile.h"
namespace Cpptraj {
namespace Cluster {

/// Abstract base class for implementing clustering algorithms.
class Algorithm {
  public:
    enum AType { HIERAGGLO = 0, DBSCAN, DPEAKS, KMEANS, UNSPECIFIED };

    Algorithm(AType t) : debug_(0), type_(t) {}
    virtual ~Algorithm() {}
    /// Set up clustering algorithm
    virtual int Setup(ArgList&) = 0;
    /// Report details on algorithm setup.
    virtual void Info() const = 0;
    /// Report brief details on algorithm setup to a file.
    virtual void Results(CpptrajFile&) const = 0;
    /// Perform clustering on specified frames using given distance matrix.
    virtual int DoClustering(List&, Cframes const&, PairwiseMatrix const&) = 0;
    /// Report any timing data
    virtual void Timing(double) const = 0;
    /// /return Algorithm-specific between-cluster distance. Default to centroid distance.
    virtual double ClusterDistance(Node const&, Node const&, PairwiseMatrix const&,
                                   bool, Cframes const&) const;
    /// /return Epsilon for density-based algorithms; intended for use with sieve restore.
    virtual double Epsilon() const { return 0.0; }
    // -------------------------------------------
    /// Set debug level for algorithm
    void SetDebug(int d) { debug_ = d; }
    /// \return Algorithm type
    AType Type() const { return type_; }
  protected:
    int debug_;
  private:
    AType type_;
};

} /* END namespace Cluster */
} /* END namespace Cpptraj */
#endif
