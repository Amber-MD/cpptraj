#ifndef INC_CLUSTER_SIEVE_H
#define INC_CLUSTER_SIEVE_H
#include "Cframes.h"
namespace Cpptraj {
namespace Cluster {

/// Hold frames to be clustered and frames to be sieved out.
class Sieve {
  public:
    enum SieveType { NONE=0, REGULAR, RANDOM };

    Sieve() : type_(NONE), sieve_(1) {}

    int SetFramesToCluster(int, std::size_t, int);

    Cframes const& FramesToCluster() const { return framesToCluster_; }
    Cframes const& SievedOut()       const { return sievedOut_;       }
    int SieveValue()                 const { return sieve_;           }
  private:
    void DetermineTypeFromSieve(int);

    Cframes framesToCluster_; ///< Frames to cluster.
    Cframes sievedOut_;       ///< Frames that will not be clustered (i.e. sieved out).
    SieveType type_;          ///< Sieveing type.
    int sieve_;               ///< Sieving value.
};

}
}
#endif
