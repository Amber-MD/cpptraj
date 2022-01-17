#ifndef INC_CLUSTER_SIEVE_H
#define INC_CLUSTER_SIEVE_H
#include "Cframes.h"
class DataSet_PairwiseCache;
namespace Cpptraj {
namespace Cluster {

/// Hold frames to be clustered and frames to be sieved out.
class Sieve {
  public:
    enum SieveType { NONE=0, REGULAR, RANDOM };
    /// CONSTRUCTOR - default no sieve
    Sieve();
    /// Set which frames to cluster based on given sieve value, max# frames, and seed
    int SetFramesToCluster(int, std::size_t, int);
    /// Set which frames to cluster based on whats in the given pairwise cache and max# frames
    int SetupFromCache( DataSet_PairwiseCache const&, std::size_t );
    /// Clear the Sieve
    void Clear();
    /// \return Array of frames to cluster (i.e. are present)
    Cframes const& FramesToCluster() const { return framesToCluster_; }
    /// \return Array of frames sieved out (i.e. not present)
    Cframes const& SievedOut()       const { return sievedOut_;       }
    /// \return Sieve value used to set up sieve
    int SieveValue()                 const { return sieve_;           }
    /// \return An array containing true for present frames, false otherwise
    std::vector<bool> GenerateFrameIsPresentArray() const;
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
