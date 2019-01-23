#ifndef INC_CLUSTER_BESTREPS_H
#define INC_CLUSTER_BESTREPS_H
#include <map>
#include "Node.h"
#include "List.h"
#include "PairwiseMatrix.h"
namespace Cpptraj {
namespace Cluster {

/// Used to find best representative structures for a cluster.
class BestReps {
  public:
    enum RepMethodType { NO_REPS = 0, CUMULATIVE, CENTROID, CUMULATIVE_NOSIEVE };

    static int FindBestRepFrames(RepMethodType, int, List&, PairwiseMatrix const&,
                                 Cframes const&, int);
  private:
    /// Used to pair representative score with frame number.
    typedef std::pair<double, int> RepPair;
    /// Used to hold pairs of representative scores to frames.
    typedef std::multimap<double, int> RepMap;

    /// Save up to maxSize of the best (lowest) representative scores/frames.
    static void SaveBestRep(RepMap&, RepPair const&, unsigned int);
    /// Set given cluster node with best representative frames/scores in reps
    static void SetBestRepFrame(Node& node, RepMap const&);
    /// Find best representative frames by shortest distance to all other frames.
    static int FindBestRepFrames_CumulativeDist(int, List&, PairwiseMatrix const&);
    /// Find best representative frames by shortest distance, ignoring sieved frames.
    static int FindBestRepFrames_NoSieve_CumulativeDist(int, List&, PairwiseMatrix const&,
                                                        Cframes const&);
    /// Find best representative frames by shortest distance to centroid.
    static int FindBestRepFrames_Centroid(int, List&, PairwiseMatrix const&);
};

}
}
#endif
