#ifndef INC_CLUSTER_BESTREPS_H
#define INC_CLUSTER_BESTREPS_H
#include <map>
namespace Cpptraj {
namespace Cluster {
class Cframes;
class List;
class MetricArray;
class Node;
/// Used to find best representative structures for a cluster.
class BestReps {
  public:
    enum RepMethodType { NO_REPS = 0, CUMULATIVE, CENTROID, CUMULATIVE_NOSIEVE };

    /// CONSTRUCTOR
    BestReps();

    /// Initialize best rep frames search with method type, # to save, and debug level
    int InitBestReps(RepMethodType, int, int);
    /// Find best rep frames for each cluster in given list
    int FindBestRepFrames(List&, MetricArray&, Cframes const&) const;
    /// Find the best rep frames for given node
    int FindBestRepFrames(Node&, MetricArray&, Cframes const&) const;
  private:
    /// Print best reps to stdout
    static void PrintBestReps(Node const&);

    /// Used to pair representative score with frame number.
    typedef std::pair<double, int> RepPair;
    /// Used to hold pairs of representative scores to frames.
    typedef std::multimap<double, int> RepMap;

    /// Save up to maxSize of the best (lowest) representative scores/frames.
    static inline void SaveBestRep(RepMap&, RepPair const&, unsigned int);
    /// Set given cluster node with best representative frames/scores in reps
    static inline void SetBestRepFrame(Node& node, RepMap const&);

    /// Find best rep frames in node by shortest distance to all other frames.
    int FindBestRepFrames_CumulativeDist(Node&, MetricArray&) const;
    /// Find best representative frames by shortest distance to all other frames.
    int FindBestRepFrames_CumulativeDist(List&, MetricArray&) const;
    /// Find best rep frames in node by shortest distance, ignoring sieved frames.
    int FindBestRepFrames_NoSieve_CumulativeDist(Node&, MetricArray&, Cframes const&) const;
    /// Find best representative frames by shortest distance, ignoring sieved frames.
    int FindBestRepFrames_NoSieve_CumulativeDist(List&, MetricArray&, Cframes const&) const;
    /// Find best rep frames in node by shortest distance to centroid.
    int FindBestRepFrames_Centroid(Node&, MetricArray&) const;
    /// Find best representative frames by shortest distance to centroid.
    int FindBestRepFrames_Centroid(List&, MetricArray&) const;

    int debug_;   ///< Debug level, set in call to FindBestRepFrames
    int nToSave_; ///< Number of representatives to find
    RepMethodType type_; ///< Method to use to find best reps
};

}
}
#endif
