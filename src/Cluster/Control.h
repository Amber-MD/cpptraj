#ifndef INC_CLUSTER_CONTROL_H
#define INC_CLUSTER_CONTROL_H
#include "List.h"
#include "PairwiseMatrix.h"
#include "Algorithm.h"
#include "Metric.h"
#include "BestReps.h"
#include "../DataSet_Coords.h"
namespace Cpptraj {
namespace Cluster {

/// Hold clusters, algorithm, and pairwise matrix.
class Control {
  public:
    Control();
    ~Control();

    static const char* PairwiseArgs;
    static const char* AlgorithmArgs;

    enum SieveRestoreType { NO_RESTORE = 0, CLOSEST_CENTROID, EPSILON_CENTROID, EPSILON_FRAME };
    enum SieveType { NONE=0, REGULAR, RANDOM };

    static int SetFramesToCluster(Cframes&, Cframes&, int, std::size_t, int);

    int SetupForCoordsDataSet(DataSet_Coords*, std::string const&, ArgList&, int);

    void Info() const;
    int Run();

    List const& Clusters()     const { return clusters_; }
    Metric const& DistMetric() const { return *metric_;  }
  private:
    static PairwiseMatrix* AllocatePairwise(PairwiseMatrix::Type, Metric*);
    int AllocatePairwise(ArgList&, Metric*);

    static Metric* AllocateMetric(Metric::Type);

    static Algorithm* AllocateAlgorithm(Algorithm::AType);
    int AllocateAlgorithm(ArgList&);    

    int Common(ArgList&);

    List clusters_;
    Metric* metric_;
    PairwiseMatrix* pmatrix_;
    Algorithm* algorithm_;
    int verbose_;

    int sieve_;                     ///< Sieve value
    int sieveSeed_;                 ///< Seed if doing random sieve
    SieveRestoreType sieveRestore_; ///< Specify if/how sieved frames should be restored.
    double restoreEpsilon_;         ///< Cutoff to use if restoring by epsilon to centroids.

    BestReps::RepMethodType bestRep_; ///< How to determine best rep frames.
    int nRepsToSave_;                 ///< How many rep frames to save.

    bool suppressInfo_;               ///< If true do not write cluster info
    std::string clusterinfo_;         ///< Cluster info file name.

};

} /** END namespace Cluster. */
} /** END namespace Cpptraj. */
#endif
