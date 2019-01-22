#ifndef INC_CLUSTER_CONTROL_H
#define INC_CLUSTER_CONTROL_H
#include "List.h"
#include "PairwiseMatrix.h"
#include "Algorithm.h"
#include "Metric.h"
#include "../DataSet_Coords.h"
namespace Cpptraj {
namespace Cluster {

/// Hold clusters, algorithm, and pairwise matrix.
class Control {
  public:
    Control();

    static const char* PairwiseArgs;
    static const char* AlgorithmArgs;

    enum SieveRestoreType { NO_RESTORE = 0, CLOSEST_CENTROID, EPSILON_CENTROID, EPSILON_FRAME };

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
};

} /** END namespace Cluster. */
} /** END namespace Cpptraj. */
#endif
