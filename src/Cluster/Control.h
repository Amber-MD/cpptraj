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
    Control() : metric_(0), pmatrix_(0), algorithm_(0), verbose_(0),
                sieve_(1), sieveSeed_(-1)
      {}

    static const char* PairwiseArgs;
    static const char* AlgorithmArgs;

    int SetupForCoordsDataSet(DataSet_Coords*, std::string const&, ArgList&, int);

    void Info() const;
    int Run();

    List const& Clusters()     const { return clusters_; }
    Metric const& DistMetric() const { return *metric_;  }
  private:
    static PairwiseMatrix* AllocatePairwise(PairwiseMatrix::Type, Metric*);
    int AllocatePairwise(ArgList&, Metric*);

    static Metric* AllocateMetric(Metric::Type);

    static Algorithm* AllocateAlgorithm(Algorithm::Type);
    int AllocateAlgorithm(ArgList&);    

    int Common(ArgList&);

    List clusters_;
    Metric* metric_;
    PairwiseMatrix* pmatrix_;
    Algorithm* algorithm_;
    int verbose_;

    int sieve_; ///< Sieve value
    int sieveSeed_; ///< Seed if doing random sieve
};

} /** END namespace Cluster. */
} /** END namespace Cpptraj. */
#endif
