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
    Control() : metric_(0), pmatrix_(0), algorithm_(0) {}

    static void Help();

    int SetupForCoordsDataSet(DataSet_Coords*, std::string const&, ArgList&, int);
  private:
    static PairwiseMatrix* AllocatePairwise(PairwiseMatrix::Type, Metric*);
    int AllocatePairwise(ArgList&, Metric*);
    static Metric* AllocateMetric(Metric::Type);

    List clusters_;
    Metric* metric_;
    PairwiseMatrix* pmatrix_;
    Algorithm* algorithm_;
};

} /** END namespace Cluster. */
} /** END namespace Cpptraj. */
#endif
