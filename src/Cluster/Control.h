#ifndef INC_CLUSTER_CONTROL_H
#define INC_CLUSTER_CONTROL_H
#include "List.h"
#include "PairwiseMatrix.h"
#include "Algorithm.h"
#include "../DataSet_Coords.h"
namespace Cpptraj {
namespace Cluster {

/// Hold clusters, algorithm, and pairwise matrix.
class Control {
  public:
    Control() : pmatrix_(0), algorithm_(0) {}

    static void Help();

    int SetupForCoordsDataSet(DataSet_Coords*, ArgList&, int);
  private:
    List clusters_;
    PairwiseMatrix* pmatrix_;
    Algorithm* algorithm_;
};

} /** END namespace Cluster. */
} /** END namespace Cpptraj. */
#endif
