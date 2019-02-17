#ifndef INC_CLUSTER_OUTPUT_H
#define INC_CLUSTER_OUTPUT_H
#include "../CpptrajFile.h"
#include "List.h"
#include "Algorithm.h"
#include "Metric.h"
#include "Cframes.h"
#include "PairwiseMatrix.h" // TODO anything that needs this calcd outside here?

namespace Cpptraj {
namespace Cluster {

/// Cluster output routines.
class Output {
  public:
    static void PrintClustersToFile(CpptrajFile&, List const&, Algorithm const&, Metric*,
                                    int, Cframes const&);
    static int PrintSilhouetteFrames(CpptrajFile&, List const&);
    static int PrintSilhouettes(CpptrajFile&, List const&);
    static int Summary(CpptrajFile&, List const&, Algorithm const&, PairwiseMatrix const&,
                        bool, Cframes const&);
    static void Summary_Part(CpptrajFile&, unsigned int, Cframes const&, List const&);
  private:
    static unsigned int DetermineNameWidth(List const&);
};

}
}
#endif
