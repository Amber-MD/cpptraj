#ifndef INC_CLUSTER_OUTPUT_H
#define INC_CLUSTER_OUTPUT_H
#include "../CpptrajFile.h"
#include "List.h"
#include "Algorithm.h"
#include "Metric.h"
#include "Cframes.h"

namespace Cpptraj {
namespace Cluster {

/// Cluster output routines.
class Output {
  public:
    static void PrintClustersToFile(CpptrajFile&, List const&, Algorithm const&, Metric*,
                                    int, Cframes const&);
    static int PrintSilhouetteFrames(CpptrajFile&, std::vector< std::vector<double> > const&);
    static int PrintSilhouettes(CpptrajFile&, std::vector<double> const&);
};

}
}
#endif
