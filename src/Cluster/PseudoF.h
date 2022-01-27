#ifndef INC_CLUSTER_PSEUDOF_H
#define INC_CLUSTER_PSEUDOF_H
namespace Cpptraj {
namespace Cluster {
class List;
class MetricArray;
/// Calculate pseudo-F
double ComputePseudoF(List const&, double&, MetricArray&, int);
}
}
#endif
