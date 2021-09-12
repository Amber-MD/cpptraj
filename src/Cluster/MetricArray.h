#ifndef CPPTRAJ_CLUSTER_METRICARRAY_H
#define CPPTRAJ_CLUSTER_METRICARRAY_H
#include <vector>
namespace Cpptraj {
namespace Cluster {
class Metric;
/// Hold Metrics of various types
class MetricArray {
  public:
    MetricArray();
  private:
    std::vector<Metric*> metrics_;
};

}
}
#endif
