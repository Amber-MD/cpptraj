#ifndef CPPTRAJ_CLUSTER_CENTROIDARRAY_H
#define CPPTRAJ_CLUSTER_CENTROIDARRAY_H
#include <vector>
namespace Cpptraj {
namespace Cluster {
class Centroid;
/// Hold Centroids of various types
class CentroidArray {
  public:
    CentroidArray();
  private:
    std::vector<Centroid*> centroids_;
};

}
}
#endif
