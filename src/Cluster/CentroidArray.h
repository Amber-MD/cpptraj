#ifndef CPPTRAJ_CLUSTER_CENTROIDARRAY_H
#define CPPTRAJ_CLUSTER_CENTROIDARRAY_H
#include <vector>
namespace Cpptraj {
namespace Cluster {
class Centroid;
/// Hold Centroids of various types
class CentroidArray {
  public:
    /// CONSTRUCTOR
    CentroidArray();
    /// DESTRUCTOR
    ~CentroidArray();
    /// COPY
    CentroidArray(CentroidArray const&);
    /// ASSIGN
    CentroidArray& operator=(CentroidArray const&);
    /// Clear the centroid array
    void Clear();
  private:
    std::vector<Centroid*> centroids_;
};

}
}
#endif
