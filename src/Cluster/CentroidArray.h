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

    /// \return number of centroids
    unsigned int size() const { return centroids_.size(); }
    /// \return true if no centroids
    bool empty() const { return centroids_.empty(); }

    /// Add given centroid pointer to array
    void push_back( Centroid* c) { centroids_.push_back( c ); }
    /// \return Centroid pointer at given position
    Centroid* operator[](int idx) const { return centroids_[idx]; }
  private:
    std::vector<Centroid*> centroids_;
};

}
}
#endif
