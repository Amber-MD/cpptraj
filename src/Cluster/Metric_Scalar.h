#ifndef INC_CPPTRAJ_CLUSTER_METRIC_SCALAR_H
#define INC_CPPTRAJ_CLUSTER_METRIC_SCALAR_H
#include "Metric.h"
class DataSet_1D;
namespace Cpptraj {
namespace Cluster {
/// Metric for delta between points in a regular 1D scalar DataSet.
class Metric_Scalar : public Metric {
  public:
    /// CONSTRUCTOR
    Metric_Scalar();
    /// Initialize with DataSet
    int Init(DataSet_1D*);

    // ----- Metric Routines ---------------------
    /// Setup
    int Setup();
    /// \return Copy of this Metric TODO const?
    Metric* Copy() { return new Metric_Scalar( *this ); }

    /// \return Absolute difference between two values indicated by given indices
    double FrameDist(int, int);
    /// \return Absolute difference between two values of given Centroids (averages)
    double CentroidDist(Centroid*, Centroid*);
    /// \return Absolute difference between value indicated by index and value of centroid (average)
    double FrameCentroidDist(int, Centroid*);

    /// Calculate centroid (average) from given frames.
    void CalculateCentroid(Centroid*, Cframes const&);
    /// \return new centroid from given frames.
    Centroid* NewCentroid(Cframes const&);
    /// Update centroid (average) by performing given operation between given frame and centroid.
    void FrameOpCentroid(int, Centroid*, double, CentOpType);

    /// \return string containing description of the distance metric
    std::string Description() const;
    /// Print Metric info to stdout.
    void Info() const;
    /// \return total number of frames.
    unsigned int Ntotal() const;
  private:
    DataSet_1D* data_; ///< DataSet to calculate distances for.
};

}
}
#endif   
