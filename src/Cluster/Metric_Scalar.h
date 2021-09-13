#ifndef INC_CPPTRAJ_CLUSTER_METRIC_SCALAR_H
#define INC_CPPTRAJ_CLUSTER_METRIC_SCALAR_H
#include "Metric.h"
class DataSet_1D;
namespace Cpptraj {
namespace Cluster {
/// Metric for regular 1D scalar DataSet
class Metric_Scalar : public Metric {
  public:
    /// CONSTRUCTOR
    Metric_Scalar();
    /// Setup
    int Setup();
    /// \return Absolute difference between two values indicated by given indices
    double FrameDist(int, int);
    /// \return Absolute difference between two values of given Centroids (averages)
    double CentroidDist(Centroid*, Centroid*);
    /// \return Absolute difference between value indicated by index and value of centroid (average)

    
    
