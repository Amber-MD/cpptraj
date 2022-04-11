#ifndef INC_CLUSTER_METRIC_H
#define INC_CLUSTER_METRIC_H
#include <string>
#include <vector>
namespace Cpptraj {
namespace Cluster {
class Centroid;
class Cframes;
/// Definition of a noise point.
const int NOISE = -1;
/// Definition of an unclassified point.
const int UNCLASSIFIED = -2;

/// Abstract base class for calculating distance between points or determining centroid.
class Metric {
  public:
    enum Type { RMS=0, DME, SRMSD, SCALAR, TORSION, MATRIX2D, UNKNOWN_METRIC };
    enum CentOpType { ADDFRAME=0, SUBTRACTFRAME };

    /// CONSTRUCTOR
    Metric(Type t) { type_ = t; }
    /// DESTRUCTOR
    virtual ~Metric() {}
    /// Set up the Metric prior to clustering.
    virtual int Setup() = 0;
    /// \return distance between given frames.
    virtual double FrameDist(int, int) = 0;
    /// \return distance between given centroids.
    virtual double CentroidDist( Centroid*, Centroid* ) = 0;
    /// \return distance between given frame and centroid.
    virtual double FrameCentroidDist(int, Centroid* ) = 0;
    /// Calculate centroid from given frames.
    virtual void CalculateCentroid(Centroid*, Cframes const&) = 0;
    /// \return new centroid from given frames.
    virtual Centroid* NewCentroid(Cframes const&) = 0;
    /// \return copy of this Metric
    virtual Metric* Copy() = 0;
    /// Update centroid by performing given operation between given frame and centroid.
    virtual void FrameOpCentroid(int, Centroid*, double, CentOpType) = 0;
    /// \return string containing description of the distance metric
    virtual std::string Description() const = 0;
    /// Print Metric info to stdout.
    virtual void Info() const = 0;
    /// \return total number of frames.
    virtual unsigned int Ntotal() const = 0;

    /// \return Metric type
    Type MetricType() const { return type_; }
  private:
    Type type_;
};

} // END namespace Cluster
} // END namepsace Cpptraj
#endif
