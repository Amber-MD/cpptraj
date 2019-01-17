#ifndef INC_CLUSTER_METRIC_H
#define INC_CLUSTER_METRIC_H
#include <vector>
#include "../DataSet.h"
#include "Centroid.h"
namespace Cpptraj {
namespace Cluster {

/// This will hold cluster frame numbers.
typedef std::vector<int> Cframes;
/// Iterator for Cframes
typedef Cframes::const_iterator Cframes_it;
/// Definition of a noise point.
const int NOISE = -1;

/// Abstract base class for calculating distance between points or determining centroid.
class Metric {
  public:
    enum CentOpType { ADDFRAME=0, SUBTRACTFRAME };
    typedef std::vector<DataSet*> DsArray; // TODO should this be here?
    virtual ~Metric() {}
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
    /// \return total number of frames.
    virtual unsigned int Ntotal() const = 0;
  protected:
    typedef double (*DistCalc)(double,double);
};

} // END namespace Cluster
} // END namepsace Cpptraj
#endif
