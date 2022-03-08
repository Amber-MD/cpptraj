#ifndef INC_CLUSTER_CENTROID_H
#define INC_CLUSTER_CENTROID_H
#include <string>
namespace Cpptraj {
namespace Cluster {

/// Abstract Base Class for Cluster centroid.
/** This class is a container for the cluster centroid type appropriate for
  * the data being clustered. For COORDS DataSets this is a frame, for other
  * DataSets this is just a number. Centroid classes must implement a Copy()
  * function.
  */
class Centroid { 
  public:
    virtual ~Centroid() {}
    virtual Centroid* Copy() = 0;
  // TODO: Should centroids remember number of frames that went into them?
  //       This would make it so FrameOpCentroid wouldnt require extra arg.
    virtual void Print(std::string const&) const {}
};


} // END namespace Cluster
} // END namepsace Cpptraj
#endif
