#ifndef INC_CLUSTER_CENTROID_COORD_H
#define INC_CLUSTER_CENTROID_COORD_H
#include "../Frame.h"
#include "Centroid.h"
namespace Cpptraj {
namespace Cluster {

/// Cluster Centroid for Coords DataSet.
class Centroid_Coord : public Centroid {
  public:
    Centroid_Coord() {}
    Centroid_Coord(Frame const& frame) : cframe_(frame) {}
    Centroid_Coord(int natom)          : cframe_(natom) {}
    Centroid* Copy() { return (Centroid*)new Centroid_Coord(cframe_); }
    void Print(std::string const&) const;
    //friend class ClusterDist_DME;
    //friend class ClusterDist_RMS;
    //friend class ClusterDist_SRMSD;
  private:
    Frame cframe_;
};


}
}
#endif
