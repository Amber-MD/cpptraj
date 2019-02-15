#ifndef INC_CLUSTER_RESULTS_COORDS_H
#define INC_CLUSTER_RESULTS_COORDS_H
#include "Results.h"
#include "../TrajectoryFile.h"
namespace Cpptraj {
namespace Cluster {

/// Class for handling cluster results specific to COORDS input data.
class Results_Coords : public Results {
  public:
    Results_Coords() : Results(COORDS) {}
    // ----- Results functions -------------------
    int GetOptions(ArgList&);
    int DoOutput() const;
  private:
    void GetClusterTrajArgs(ArgList&, const char*, const char*, std::string&,
                            TrajectoryFile::TrajFormatType&) const;

    static const TrajectoryFile::TrajFormatType DEF_TRAJ_FMT_;
};

}
}
#endif
