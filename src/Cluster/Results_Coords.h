#ifndef INC_CLUSTER_RESULTS_COORDS_H
#define INC_CLUSTER_RESULTS_COORDS_H
#include "Results.h"
#include "../TrajectoryFile.h"
class DataSet_Coords;
namespace Cpptraj {
namespace Cluster {

/// Class for handling cluster results specific to COORDS input data.
class Results_Coords : public Results {
  public:
    Results_Coords(DataSet_Coords*);
    // ----- Results functions -------------------
    int GetOptions(ArgList&);
    void Info() const;
    int DoOutput(List const&) const;
  private:
    void GetClusterTrajArgs(ArgList&, const char*, const char*, std::string&,
                            TrajectoryFile::TrajFormatType&) const;
    void WriteClusterTraj( List const& ) const;
    void WriteAvgStruct( List const& ) const;
    void WriteSingleRepTraj( List const& ) const;
    void WriteRepTraj( List const& ) const;

    static const TrajectoryFile::TrajFormatType DEF_TRAJ_FMT_;

    DataSet_Coords* coords_;                      ///< The COORDS data set
    bool writeRepFrameNum_;                       ///< If true frame #s will be in rep file names.
    TrajectoryFile::TrajFormatType clusterfmt_;   ///< Cluster trajectory format.
    TrajectoryFile::TrajFormatType singlerepfmt_; ///< Cluster all rep single trajectory format.
    TrajectoryFile::TrajFormatType reptrajfmt_;   ///< Cluster rep to separate trajectory format.
    TrajectoryFile::TrajFormatType avgfmt_;       ///< Cluster traj average structure file format.
    std::string clusterfile_;   ///< Cluster trajectory base filename.
    std::string singlerepfile_; ///< Cluster all rep single trajectory filename.
    std::string reptrajfile_;   ///< Cluster rep to separate trajectory filename.
    std::string avgfile_;       ///< Cluster traj average structure filename.

};

}
}
#endif