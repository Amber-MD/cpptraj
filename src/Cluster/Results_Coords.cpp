#include "Results_Coords.h"
#include "../CpptrajStdio.h"

/** The default output trajectory format. */
const TrajectoryFile::TrajFormatType Cpptraj::Cluster::Results_Coords::DEF_TRAJ_FMT_ =
  TrajectoryFile::AMBERTRAJ;


void Cpptraj::Cluster::Results_Coords::GetClusterTrajArgs(ArgList& argIn,
                                             const char* trajKey, const char* fmtKey,
                                             std::string& trajName,
                                             TrajectoryFile::TrajFormatType& fmt) const
{
  trajName = argIn.GetStringKey( trajKey );
  fmt = TrajectoryFile::WriteFormatFromString( argIn.GetStringKey(fmtKey), fmt );
  // If file name specified but not format, try to guess from name
  if (!trajName.empty() && fmt == TrajectoryFile::UNKNOWN_TRAJ)
    fmt = TrajectoryFile::WriteFormatFromFname( trajName, DEF_TRAJ_FMT_ );
}

int Cpptraj::Cluster::Results_Coords::GetOptions(ArgList& analyzeArgs) {
  writeRepFrameNum_ = analyzeArgs.hasKey("repframe");
  GetClusterTrajArgs(analyzeArgs, "clusterout",   "clusterfmt",   clusterfile_,   clusterfmt_);
  GetClusterTrajArgs(analyzeArgs, "singlerepout", "singlerepfmt", singlerepfile_, singlerepfmt_);
  GetClusterTrajArgs(analyzeArgs, "repout",       "repfmt",       reptrajfile_,   reptrajfmt_);
  GetClusterTrajArgs(analyzeArgs, "avgout",       "avgfmt",       avgfile_,       avgfmt_);

  return 0;
}

void Cpptraj::Cluster::Results_Coords::Info() const {
    if (!clusterfile_.empty())
    mprintf("\tCluster trajectories will be written to %s, format %s\n",
            clusterfile_.c_str(), TrajectoryFile::FormatString(clusterfmt_));
  if (!singlerepfile_.empty())
    mprintf("\tCluster representatives will be written to 1 traj (%s), format %s\n",
            singlerepfile_.c_str(), TrajectoryFile::FormatString(singlerepfmt_));
  if (!reptrajfile_.empty()) {
    mprintf("\tCluster representatives will be written to separate trajectories,\n");
    mprintf("\t\tprefix (%s), format %s",reptrajfile_.c_str(), 
            TrajectoryFile::FormatString(reptrajfmt_));
    if (writeRepFrameNum_) mprintf(", with frame #s");
    mprintf("\n");
  }
  if (!avgfile_.empty())
    mprintf("\tAverage structures for clusters will be written to %s, format %s\n",
            avgfile_.c_str(), TrajectoryFile::FormatString(avgfmt_));
}
