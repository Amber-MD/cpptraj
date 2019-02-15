#include "Results_Coords.h"

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

