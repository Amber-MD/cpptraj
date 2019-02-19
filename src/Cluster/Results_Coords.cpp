#include "Results_Coords.h"
#include "../CpptrajStdio.h"
#include "../StringRoutines.h"
#include "../Trajout_Single.h"

/** The default output trajectory format. */
const TrajectoryFile::TrajFormatType Cpptraj::Cluster::Results_Coords::DEF_TRAJ_FMT_ =
  TrajectoryFile::AMBERTRAJ;

/// CONSTRUCTOR
Cpptraj::Cluster::Results_Coords::Results_Coords(DataSet_Coords* ds) :
  Results(COORDS),
  coords_(ds),
  writeRepFrameNum_(false),
  clusterfmt_(TrajectoryFile::UNKNOWN_TRAJ),
  singlerepfmt_(TrajectoryFile::UNKNOWN_TRAJ),
  reptrajfmt_(TrajectoryFile::UNKNOWN_TRAJ),
  avgfmt_(TrajectoryFile::UNKNOWN_TRAJ)
{}


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

int Cpptraj::Cluster::Results_Coords::DoOutput(List const& CList) const {
  // Write clusters to trajectories
  if (!clusterfile_.empty())
    WriteClusterTraj( CList ); 
  // Write all representative frames to a single traj
  if (!singlerepfile_.empty())
    WriteSingleRepTraj( CList );
  // Write all representative frames to separate trajs
  if (!reptrajfile_.empty())
    WriteRepTraj( CList );
  // Write average structures for each cluster to separate files.
  if (!avgfile_.empty())
    WriteAvgStruct( CList );
  return 0;
}

/** Write frames in each cluster to a trajectory file.  */
void Cpptraj::Cluster::Results_Coords::WriteClusterTraj( List const& CList ) const {
  Topology* clusterparm = coords_->TopPtr();
  // Loop over all clusters
  for (List::cluster_iterator C = CList.begincluster();
                                     C != CList.endcluster(); ++C)
  {
    // Create filename based on cluster number.
    int cnum = C->Num();
    std::string cfilename =  clusterfile_ + ".c" + integerToString( cnum );
    // Set up trajectory file 
    Trajout_Single clusterout;
    if (clusterout.PrepareTrajWrite(cfilename, ArgList(), clusterparm,
                                    coords_->CoordsInfo(), C->Nframes(),
                                    clusterfmt_)) 
    {
      mprinterr("Error: Could not set up cluster trajectory %s for write.\n",
                cfilename.c_str());
      return;
    }
    // Loop over all frames in cluster
    int set = 0;
    Frame clusterframe = coords_->AllocateFrame();
    for (Node::frame_iterator fnum = C->beginframe();
                                     fnum != C->endframe(); ++fnum)
    {
      coords_->GetFrame( *fnum, clusterframe );
      clusterout.WriteSingle(set++, clusterframe);
    }
    // Close traj
    clusterout.EndTraj();
  }
}

/** Write average of clusters to files. */
void Cpptraj::Cluster::Results_Coords::WriteAvgStruct( List const& CList ) const {
  Topology avgparm = coords_->Top();
  // Get extension for representative frame format 
  std::string tmpExt = TrajectoryFile::WriteFormatExtension(avgfmt_);
  // Loop over all clusters
  for (List::cluster_iterator C = CList.begincluster();
                                     C != CList.endcluster(); ++C)
  {
    // Create filename based on cluster number.
    int cnum = C->Num();
    std::string cfilename = avgfile_ + ".c" + integerToString( cnum ) + tmpExt;
    // Set up trajectory file
    Trajout_Single clusterout; // FIXME CoordinateInfo OK for just coords?
    if (clusterout.PrepareTrajWrite(cfilename, ArgList(), &avgparm,
                                    CoordinateInfo(), 1, avgfmt_))
    {
      mprinterr("Error: Could not set up cluster average file %s for write.\n",
                cfilename.c_str());
      return;
    }
    // Get rep frame for rms fitting.
    Frame repframe = coords_->AllocateFrame();
    coords_->GetFrame( C->BestRepFrame(), repframe );
    Vec3 reftrans = repframe.CenterOnOrigin(false);
    // Loop over all frames in cluster
    Frame clusterframe = coords_->AllocateFrame();
    Frame avgframe = clusterframe;
    avgframe.ZeroCoords();
    for (Node::frame_iterator fnum = C->beginframe();
                                     fnum != C->endframe(); ++fnum)
    {
      coords_->GetFrame( *fnum, clusterframe );
      clusterframe.RMSD_FitToRef( repframe, reftrans );
      avgframe += clusterframe;
    }
    avgframe.Divide( (double)C->Nframes() );
    clusterout.WriteSingle(0, avgframe);
    clusterout.EndTraj();
  }
}
 
/** Write representative frame of each cluster to a trajectory file.  */
void Cpptraj::Cluster::Results_Coords::WriteSingleRepTraj( List const& CList ) const {
  Trajout_Single clusterout;
  // Set up trajectory file. Use parm from COORDS DataSet. 
  Topology *clusterparm = coords_->TopPtr();
  int nRepsToSave = CList.front().BestReps().size();
  if (clusterout.PrepareTrajWrite(singlerepfile_, ArgList(), clusterparm,
                                  coords_->CoordsInfo(), CList.Nclusters() * nRepsToSave,
                                  singlerepfmt_)) 
  {
    mprinterr("Error: Could not set up single trajectory for represenatatives %s for write.\n",
                singlerepfile_.c_str());
     return;
  }
  // Set up frame to hold cluster rep coords. 
  Frame clusterframe = coords_->AllocateFrame();
  int framecounter = 0;
  // Write rep frames from all clusters.
  for (List::cluster_iterator cluster = CList.begincluster(); 
                                     cluster != CList.endcluster(); ++cluster) 
  {
    for (Node::RepPairArray::const_iterator rep = cluster->BestReps().begin();
                                                   rep != cluster->BestReps().end(); ++rep)
    {
      coords_->GetFrame( rep->first, clusterframe );
      clusterout.WriteSingle(framecounter++, clusterframe);
    }
  }
  // Close traj
  clusterout.EndTraj();
}

/** Write representative frame of each cluster to a separate trajectory file,
  * repfile.REPNUM.FMT
  */
void Cpptraj::Cluster::Results_Coords::WriteRepTraj( List const& CList ) const {
  // Get extension for representative frame format 
  std::string tmpExt = TrajectoryFile::WriteFormatExtension(reptrajfmt_);
  // Use Topology from COORDS DataSet to set up input frame
  Topology* clusterparm = coords_->TopPtr();
  Frame clusterframe = coords_->AllocateFrame();
  // Loop over all clusters
  for (List::cluster_iterator C = CList.begincluster();
                                     C != CList.endcluster(); ++C)
  {
    if (writeRepFrameNum_) {
      // Each rep from cluster to separate file.
      for (Node::RepPairArray::const_iterator rep = C->BestReps().begin();
                                                     rep != C->BestReps().end(); ++rep)
      {
        Trajout_Single clusterout;
        // Get best rep frame # 
        int framenum = rep->first;
        // Create filename based on cluster number and frame #
        std::string cfilename = reptrajfile_ + ".c" + integerToString(C->Num()) +
                                ("." + integerToString(framenum+1)) + tmpExt;
        // Set up trajectory file.
        if (clusterout.PrepareTrajWrite(cfilename, ArgList(), clusterparm,
                                        coords_->CoordsInfo(), 1, reptrajfmt_))
        {
          mprinterr("Error: Could not set up representative trajectory file %s for write.\n",
                    cfilename.c_str());
          return;
        }
        // Write cluster rep frame
        coords_->GetFrame( framenum, clusterframe );
        clusterout.WriteSingle(framenum, clusterframe);
        // Close traj
        clusterout.EndTraj();
      }
    } else {
      // Each rep from cluster to single file.
      Trajout_Single clusterout;
      // Create filename based on cluster number
      std::string cfilename = reptrajfile_ + ".c" + integerToString(C->Num()) + tmpExt;
      // Set up trajectory file.
      int nRepsToSave = C->BestReps().size();
      if (clusterout.PrepareTrajWrite(cfilename, ArgList(), clusterparm,
                                      coords_->CoordsInfo(), nRepsToSave, reptrajfmt_))
      {
        mprinterr("Error: Could not set up representative trajectory file %s for write.\n",
                  cfilename.c_str());
        return;
      }
      int framecounter = 0;
      for (Node::RepPairArray::const_iterator rep = C->BestReps().begin();
                                                     rep != C->BestReps().end(); ++rep)
      {
        // Write cluster rep frame
        coords_->GetFrame( rep->first, clusterframe );
        clusterout.WriteSingle( framecounter++, clusterframe );
      }
      // Close traj
      clusterout.EndTraj();
    }
  }
}
