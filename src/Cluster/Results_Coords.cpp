#include "Results_Coords.h"
#include "List.h"
#include "Metric_DME.h"
#include "Metric_RMS.h"
#include "Metric_SRMSD.h"
#include "MetricArray.h"
#include "Node.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../DataSet_Coords.h"
#include "../DataSetList.h" // For PrepareTrajWrite
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

/** Help for Results_Coords keywords. */
void Cpptraj::Cluster::Results_Coords::Help() {
  mprintf("\t[clusterout <trajfileprefix> [clusterfmt <trajformat>]]\n"
          "\t[singlerepout <trajfilename> [singlerepfmt <trajformat>]]\n"
          "\t[repout <repprefix> [repfmt <trajformat>] [repframe]]\n"
          "\t[avgout <avgprefix> [avgfmt <trajformat>]]\n"
          "\t[assignrefs [refcut <rms>] [refmask <mask>]]\n");
}

/** Get arguments related to writing cluster data to trajectories. */
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

/** Get user specified options. */
int Cpptraj::Cluster::Results_Coords::GetOptions(ArgList& analyzeArgs, DataSetList const& DSL,
                                                 MetricArray const& metricIn)
{
  writeRepFrameNum_ = analyzeArgs.hasKey("repframe");
  GetClusterTrajArgs(analyzeArgs, "clusterout",   "clusterfmt",   clusterfile_,   clusterfmt_);
  GetClusterTrajArgs(analyzeArgs, "singlerepout", "singlerepfmt", singlerepfile_, singlerepfmt_);
  GetClusterTrajArgs(analyzeArgs, "repout",       "repfmt",       reptrajfile_,   reptrajfmt_);
  GetClusterTrajArgs(analyzeArgs, "avgout",       "avgfmt",       avgfile_,       avgfmt_);

  if (analyzeArgs.hasKey("assignrefs")) {
    refSets_ = DSL.GetSetsOfType("*", DataSet::REF_FRAME);
    if (refSets_.empty()) {
      mprinterr("Error: 'assignrefs' specified but no references loaded.\n");
      return 1;
    }
    refCut_ = analyzeArgs.getKeyDouble("refcut", 1.0);
    refmaskexpr_ = analyzeArgs.GetStringKey("refmask");
    useMass_ = false;
    // Attempt to set defaults from Metric if applicable
    Metric const* coordsMetric = metricIn.CoordsMetric();
    if (coordsMetric != 0) {
      if (coordsMetric->MetricType() == Metric::RMS) {
        Metric_RMS const& met = static_cast<Metric_RMS const&>( *coordsMetric );
        useMass_ = met.UseMass();
        if (refmaskexpr_.empty()) refmaskexpr_ = met.Mask().MaskString();
      } else if (coordsMetric->MetricType() == Metric::SRMSD) {
        Metric_SRMSD const& met = static_cast<Metric_SRMSD const&>( *coordsMetric );
        useMass_ = met.UseMass();
        if (refmaskexpr_.empty()) refmaskexpr_ = met.Mask().MaskString();
      } else if (coordsMetric->MetricType() == Metric::DME) {
        Metric_DME const& met = static_cast<Metric_DME const&>( *coordsMetric );
        if (refmaskexpr_.empty()) refmaskexpr_ = met.Mask().MaskString();
      }
    }
    // Set a default mask if needed
    if (refmaskexpr_.empty()) {
      refmaskexpr_.assign("!@H=");
      mprintf("Warning: 'assignrefs' specified but no 'refmask' given.\n"
              "Warning:   Using default mask expression: '%s'\n", refmaskexpr_.c_str());
    }
  }

  return 0;
}

/** Write info on what results will be calculated/written. */
void Cpptraj::Cluster::Results_Coords::Info() const {
  if (coords_ == 0) {
    mprintf("\tNo coordinates set provided for cluster results.\n");
    return;
  }

  mprintf("\tCoordinates set for cluster results: %s\n", coords_->legend());
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

  if (!refSets_.empty())
    mprintf("\tClusters will be identified with loaded reference structures if RMSD\n"
            "\t  (mask '%s') to representative frame is < %g Ang.\n",
            refmaskexpr_.c_str(), refCut_);
}

/** Do any output needed. */
int Cpptraj::Cluster::Results_Coords::DoOutput(List const& CList) const {
  int err = 0;
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
  return err;
}

/** Calculate any results that belong in clusters. */
int Cpptraj::Cluster::Results_Coords::CalcResults(List& CList) const {
  int err = 0;
  // Assign reference structures to clusters
  if (!refSets_.empty())
    err += AssignRefsToClusters( CList );
  return err;
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
    if (clusterout.PrepareTrajWrite(cfilename, ArgList(), DataSetList(), clusterparm,
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
    if (clusterout.PrepareTrajWrite(cfilename, ArgList(), DataSetList(), &avgparm,
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
  if (clusterout.PrepareTrajWrite(singlerepfile_, ArgList(), DataSetList(), clusterparm,
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
        if (clusterout.PrepareTrajWrite(cfilename, ArgList(), DataSetList(), clusterparm,
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
      if (clusterout.PrepareTrajWrite(cfilename, ArgList(), DataSetList(), clusterparm,
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

/** Assign reference frames to clusters via closest RMSD. */
int Cpptraj::Cluster::Results_Coords::AssignRefsToClusters( List& CList )
const
{
  // Pre-center all reference coords at the origin. No need to store trans vectors.
  std::vector<Frame> refFrames;
  refFrames.reserve( refSets_.size() );
  for (unsigned int idx = 0; idx != refSets_.size(); idx++) {
    AtomMask rMask( refmaskexpr_ );
    DataSet_Coords_REF* REF_ds = (DataSet_Coords_REF*)refSets_[idx];
    if ( REF_ds->Top().SetupIntegerMask( rMask, REF_ds->RefFrame() ) ) {
      mprintf("Warning: Could not set up mask for reference '%s'\n", REF_ds->legend());
      continue;
    }
    refFrames.push_back( Frame(REF_ds->RefFrame(), rMask) );
    refFrames.back().CenterOnOrigin( useMass_ );
  }
  // For each cluster, assign the reference name with the lowest RMSD
  // to the representative frame that is below the cutoff.
  AtomMask tMask( refmaskexpr_ );
  if (coords_->Top().SetupIntegerMask( tMask )) {
    mprinterr("Error: Could not set up mask for assigning references.\n");
    return 1;
  }
  Frame TGT( coords_->AllocateFrame(), tMask );
  unsigned int cidx = 0;
  for (List::cluster_it cluster = CList.begin();
                        cluster != CList.end(); ++cluster, ++cidx)
  {
    coords_->GetFrame( cluster->BestRepFrame(), TGT, tMask );
    double minRms = TGT.RMSD_CenteredRef( refFrames[0], useMass_ );
    unsigned int minIdx = 0;
    for (unsigned int idx = 1; idx < refSets_.size(); idx++) {
      double rms = TGT.RMSD_CenteredRef( refFrames[idx], useMass_ );
      if (rms < minRms) {
        minRms = rms;
        minIdx = idx;
      }
    }
    if (minRms < refCut_) {
      //mprintf("DEBUG: Assigned cluster %i to reference \"%s\" (%g)\n", cidx,
      //        refSets_[minIdx]->Meta().Name().c_str(), minRms);
      cluster->SetNameAndRms( refSets_[minIdx]->Meta().Name(), minRms );
    } else {
      //mprintf("DEBUG: Cluster %i was closest to reference \"(%s)\" (%g)\n", cidx,
      //        refSets_[minIdx]->Meta().Name().c_str(), minRms);
      cluster->SetNameAndRms( "(" + refSets_[minIdx]->Meta().Name() + ")", minRms );
    }
  }
  return 0;
}
