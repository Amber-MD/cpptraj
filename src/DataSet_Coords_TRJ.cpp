#include "DataSet_Coords_TRJ.h"
#include "Trajin_Single.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataSet_Coords_TRJ::DataSet_Coords_TRJ() :
  DataSet_Coords(TRAJ),
  Traj_(0), 
  deleteTrajectories_(false)
{}

// DESTRUCTOR
DataSet_Coords_TRJ::~DataSet_Coords_TRJ() {
  if (deleteTrajectories_) {
    for (ListType::const_iterator trj = trajinList_.begin(); trj != trajinList_.end(); ++trj)
      delete *trj;
  }
}

int DataSet_Coords_TRJ::CoordsSetup(Topology const& topIn, CoordinateInfo const& cInfoIn ) {
  if (trajinList_.empty()) {
    top_ = topIn;
    cInfo_ = cInfoIn;
  } else {
    if ( topIn.Natom() != top_.Natom() ) {
      mprinterr("Error: For TRAJ data set currently all trajectories must have same number\n"
                "Error:  of atoms: %i != %i. Recommended course of action is to\n"
                "Error:  create a trajectory where all frames have been stripped to the same\n"
                "Error:  number of atoms first.\n", topIn.Natom(), top_.Natom());
      return 1;
    }
    // Since velocity info is not always allocated in Frame, if one traj
    // has velocity info ensure that all do.
    if (cInfoIn.HasVel()) cInfo_.SetVelocity( true );
  }
  return 0;
}

// DataSet_Coords_TRJ::UpdateTrjFrames()
int DataSet_Coords_TRJ::UpdateTrjFrames(TrajFrameCounter const& count) {
  if (count.TotalReadFrames() > 0)
    IDX_.AddTraj( count.TotalReadFrames(), count.Start(), count.Offset() );
  else {
    mprinterr("Error: Cannot use trajectories with unknown # of frames as data set.\n");
    return 1;
  }
  return 0;
}

/** Add a new Trajin_Single class. */
int DataSet_Coords_TRJ::AddSingleTrajin(std::string const& fname, ArgList& argIn,
                                        Topology* parmIn) // TODO: const&
{
  if (parmIn==0) return 1;
  if (!trajinList_.empty() && !deleteTrajectories_) {
    mprinterr("Internal Error: This DataSet_Coords_TRJ class set up for copies.\n");
    return 1;
  }
  Trajin_Single* trajin = new Trajin_Single();
  if (trajin->SetupTrajRead(fname, argIn, parmIn)) {
    mprinterr("Error: Could not set up trajectory '%s'\n", fname.c_str());
    return 1;
  }
  if (CoordsSetup(*parmIn, trajin->TrajCoordInfo())) return 1;
  if (UpdateTrjFrames( trajin->Traj().Counter() )) return 1;
  trajinList_.push_back( trajin );
  deleteTrajectories_ = true;
  return 0;
} 

/** Add an existing and already set-up Trajin class; this class will not be
  * deleted upon destruction. Intended for use with TrajinList.
  */
int DataSet_Coords_TRJ::AddInputTraj(Trajin* tIn) {
  if (!trajinList_.empty() && deleteTrajectories_) {
    mprinterr("Internal Error: This DataSet_Coords_TRJ class not set up for copies.\n");
    return 1;
  }
  if (tIn == 0) return 1;
  if (CoordsSetup( *(tIn->Traj().Parm()), tIn->TrajCoordInfo() )) return 1;
  if (UpdateTrjFrames( tIn->Traj().Counter() )) return 1;
  trajinList_.push_back( tIn );
  deleteTrajectories_ = false;
  return 0;
}

void DataSet_Coords_TRJ::GetFrame(int idx, Frame& fIn) {
# ifdef _OPENMP
# pragma omp critical
  {
# endif
  // Determine which trajectory has the desired index
  int internalIdx = IDX_.FindIndex( idx );
  if (internalIdx < 0) {
    mprinterr("Internal Error: Global index %i is out of range.\n", idx);
    return;
  }
  // If desired traj is different than current, open desired traj
  if (IDX_.TrajHasChanged()) {
    if (Traj_ != 0) Traj_->EndTraj();
    Trajin* previousTraj = Traj_;
    Traj_ = trajinList_[ IDX_.CurrentTrajNum() ];
    // NOTE: Currently enforcing all traj have same # atoms, no need to check topology.
    // See if frame needs (re-)allocation.
    if (previousTraj == 0 || previousTraj->TrajCoordInfo() != Traj_->TrajCoordInfo())
      readFrame_.SetupFrameV( Traj_->Traj().Parm()->Atoms(), Traj_->TrajCoordInfo() );
    // Open traj.
    if (Traj_->BeginTraj()) {
      mprinterr("Error: Could not open trajectory %i '%s'\n", IDX_.CurrentTrajNum(),
                Traj_->Traj().Filename().full());
      return;
    }
  }
  // Read the frame
  // TODO: May need to use readFrame here as well...
  if (Traj_->ReadTrajFrame( internalIdx, fIn )) {
    mprinterr("Error: Could not read '%s' frame %i\n",
                Traj_->Traj().Filename().full(), internalIdx + 1);
  }
# ifdef _OPENMP
  }
# endif
}

void DataSet_Coords_TRJ::GetFrame(int idx, Frame& fIn, AtomMask const& mask) {
  GetFrame( idx, readFrame_ );
  fIn.SetFrame( readFrame_, mask );
}

void DataSet_Coords_TRJ::Info() const {
  if (trajinList_.size() == 1)
    mprintf(" (1 trajectory)");
  else
    mprintf(" (%zu trajectories)", trajinList_.size());
  CommonInfo();
}
