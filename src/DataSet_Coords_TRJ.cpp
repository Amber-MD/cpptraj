#include "DataSet_Coords_TRJ.h"
#include "CpptrajStdio.h"

DataSet_Coords_TRJ::DataSet_Coords_TRJ() :
  DataSet_Coords(TRAJ),
  Traj_(0), 
  currentTrajNum_(-1),
  globalOffset_(0),
  maxFrames_(0)
{}

int DataSet_Coords_TRJ::AddInputTraj(Trajin* tIn) {
  if (tIn == 0) return 1;
  if (trajinList_.empty())
    SetTopology( *(tIn->TrajParm()) );
  else {
    if ( tIn->TrajParm()->Natom() != Topology().Natom() ) {
      mprinterr("Error: For TRAJ data set currently all trajectories must have same number\n"
                "Error:  of atoms (%i != %i). Recommended course of action is to create a\n"
                "Error:  trajectory where all frames have been stripped to the same number of"
                "Error:  atoms first.\n");
      return 1;
    }
  }
  // TODO: Need some way of enforcing topology sameness
  if (tIn->TotalReadFrames() > 0)
    maxFrames_ += tIn->TotalReadFrames();
  else {
    mprinterr("Error: Cannot use trajectories with unknown # of frames as data set.\n");
    return 1;
  }
  trajinList_.push_back( tIn );
  return 0;
}

void DataSet_Coords_TRJ::GetFrame(int idx, Frame& fIn) {
  // Determine which trajectory has the desired index
  globalOffset_ = 0;
  int currentMax = 0;
  int desiredTrajNum = 0;
  for (; desiredTrajNum < (int)trajinList_.size(); ++desiredTrajNum) {
    currentMax += trajinList_[desiredTrajNum]->TotalReadFrames();
    if (idx < currentMax) break;
    globalOffset_ += trajinList_[desiredTrajNum]->TotalReadFrames();
  }
  // If desired traj is different than current, open desired traj
  if ( desiredTrajNum != currentTrajNum_ ) {
    if (Traj_ != 0) Traj_->EndTraj();
    Trajin* prevTraj = Traj_;
    currentTrajNum_ = desiredTrajNum;
    Traj_ = trajinList_[currentTrajNum_];
    // Set up frame for reading
    bool parmHasChanged;
    if (prevTraj == 0)
      parmHasChanged = true;
    else
      parmHasChanged = ( prevTraj->TrajParm()->Natom() != 
                            Traj_->TrajParm()->Natom()    );
    if ( parmHasChanged || (readFrame_.HasVelocity() != 
                            Traj_->HasVelocity()) )
      readFrame_.SetupFrameV(Traj_->TrajParm()->Atoms(), Traj_->HasVelocity(),
                             Traj_->NreplicaDimension());
    if (Traj_->BeginTraj(false)) {
      mprinterr("Error: Could not open trajectory %i '%s'\n", desiredTrajNum,
                Traj_->TrajFilename().full());
      return;
    }
  }
  // Convert desired index into trajectory internal index
  int internalIdx = ((idx - globalOffset_) * Traj_->Offset()) + Traj_->Start();
  // Read desired index
  // TODO: May need to use readFrame here as well...
  if (Traj_->ReadTrajFrame( internalIdx, fIn )) {
    mprinterr("Error: Could not read '%s' frame %i\n", 
              Traj_->TrajFilename().full(), internalIdx + 1);
  }
}

void DataSet_Coords_TRJ::GetFrame(int idx, Frame& fIn, AtomMask const& mask) {
  GetFrame( idx, readFrame_ );
  fIn.SetFrame( readFrame_, mask );
}

void DataSet_Coords_TRJ::Info() const {
  if (trajinList_.size() == 1)
    mprintf(" (%zu trajectories)", trajinList_.size());
  else
    mprintf(" (1 trajectory)");
}
