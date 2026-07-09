#include "EnsembleOutList.h"
#include "Topology.h"
#include "CpptrajStdio.h"
#include "EnsembleOut_Single.h"
#include "EnsembleOut_Multi.h"
#include "ArgList.h"

/// CONSTRUCTOR
EnsembleOutList::EnsembleOutList() : debug_(0) {}

/// DESTRUCTOR
EnsembleOutList::~EnsembleOutList() { Clear(); }

/** Set the list debug level. Will apply to new output ensembles. */
void EnsembleOutList::SetDebug(int d) { debug_ = d; }

/** Clear the output ensemble list and free memory. */
void EnsembleOutList::Clear() {
  for (EnsArray::const_iterator ens = ensout_.begin(); ens != ensout_.end(); ++ens)
    delete *ens;
  ensout_.clear();
  ensTops_.clear();
  active_.clear();
  open_.clear();
}

// TODO Pass in more ensemble information, maps etc?
int EnsembleOutList::AddEnsembleOut(std::string const& fname, ArgList const& args,
                                    DataSetList const& DSLin,
                                    Topology* eParm, int ensembleSize)
{
  if (eParm == 0) {
    mprinterr("Error: No topology information.\n");
    return 1;
  }
  if (fname.empty()) {
    mprinterr("Internal Error: EnsembleOutList::AddEnsembleOut() called with empty filename.\n");
    return 1;
  }
  // Determine if this filename is in use in order to prevent overwrites
  for (EnsArray::const_iterator eo = ensout_.begin();
                                eo != ensout_.end(); ++eo)
  {
    if ( (*eo)->Traj().Filename().Full() == fname ) {
      mprinterr("Error: Output ensemble filename %s already in use.\n", fname.c_str());
      return 1;
    }
  }
  ArgList argIn = args;
  EnsembleOut* ens = 0;
# ifdef ENABLE_SINGLE_ENSEMBLE
  // See if single ensemble output desired. // FIXME: Should not depend on keyword
  if (argIn.hasKey("ensemble"))
    ens = new EnsembleOut_Single();
  else
# endif
    // Create new multi output trajectory
    ens = new EnsembleOut_Multi();
  if (ens == 0) return 1;
  if (ens->InitEnsembleWrite(fname, argIn, DSLin, ensembleSize, TrajectoryFile::UNKNOWN_TRAJ)) {
    delete ens;
    return 1;
  }
  ensout_.push_back( ens );
  ensTops_.push_back( eParm );
  open_.push_back(false);
  return 0;
}

// EnsembleOutList::SetupEnsembleOut()
int EnsembleOutList::SetupEnsembleOut(Topology* CurrentParm, CoordinateInfo const& cInfo,
                                      int Nframes)
{
  active_.clear();
  for (unsigned int i = 0; i != ensout_.size(); i++) {
    // Check that input parm matches setup parm - if not, skip
    if (CurrentParm->Pindex() == ensTops_[i]->Pindex()) {
      if (!open_[i]) {
        if ( ensout_[i]->SetupEnsembleWrite( CurrentParm, cInfo, Nframes ) )
        {
          mprinterr("Error: Setting up output ensemble %s\n", ensout_[i]->Traj().Filename().full());
          return 1;
        }
        open_[i] = true;
      }
      active_.push_back( ensout_[i] );
    }
  }
  ListActive();
  return 0;
}

/** Check active output ensembles against the current coordinate info. */
int EnsembleOutList::CheckEnsembleOutCoordInfo(CoordinateInfo const& currentCoordInfo,
                                               FramePtrArray const& currentFrames)
const
{
  for (EnsArray::const_iterator it = active_.begin(); it != active_.end(); ++it)
  {
    EnsembleOut* activeTraj = *it;
    CoordinateInfo const& ensembleoutCoordInfo = activeTraj->Traj().CoordInfo();
    // Velocities
    if (currentCoordInfo.HasVel() != ensembleoutCoordInfo.HasVel()) {
      if (currentCoordInfo.HasVel()) {
        mprintf("Warning: Input ensemble has velocity information but output ensemble was set up without velocities.\n");
        mprintf("Warning: Velocity information will not be written to output ensemble.\n");
      } else {
        mprintf("Warning: Input ensemble has no velocity information but output ensemble was set up with velocities.\n");
        mprintf("Warning: All zeroes will be written for velocities.\n");
        for (FramePtrArray::const_iterator frm = currentFrames.begin(); frm != currentFrames.end(); ++frm)
          (*frm)->AddVelocities(Frame::Darray((*frm)->size(), 0.0));
      }
    }
    // Forces
    if (currentCoordInfo.HasForce() != ensembleoutCoordInfo.HasForce()) {
      if (currentCoordInfo.HasForce()) {
        mprintf("Warning: Input ensemble has force information but output ensemble was set up without forces.\n");
        mprintf("Warning: Force information will not be written to output ensemble.\n");
      } else {
        mprintf("Warning: Input ensemble has no force information but output ensemble was set up with forces.\n");
        mprintf("Warning: All zeroes will be written for forces.\n");
        for (FramePtrArray::const_iterator frm = currentFrames.begin(); frm != currentFrames.end(); ++frm)
          (*frm)->AddForces(Frame::Darray((*frm)->size(), 0.0));
      }
    }
  }
  return 0;
}

/** List only active output trajectories. */
void EnsembleOutList::ListActive() const {
  if (!ensout_.empty()) {
    mprintf(".....................................................\n");
    if (!active_.empty()) {
      mprintf("ACTIVE OUTPUT ENSEMBLES (%zu):\n", active_.size());
      for (EnsArray::const_iterator it = active_.begin(); it != active_.end(); ++it)
      {
        mprintf("  %s", (*it)->Traj().Filename().full());
        std::string meta = (*it)->Traj().CoordInfo().InfoString();
        if (!meta.empty()) mprintf(" (%s)", meta.c_str());
        mprintf("\n");
      }
    } else
      mprintf("NO ACTIVE OUTPUT ENSEMBLES.\n");
  }
}

/** Go through each active output traj, call write. */
int EnsembleOutList::WriteEnsembleOut(int set, FramePtrArray const& Farray)
{
  for (EnsArray::const_iterator ens = active_.begin(); ens != active_.end(); ++ens) {
    if ( (*ens)->WriteEnsemble(set, Farray) ) {
      mprinterr("Error writing output ensemble, frame %i.\n", set+1);
      return 1;
    }
  }
  return 0;
}

/** Close output trajectories. Called after input traj processing completed. */
void EnsembleOutList::CloseEnsembleOut() {
  for (EnsArray::const_iterator ens = ensout_.begin(); ens != ensout_.end(); ++ens)
    (*ens)->EndEnsemble();
  Clear();
}

void EnsembleOutList::List(std::vector<int> const& PindexFrames) const {
  if (!ensout_.empty()) {
    //mprintf("\nOUTPUT ENSEMBLE:\n");
    mprintf("\nENSEMBLE OUTPUT TRAJECTORIES (Numerical filename"
            " suffix corresponds to above map):\n");
    if (PindexFrames.empty())
      for (unsigned int i = 0; i != ensout_.size(); i++)
        ensout_[i]->PrintInfo( 0 );
    else
      for (unsigned int i = 0; i != ensout_.size(); i++)
        ensout_[i]->PrintInfo( PindexFrames[ensTops_[i]->Pindex()] );
  }
}
#ifdef MPI
// -----------------------------------------------------------------------------
int EnsembleOutList::ParallelSetupEnsembleOut(Topology* CurrentParm,
                                              CoordinateInfo const& cInfo, int Nframes,
                                              Parallel::Comm const& commIn)
{
  active_.clear();
  for (unsigned int i = 0; i != ensout_.size(); i++) {
    // Check that input parm matches setup parm - if not, skip
    if (CurrentParm->Pindex() == ensTops_[i]->Pindex()) {
      if (!open_[i]) {
        ensout_[i]->SetTrajComm( commIn );
        if ( ensout_[i]->SetupEnsembleWrite( CurrentParm, cInfo, Nframes ) )
        {
          mprinterr("Error: Setting up output ensemble '%s' in parallel\n",
                    ensout_[i]->Traj().Filename().full());
          return 1;
        }
        open_[i] = true;
      }
      active_.push_back( ensout_[i] );
    }
  }
  ListActive();
  return 0;
}
#endif
