#include "TrajoutList.h"
#include "CpptrajStdio.h"
#include "EnsembleOut_Single.h"
#include "EnsembleOut_Multi.h"

void EnsembleOutList::Clear() {
  for (EnsArray::const_iterator ens = ensout_.begin(); ens != ensout_.end(); ++ens)
    delete *ens;
  ensout_.clear();
  ensTops_.clear();
  active_.clear();
  open_.clear();
}

int EnsembleOutList::AddEnsembleOut(std::string const& fname, ArgList const& args,
                                    Topology* eParm, int ensembleSize,
                                    TrajectoryFile::TrajFormatType fmt)
{
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
  if (ens->InitEnsembleWrite(fname, argIn, ensembleSize, fmt)) {
    delete ens;
    return 1;
  }
  ensout_.push_back( ens );
  ensTops_.push_back( eParm );
  open_.push_back(false);
  return 0;
}

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
  return 0;
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
    mprintf("\nOUTPUT ENSEMBLE:\n");
    if (PindexFrames.empty())
      for (unsigned int i = 0; i != ensout_.size(); i++)
        ensout_[i]->PrintInfo( 0 );
    else
      for (unsigned int i = 0; i != ensout_.size(); i++)
        ensout_[i]->PrintInfo( PindexFrames[ensTops_[i]->Pindex()] );
  }
}

// =============================================================================
void TrajoutList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_ > 0)
    mprintf("TrajoutList debug level set to %i\n", debug_);
}

void TrajoutList::Clear() {
  for (ListType::iterator traj = trajout_.begin(); traj != trajout_.end(); ++traj)
    delete *traj;
  trajout_.clear();
  trajoutArgs_.clear();
  trajoutTops_.clear();
  trajoutNames_.clear();
  active_.clear();
  open_.clear();
}

/** Add output trajectory to list as single output trajectory. Associate it
  * with the given Topology but no Topology-dependent setup will occur. This
  * is because during the course of a Run the Topology may be modified, by
  * e.g. a 'strip' command.
  */
int TrajoutList::AddTrajout(std::string const& filename, ArgList const& argIn, Topology* tParm)
{
  if (tParm == 0) {
    mprinterr("Error: No topology information.\n");
    return 1;
  }
  if (filename.empty()) {
    mprinterr("Internal Error: TrajoutList::AddTrajout() called with empty filename.\n");
    return 1;
  }
  // Determine if this filename is in use in order to prevent overwrites
  for (ListType::const_iterator to = trajout_.begin();
                                to != trajout_.end(); ++to)
  {
    if ( (*to)->Traj().Filename().Full() == filename ) {
      mprinterr("Error: Output trajectory filename %s already in use.\n",filename.c_str());
      return 1;
    }
  }
  // Create Trajout_Single
  Trajout_Single* to = new Trajout_Single();
  to->SetDebug( debug_ );
  // Initialize output trajectory
  ArgList args = argIn;
  if (to->InitTrajWrite(filename, args, TrajectoryFile::UNKNOWN_TRAJ)) {
    mprinterr("Error: Could not set up output trajectory.\n");
    delete to;
    return 1;
  }
  trajout_.push_back( to );
  // For potentially setting up ensemble later, save trajout arg.
  trajoutArgs_.push_back( argIn );
  trajoutTops_.push_back( tParm );
  trajoutNames_.push_back( filename );
  open_.push_back( false );
  return 0;
}

// TODO Pass in more ensemble information, maps etc?
int TrajoutList::MakeEnsembleTrajout(EnsembleOutList& ensembleList,
                                     int ensembleSize) const
{
  ensembleList.Clear();
  for (unsigned int i = 0; i != trajoutArgs_.size(); i++) {
    if (ensembleList.AddEnsembleOut(trajoutNames_[i], trajoutArgs_[i], trajoutTops_[i],
                                    ensembleSize, trajout_[i]->Traj().WriteFormat()))
      return 1;
  }
  return 0;
}

// TrajoutList::SetupTrajout()
int TrajoutList::SetupTrajout(Topology* CurrentParm, CoordinateInfo const& cInfo,
                              int Nframes)
{
  active_.clear();
  for (unsigned int i = 0; i != trajout_.size(); i++) {
    // Check that input parm matches setup parm - if not, skip
    if (CurrentParm->Pindex() == trajoutTops_[i]->Pindex()) {
      if (!open_[i]) { // Only set up if not already open.
        if ( trajout_[i]->SetupTrajWrite( CurrentParm, cInfo, Nframes) )
        {
          mprinterr("Error: Setting up output trajectory %s\n", trajoutNames_[i].c_str());
          return 1;
        }
        open_[i] = true;
      }
      active_.push_back( trajout_[i] );
    }
  }
  return 0;
}

// TrajoutList::WriteTrajout()
/** Go through each active output traj, call write. */
int TrajoutList::WriteTrajout(int set, Frame const& CurrentFrame)
{
  for (ListType::const_iterator traj = active_.begin();
                                traj != active_.end(); ++traj)
  {
    if ( (*traj)->WriteSingle(set, CurrentFrame) ) {
      mprinterr("Error writing output trajectory, frame %i.\n", set+1);
      return 1;
    }
  }
  return 0;
}

// TrajoutList::CloseTrajout()
/** Close output trajectories. Called after input traj processing completed. */
void TrajoutList::CloseTrajout() {
  for (ListType::const_iterator traj = trajout_.begin();
                                traj != trajout_.end(); ++traj)
    (*traj)->EndTraj();
  Clear();
}

// TrajoutList::List()
void TrajoutList::List(std::vector<int> const& PindexFrames) const {
  if (!trajout_.empty()) {
    mprintf("\nOUTPUT TRAJECTORIES (%zu total):\n", trajout_.size());
    if (PindexFrames.empty())
      for (unsigned int i = 0; i != trajout_.size(); i++)
        trajout_[i]->PrintInfo( 0 );
    else
      for (unsigned int i = 0; i != trajout_.size(); i++)
        trajout_[i]->PrintInfo( PindexFrames[trajoutTops_[i]->Pindex()] );
  }
}
#ifdef MPI
// -----------------------------------------------------------------------------
int TrajoutList::ParallelSetupTrajout(Topology* CurrentParm,
                                      CoordinateInfo const& cInfo, int Nframes,
                                      Parallel::Comm const& commIn)
{
  active_.clear();
  for (unsigned int i = 0; i != trajout_.size(); i++) {
    // Check that input parm matches setup parm - if not, skip
    if (CurrentParm->Pindex() == trajoutTops_[i]->Pindex()) {
      if (!open_[i]) { // Only set up if not already open.
        if ( trajout_[i]->ParallelSetupTrajWrite( CurrentParm, cInfo, Nframes, commIn) )
        {
          mprinterr("Error: Setting up output trajectory '%s' in parallel.\n",
                    trajoutNames_[i].c_str());
          return 1;
        }
        open_[i] = true;
      }
      active_.push_back( trajout_[i] );
    } else {
      mprintf("Warning: Output traj '%s' was set up for topology '%s', but\n"
              "Warning:   parallel run topology is '%s' - skipping.\n",
              trajout_[i]->Traj().Filename().full(), trajoutTops_[i]->c_str(),
              CurrentParm->c_str());
    }
  }
  return 0;
}

int TrajoutList::ParallelWriteTrajout(int set, Frame const& CurrentFrame)
{
  for (ListType::const_iterator traj = active_.begin();
                                traj != active_.end(); ++traj)
  {
    if ( (*traj)->ParallelWriteSingle(set, CurrentFrame) ) {
      mprinterr("Error writing output trajectory in parallel, frame %i.\n", set+1);
      return 1;
    }
  }
  return 0;
}

void TrajoutList::ParallelCloseTrajout() {
  for (ListType::const_iterator traj = trajout_.begin();
                                traj != trajout_.end(); ++traj)
    (*traj)->ParallelEndTraj();
  Clear();
}
#endif
