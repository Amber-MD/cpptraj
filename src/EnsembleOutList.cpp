#include "EnsembleOutList.h"
#include "CpptrajStdio.h"
#include "EnsembleOut_Single.h"
#include "EnsembleOut_Multi.h"

void EnsembleOutList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_ > 0)
    mprintf("EnsembleOutList debug level set to %i\n", debug_);
}

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
  if (ens->InitEnsembleWrite(fname, argIn, ensembleSize, TrajectoryFile::UNKNOWN_TRAJ)) {
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
