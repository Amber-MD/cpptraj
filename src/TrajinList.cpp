#include "TrajinList.h"
#include "CpptrajStdio.h"
#include "TrajectoryFile.h"
#include "Trajin_Single.h"
#include "Trajin_Multi.h"
#include "EnsembleIn_Single.h"
#include "EnsembleIn_Multi.h"
#include "StringRoutines.h" // ExpandToFilenames

TrajinList::TrajinList() : debug_(0), maxframes_(0), ensembleSize_(-1) {}

TrajinList::~TrajinList() { Clear(); }

void TrajinList::Clear() {
  for (tListType::iterator traj = trajin_.begin(); traj != trajin_.end(); ++traj)
    delete *traj;
  trajin_.clear();
  for (eListType::iterator ens = ensemble_.begin(); ens != ensemble_.end(); ++ens)
    delete *ens;
  ensemble_.clear();
  maxframes_ = 0;
  topFrames_.clear();
  ensembleSize_ = -1;
}

/** Update max # of frames in the list and # frames associated with Topologies
  * used by Trajectories in the list.
  */
void TrajinList::UpdateMaxFrames(InputTrajCommon const& traj) {
  int trajFrames = traj.Counter().TotalReadFrames();
  int pindex = traj.Parm()->Pindex();
  if (pindex >= (int)topFrames_.size())
    topFrames_.resize( pindex + 1 );
  // If < 0 frames this indicates the number of frames could not be determined. 
  if (trajFrames < 0) {
    maxframes_ = -1;
    topFrames_[pindex] = 0; // TODO should be -1?
  } else if (maxframes_ != -1) {
    // Only update # of frames if total # of frames so far is known.
    topFrames_[pindex] += trajFrames;
    maxframes_ += trajFrames;
  }
}

// TrajinList::AddEnsembleIn()
int TrajinList::AddEnsembleIn(std::string const& fname, Topology* topIn, ArgList const& argIn)
{
  if (topIn == 0) {
    mprinterr("Error: No topology for input ensemble '%s'\n", fname.c_str());
    return 1;
  }
  int err = 0;
  File::NameArray fnames = File::ExpandToFilenames( fname );
  if (fnames.empty()) return 1;
  TrajectoryFile::TrajFormatType trajinFmt;
  TrajectoryIO* tio = 0;
  for (File::NameArray::const_iterator fn = fnames.begin(); fn != fnames.end(); ++fn) {
    ArgList args = argIn;
    // Determine whether this file is multiple file or single file ensemble.
    tio = TrajectoryFile::DetectFormat( *fn, trajinFmt );
    if (tio == 0) {
      mprinterr("Error: Could not determine trajectory %s format\n", fn->full());
      err++;
      continue;
    }
    EnsembleIn* ensemble = 0;
#   ifdef ENABLE_SINGLE_ENSEMBLE
    if (tio->CanProcessEnsemble())
      ensemble = new EnsembleIn_Single();
    else
#   endif
      ensemble = new EnsembleIn_Multi();
    if (ensemble == 0) {
      mprinterr("Error: Memory allocation for input ensemble failed.\n");
      delete tio;
      return 1;
    }
    ensemble->SetDebug( debug_ );
    // CRDIDXARG: Append coordinate indices arg if there is one
    args.AddArg( finalCrdIndicesArg_ );
    if ( ensemble->SetupEnsembleRead(*fn, args, topIn) ) {
      mprinterr("Error: Could not set up input ensemble '%s'.\n", fname.c_str());
      delete ensemble;
      delete tio;
      err++;
      continue;
    }
    // Currently all input ensembles must be same size.
    if (ensembleSize_ == -1)
      ensembleSize_ = ensemble->EnsembleCoordInfo().EnsembleSize();
    else if (ensembleSize_ != ensemble->EnsembleCoordInfo().EnsembleSize()) {
      mprinterr("Error: Ensemble size (%i) does not match first ensemble size (%i).\n",
                ensemble->EnsembleCoordInfo().EnsembleSize(), ensembleSize_);
      return 1;
    }
    // CRDIDXARG: If trajectory is REMD ensemble and sorting by CRDIDX, need to
    //            save final CRDIDX for next ensemble command.
    // TODO: This is very clunky - remlog dataset should contain all exchanges
    //       so trajin doesnt have to worry about it.
#   ifdef ENABLE_SINGLE_ENSEMBLE
    if ( !tio->CanProcessEnsemble() ) {
#   endif
      EnsembleIn_Multi const& mTraj = static_cast<EnsembleIn_Multi const&>( *ensemble );
      if ( mTraj.TargetMode() == ReplicaInfo::CRDIDX ) {
        finalCrdIndicesArg_ = mTraj.FinalCrdIndices();
        if (finalCrdIndicesArg_.empty()) {
          mprinterr("Error: Could not obtain final remlog indices.\n");
          delete ensemble;
          delete tio;
          err++;
          continue;
        }
        //mprintf("DEBUG: Final crd indices arg: %s\n", finalCrdIndicesArg_.c_str());
      }
#   ifdef ENABLE_SINGLE_ENSEMBLE
    } else
      mprintf("Warning: Single ensemble cannot process crdidx.\n");
#   endif
    // Add to ensemble list and update # of frames.
    ensemble_.push_back( ensemble );
    UpdateMaxFrames( ensemble->Traj() );
    delete tio;
  }
  if (err > 0) return 1;
  // FIXME: For backwards compat. overwrite Topology box info with traj box info.
  topIn->SetBoxFromTraj( ensemble_.back()->EnsembleCoordInfo().TrajBox() );
  return 0;
}

int TrajinList::AddTrajin(std::string const& fname, Topology* topIn, ArgList const& argIn)
{
  if (topIn == 0) {
    mprinterr("Error: No topology for input trajectory '%s'\n", fname.c_str());
    return 1;
  }
  // CRDIDXARG
  finalCrdIndicesArg_.clear();
  ArgList trajin_args = argIn;
  bool isRemdtraj = trajin_args.hasKey("remdtraj");
  int err = 0;
  File::NameArray fnames = File::ExpandToFilenames( fname );
  if (fnames.empty()) return 1;
  Trajin* traj = 0;
  for (File::NameArray::const_iterator fn = fnames.begin(); fn != fnames.end(); ++fn) {
    ArgList args = trajin_args;
    if (isRemdtraj)
      traj = new Trajin_Multi();
    else
      traj = new Trajin_Single();
    if (traj == 0) {
      mprinterr("Error: Memory allocation for input trajectory failed.\n");
      return 1;
    }
    traj->SetDebug(debug_);
    if ( traj->SetupTrajRead(*fn, args, topIn) ) {
      mprinterr("Error: Could not set up input trajectory '%s'.\n", fn->full());
      delete traj;
      err++;
      continue;
    }
    // Add to trajin list and update # of frames.
    trajin_.push_back( traj );
    UpdateMaxFrames( traj->Traj() );
  }
  if (err > 0) return 1;
  // FIXME: For backwards compat. overwrite Topology box info with traj box info.
  topIn->SetBoxFromTraj( trajin_.back()->TrajCoordInfo().TrajBox() );
  return 0;
}

void TrajinList::List() const {
  if (!trajin_.empty()) {
    mprintf("\nINPUT TRAJECTORIES (%zu total):\n", trajin_.size());
    unsigned int trajnum = 0;
    for (trajin_it traj = trajin_.begin(); traj != trajin_.end(); ++traj) {
      mprintf(" %u: ", trajnum++);
      (*traj)->PrintInfo( 1 );
    }
  }
  if (!ensemble_.empty()) {
    mprintf("\nINPUT ENSEMBLES (%zu total):\n", ensemble_.size());
    for (unsigned int en = 0; en != ensemble_.size(); ++en) {
      mprintf(" %u: ", en);
      ensemble_[en]->EnsembleInfo( 1 );
    }
  }
  if (maxframes_ < 0)
    mprintf("  Total number of frames is unknown.\n");
  else if (maxframes_ > 0)
    mprintf("  Coordinate processing will occur on %i frames.\n", maxframes_);

}
