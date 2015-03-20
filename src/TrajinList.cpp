// TrajinList
#include "TrajinList.h"
#include "CpptrajStdio.h"
#include "Trajin_Single.h"
#include "Trajin_Multi.h"
#include "Trajin_Ensemble.h"
#include "StringRoutines.h" // ExpandToFilenames

TrajinList::TrajinList() : 
  debug_(0),
  maxframes_(0),
  mode_(UNDEFINED) 
{}

TrajinList::~TrajinList() {
  Clear();
}

void TrajinList::Clear() {
  for (ListType::iterator traj = trajin_.begin(); traj != trajin_.end(); ++traj)
    delete *traj;
  trajin_.clear();
  mode_ = UNDEFINED;
  maxframes_ = 0;
}

// TrajinList::AddToList()
/** Add given Trajin to trajectory file list and update # of frames in topology.
  */
void TrajinList::AddToList(Trajin* traj) {
  trajin_.push_back(traj);
  int trajFrames = traj->TotalReadFrames();
  // If < 0 frames this indicates the number of frames could not be determined. 
  if (trajFrames < 0)
    maxframes_ = -1;
  else if (maxframes_ != -1) {
    // Only update # of frames if total # of frames so far is known.
    traj->TrajParm()->IncreaseFrames( trajFrames );
    maxframes_ += trajFrames;
  }
}

// TrajinList::AddEnsemble()
int TrajinList::AddEnsemble(std::string const& fname, ArgList const& argIn,
                            TopologyList const& topListIn)
{
  if (mode_ == NORMAL) {
    mprinterr("Error: 'ensemble' and 'trajin' are mutually exclusive.\n");
    return 1;
  }
  mode_ = ENSEMBLE;
  int err = 0;
  StrArray fnames = ExpandToFilenames( fname );
  if (fnames.empty()) return 1;
  ArgList trajin_args = argIn;
  // Get parm from TopologyList based on args
  Topology* tempParm = topListIn.GetParm( trajin_args );
  TrajectoryFile::TrajFormatType trajinFmt;
  TrajectoryIO* tio = 0;
  for (StrArray::const_iterator fn = fnames.begin(); fn != fnames.end(); ++fn) {
    ArgList args = trajin_args;
    // Determine whether this file is multiple file or single file ensemble.
    tio = TrajectoryFile::DetectFormat( *fn, trajinFmt );
    if (tio == 0) {
      mprinterr("Error: Could not determine trajectory %s format\n", fn->c_str());
      err++;
      continue;
    }
    Trajin* traj = 0;
#   ifdef ENABLE_SINGLE_ENSEMBLE
    if (tio->CanProcessEnsemble())
      traj = new Trajin_Ensemble();
    else
#   endif
      traj = new Trajin_Multi();
    if (traj == 0) {
      mprinterr("Error: Memory allocation for input trajectory failed.\n");
      delete tio;
      return 1;
    }
    traj->SetDebug(debug_);
    traj->SetEnsemble(true); // TODO: Obsolete; split Trajin_Multi up.
    // CRDIDXARG: Append coordinate indices arg if there is one
    args.AddArg( finalCrdIndicesArg_ );
    if ( traj->SetupTrajRead(fname, args, tempParm) ) {
      mprinterr("Error: Could not set up input trajectory '%s'.\n", fname.c_str());
      delete traj;
      delete tio;
      err++;
      continue;
    }
    // CRDIDXARG: If trajectory is REMD ensemble and sorting by CRDIDX, need to
    //            save final CRDIDX for next ensemble command.
    // TODO: This is very clunky - remlog dataset should contain all exchanges
    //       so trajin doesnt have to worry about it.
#   ifdef ENABLE_SINGLE_ENSEMBLE
    if ( !tio->CanProcessEnsemble() ) {
#   endif
      Trajin_Multi const& mTraj = static_cast<Trajin_Multi const&>( *traj );
      if ( mTraj.TargetMode() == ReplicaInfo::CRDIDX ) {
        finalCrdIndicesArg_ = mTraj.FinalCrdIndices();
        if (finalCrdIndicesArg_.empty()) {
          mprinterr("Error: Could not obtain final remlog indices.\n");
          delete traj;
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
    AddToList(traj);
    //err += AddInputTraj( *fn, traj, argIn, topListIn );
    delete tio;
  }
  if (err > 0) return 1;
  return 0;
}

// TrajinList::AddTrajin()
int TrajinList::AddTrajin(std::string const& fname, ArgList const& argIn,
                          TopologyList const& topListIn)
{
  Trajin* traj = 0;
  if (mode_ == ENSEMBLE) {
    mprinterr("Error: 'trajin' and 'ensemble' are mutually exclusive.\n");
    return 1;
  }
  mode_ = NORMAL;
  finalCrdIndicesArg_.clear();
  ArgList trajin_args = argIn;
  bool isRemdtraj = trajin_args.hasKey("remdtraj");
  int err = 0;
  StrArray fnames = ExpandToFilenames( fname );
  if (fnames.empty()) return 1;
  // Get parm from TopologyList based on args
  Topology* tempParm = topListIn.GetParm( trajin_args );
  for (StrArray::const_iterator fn = fnames.begin(); fn != fnames.end(); ++fn) {
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
    if ( traj->SetupTrajRead(*fn, args, tempParm) ) {
      mprinterr("Error: Could not set up input trajectory '%s'.\n", fname.c_str());
      delete traj;
      err++;
      continue;
    }
    AddToList(traj);
    //err += AddInputTraj( *fn, traj, argIn, topListIn );
  }
  if (err > 0) return 1;
  return 0;
}
/*
// TrajinList::AddInputTraj()
// NOTE: Pass in a copy of ArgList so the arguments arent modified in case
//       multiple files passed to trajin/ensemble.
int TrajinList::AddInputTraj(std::string const& fname, Trajin* traj, ArgList argIn, 
                             TopologyList const& topListIn)
{
  if (traj==0) {
    mprinterr("Error: Could not allocate memory for input trajectory.\n");
    return 1;
  }
  // Get parm from TopologyList based on args
  Topology* tempParm = topListIn.GetParm( argIn );
  traj->SetDebug(debug_);
  // CRDIDXARG: Append coordinate indices arg if there is one
  argIn.AddArg( finalCrdIndicesArg_ ); 
  if ( traj->SetupTrajRead(fname, argIn, tempParm) ) {
    mprinterr("Error: Could not set up input trajectory '%s'.\n", fname.c_str());
    delete traj;
    return 1;
  }
  // CRDIDXARG: If trajectory is REMD ensemble and sorting by CRDIDX, need to
  //            save final CRDIDX for next ensemble command.
  if ( traj->IsEnsemble() ) {
    Trajin_Multi const& mTraj = static_cast<Trajin_Multi const&>( *traj );
    if ( mTraj.TargetMode() == Trajin_Multi::CRDIDX ) {
      finalCrdIndicesArg_ = mTraj.FinalCrdIndices();
      if (finalCrdIndicesArg_.empty()) {
        mprinterr("Error: Could not obtain final remlog indices.\n");
        delete traj;
        return 1;
      }
      //mprintf("DEBUG: Final crd indices arg: %s\n", finalCrdIndicesArg_.c_str());
    }
  } else
    finalCrdIndicesArg_.clear();
  // Add to trajectory file list
  trajin_.push_back(traj);
  // Update total # of frames
  int trajFrames = traj->TotalReadFrames();
  // If < 0 frames this indicates the number of frames could not be determined. 
  if (trajFrames < 0) 
    maxframes_ = -1;
  else if (maxframes_ != -1) {
    // Only update # of frames if total # of frames so far is known.
    traj->TrajParm()->IncreaseFrames( trajFrames );
    maxframes_ += trajFrames;
  }

  return 0;
}
*/
void TrajinList::List() const {
  if (!trajin_.empty()) {
    mprintf("\nINPUT TRAJECTORIES:\n");
    unsigned int trajnum = 0;
    for (ListType::const_iterator traj = trajin_.begin(); traj != trajin_.end(); ++traj) {
      mprintf(" %u: ", trajnum++);
      (*traj)->PrintInfo( 1 );
    }
    if (maxframes_ < 0)
      mprintf("  Total number of frames is unknown.\n");
    else
      mprintf("  Coordinate processing will occur on %i frames.\n", maxframes_);
  }
}
