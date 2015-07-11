#include "TrajoutList.h"
#include "CpptrajStdio.h"
#include "EnsembleOut_Single.h"
#include "EnsembleOut_Multi.h"

void TrajoutList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_ > 0)
    mprintf("TrajoutList debug level set to %i\n", debug_);
}

void TrajoutList::Clear() {
  trajout_.clear();
  trajoutArgs_.clear();
  trajoutTops_.clear();
  trajoutNames_.clear();
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
    if ( to->TrajFilename().Full() == filename ) {
      mprinterr("Error: Output trajectory filename %s already in use.\n",filename.c_str());
      return 1;
    }
  }
  // Create Trajout_Single
  trajout_.push_back( Trajout_Single() );
  // Initialize output trajectory
  ArgList args = argIn;
  if (trajout_.back().InitTrajWrite(filename, args, TrajectoryFile::UNKNOWN_TRAJ)) {
    mprinterr("Error: Could not set up output trajectory.\n");
    trajout_.resize( trajout_.size() - 1 );
    return 1;
  }
  // For potentially setting up ensemble later, save trajout arg.
  trajoutArgs_.push_back( argIn );
  trajoutTops_.push_back( tParm );
  trajoutNames_.push_back( filename );
  return 0;
}

TrajoutList::EnsembleArray TrajoutList::MakeEnsembleTrajout() const {
  EnsembleArray ensembleList;
  for (unsigned int i = 0; i != trajoutArgs_.size(); i++) {
    ArgList argIn = trajoutArgs_[i];
    EnsembleOut* ens = 0;
#   ifdef ENABLE_SINGLE_ENSEMBLE
    // See if single ensemble output desired.
    if (argIn.hasKey("ensemble"))
      ens = new EnsembleOut_Single();
    else
#   endif
      // Create new multi output trajectory
      ens = new EnsembleOut_Multi();
    if (ens == 0) return EnsembleArray();;
    if (ens->InitEnsembleWrite(trajoutNames_[i], argIn,
                               trajoutTops_[i]->ParmCoordInfo().EnsembleSize(),
                               trajout_[i].WriteFormat()))
    {
      delete ens;
      return EnsembleArray();
    }
    ensembleList.push_back( ens );
  }
  return ensembleList;
}
