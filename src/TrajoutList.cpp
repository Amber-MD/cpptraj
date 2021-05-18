#include "TrajoutList.h"
#include "Trajout_Single.h"
#include "Topology.h"
#include "CpptrajStdio.h"
#include "ArgList.h"
// NOTE for AddTrajout(), pytraj compatibility
#include "DataSetList.h"

/// CONSTRUCTOR
TrajoutList::TrajoutList() : debug_(0) {}

/// DESTRUCTOR
TrajoutList::~TrajoutList() { Clear(); }

/** Set the list debug level; passed to new output trajectories. */
void TrajoutList::SetDebug(int d) { debug_ = d; }

/** Clear the output trajectory list, free all memory. */
void TrajoutList::Clear() {
  for (ListType::iterator traj = trajout_.begin(); traj != trajout_.end(); ++traj)
    delete *traj;
  trajout_.clear();
  active_.clear();
  trajoutTops_.clear();
  open_.clear();
}

// NOTE This is a legacy function to maintain pytraj compatibility
int TrajoutList::AddTrajout(std::string const& filename, ArgList const& argIn, Topology* tParm)
{
  DataSetList blank;
  return AddTrajout(filename, argIn, blank, tParm);
}

/** Add output trajectory to list as single output trajectory. Associate it
  * with the given Topology but no Topology-dependent setup will occur. This
  * is because during the course of a Run the Topology may be modified, by
  * e.g. a 'strip' command.
  */
int TrajoutList::AddTrajout(std::string const& filename, ArgList const& argIn, DataSetList const& DSLin, Topology* tParm)
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
  if (to->InitTrajWrite(filename, args, DSLin, TrajectoryFile::UNKNOWN_TRAJ)) {
    mprinterr("Error: Could not set up output trajectory.\n");
    delete to;
    return 1;
  }
  trajout_.push_back( to );
  trajoutTops_.push_back( tParm ); 
  open_.push_back( false );
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
          mprinterr("Error: Setting up output trajectory '%s'\n",
                     trajout_[i]->Traj().Filename().full());
          return 1;
        }
        open_[i] = true;
      }
      active_.push_back( trajout_[i] );
    }
  }
  ListActive();
  return 0;
}

/** List only active output trajectories. */
void TrajoutList::ListActive() const {
  if (!trajout_.empty()) {
    mprintf(".....................................................\n");
    if (!active_.empty()) {
      mprintf("ACTIVE OUTPUT TRAJECTORIES (%zu):\n", active_.size());
      for (ListType::const_iterator it = active_.begin(); it != active_.end(); ++it)
      {
        mprintf("  %s", (*it)->Traj().Filename().full());
        std::string meta = (*it)->Traj().CoordInfo().InfoString();
        if (!meta.empty()) mprintf(" (%s)", meta.c_str());
        mprintf("\n");
      }
    } else
      mprintf("NO ACTIVE OUTPUT TRAJECTORIES.\n");
  }
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
        trajout_[i]->SetTrajComm( commIn );
        if ( trajout_[i]->SetupTrajWrite( CurrentParm, cInfo, Nframes) )
        {
          mprinterr("Error: Setting up output trajectory '%s' in parallel.\n",
                    trajout_[i]->Traj().Filename().full());
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
  ListActive();
  return 0;
}
#endif
