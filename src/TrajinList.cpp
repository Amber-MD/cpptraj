// TrajinList
#include "TrajinList.h"
#include "CpptrajStdio.h"
#include "Trajin_Single.h"
#include "Trajin_Multi.h"

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

// TrajinList::AddTrajin()
/** Add trajectory to the trajectory list as an input trajectory. 
  * Associate the trajectory with one of the parm files in the 
  * TopologyList. 
  */
int TrajinList::AddTrajin(ArgList& argIn, TopologyList& topListIn) {
  Trajin* traj = 0;
  if ( argIn.CommandIs("trajin") ) {
    // trajin and ensemble are currently mutually exclusive
    if (mode_ == ENSEMBLE) {
      mprinterr("Error: 'trajin' and 'ensemble' are mutually exclusive.\n");
      return 1;
    }
    if (argIn.hasKey("remdtraj"))
      traj = new Trajin_Multi();
    else
      traj = new Trajin_Single();
    mode_ = NORMAL;
  } else if ( argIn.CommandIs("ensemble") ) {
    if (mode_ == NORMAL) {
      mprinterr("Error: 'ensemble' and 'trajin' are mutually exclusive.\n");
      return 1;
    }
    traj = new Trajin_Multi();
    mode_ = ENSEMBLE;
  } else
    return 1; // Command does not pertain to TrajinList
  if (traj==0) {
    mprinterr("Error: TrajinList::Add: Could not allocate memory for traj.\n");
    return 1;
  }
  // Get parm from TopologyList based on args
  Topology* tempParm = topListIn.GetParm( argIn );
  traj->SetDebug(debug_);
  if ( traj->SetupTrajRead(argIn.GetStringNext(), &argIn, tempParm) ) {
    mprinterr("Error: trajin: Could not set up trajectory.\n");
    delete traj;
    return 1;
  }

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

void TrajinList::List() const {
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
