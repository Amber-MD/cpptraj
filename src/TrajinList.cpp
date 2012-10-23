// TrajinList
#include <cstddef> //NULL
#include "TrajinList.h"
#include "CpptrajStdio.h"
//#include "MpiRoutines.h" // worldrank,worldsize

TrajinList::TrajinList() : 
  debug_(0),
  maxframes_(0)
{ 
  currentTraj_ = trajin_.end();
}

TrajinList::~TrajinList() {
  for (std::vector<TrajectoryFile*>::iterator traj = trajin_.begin();
                                              traj != trajin_.end(); traj++)
    delete *traj;
}

void TrajinList::Help() { 
  mprintf("trajin <filename> [start] [stop] [offset] [parm <parmfile> | parmindex <#>]\n");
  mprintf("       [remdtraj remdtrajtemp <T>]\n");
}

// TrajinList::AddTrajin()
/** Add trajectory to the trajectory list as an input trajectory. 
  * Associate the trajectory with one of the parm files in the 
  * TopologyList. 
  */
int TrajinList::AddTrajin(ArgList *argIn, Topology *parmIn) {
  TrajectoryFile *traj = new TrajectoryFile(); 
  if (traj==NULL) {
    mprinterr("Error: TrajinList::Add: Could not allocate memory for traj.\n");
    return 1;
  }
  traj->SetDebug(debug_);
  if ( traj->SetupTrajRead(argIn->GetStringNext(),argIn,parmIn) ) {
    mprinterr("Error: trajin: Could not set up trajectory.\n");
    delete traj;
    return 1;
  }

  // Add to trajectory file list
  trajin_.push_back(traj);
  // Update total # of frames
  int trajFrames = traj->Total_Read_Frames();
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

// TrajinList::Begin()
void TrajinList::Begin() {
  currentTraj_ = trajin_.begin();
}

// TrajinList::NextTraj()
/** /return Trajectory pointed to by iterator, or NULL if no more trajectories.
  */
TrajectoryFile *TrajinList::NextTraj() {
  if (currentTraj_ == trajin_.end()) 
    return NULL;
  TrajectoryFile *trj = *currentTraj_;
  ++currentTraj_;
  return trj;
}

void TrajinList::List() {
  mprintf("\nINPUT TRAJECTORIES:\n");
  for (std::vector<TrajectoryFile*>::iterator traj = trajin_.begin();
                                              traj != trajin_.end(); traj++)
    (*traj)->PrintInfo( 1 );
  if (maxframes_ < 0)
    mprintf("  Total number of frames is unknown.\n");
  else
    mprintf("  Total number of frames is %i.\n", maxframes_);
}
