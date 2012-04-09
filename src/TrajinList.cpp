// TrajinList
#include <cstddef> //NULL
#include "TrajinList.h"
#include "CpptrajStdio.h"
//#include "MpiRoutines.h" // worldrank,worldsize

TrajinList::TrajinList() : debug_(0)
{ 
  currentTraj_ = trajin_.end();
}

TrajinList::~TrajinList() {
  for (std::vector<TrajectoryFile*>::iterator traj = trajin_.begin();
                                              traj != trajin_.end(); traj++)
    delete *traj;
}

// TrajinList::AddTrajin()
/** Add trajectory to the trajectory list as an input trajectory. 
  * Associate the trajectory with one of the parm files in the 
  * TopologyList. 
  * trajin <filename> [start] [stop] [offset] [parm <parmfile> | parmindex <#>]
  *        [remdtraj remdtrajtemp <T>]
  */
int TrajinList::AddTrajin(char *filename, ArgList *A, Topology *parmIn) {
  TrajectoryFile *traj;

  traj = new TrajectoryFile(); 
  if (traj==NULL) {
    mprinterr("Error: TrajinList::Add: Could not allocate memory for traj.\n");
    return 1;
  }
  traj->SetDebug(debug_);
  if ( traj->SetupRead(filename,A,parmIn) ) {
    mprinterr("Error: trajin: Could not set up trajectory.\n");
    delete traj;
    return 1;
  }

  // Add to trajectory file list
  trajin_.push_back(traj);

  return 0;
}

// TrajinList::SetupFrames()
/** Only called for input trajectories.
  * Loop over all trajectories and call their setup frames routine to calc
  * actual start and stop and how many frames total will be processed. Update 
  * the number of frames that will be read for the associated traj parm.
  * Return the total number of frames to be processed across all trajins.
  */
int TrajinList::SetupFrames() {
  std::vector<TrajectoryFile*>::iterator traj;
  int maxFrames = 0;

  for (std::vector<TrajectoryFile*>::iterator traj = trajin_.begin(); 
                                              traj != trajin_.end(); traj++) 
  {
    int trajFrames = (*traj)->Total_Read_Frames();
    // If < 0 frames this indicates the number of frames could not be determined. 
    if (trajFrames < 0) {
      maxFrames = -1;
    } else if (maxFrames != -1) {
      // Only update # of frames if total # of frames so far is known.
      (*traj)->TrajParm()->IncreaseFrames( trajFrames );
      maxFrames += trajFrames;
    }
    // Print input traj information
    (*traj)->PrintInfo( 1 );
  }

  return maxFrames;
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

