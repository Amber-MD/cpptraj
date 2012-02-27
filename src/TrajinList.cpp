// TrajinList
#include <cstddef> //NULL
#include "TrajinList.h"
#include "CpptrajStdio.h"
//#include "MpiRoutines.h" // worldrank,worldsize

// CONSTRUCTOR
TrajinList::TrajinList() {
  fileAccess=READ;
}

// TrajinList::AddTrajin()
/** Add trajectory to the trajectory list as an input trajectory. 
  * Associate the trajectory with one of the parm files in the 
  * ParmFileList. 
  * trajin <filename> [start] [stop] [offset] [parm <parmfile> | parmindex <#>]
  *        [remdtraj remdtrajtemp <T>]
  */
int TrajinList::AddTrajin(char *filename, ArgList *A, AmberParm *parmIn) {
  TrajectoryFile *traj;

  traj = new TrajectoryFile(); 
  if (traj==NULL) {
    mprinterr("Error: TrajinList::Add: Could not allocate memory for traj.\n");
    return 1;
  }
  traj->SetDebug(debug);
  if ( traj->SetupRead(filename,A,parmIn) ) {
    mprinterr("Error: trajin: Could not set up trajectory.\n");
    delete traj;
    return 1;
  }

  // Add to trajectory file list
  trajList.push_back(traj);

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
  std::list<TrajectoryFile*>::iterator traj;
  int maxFrames, trajFrames;

  maxFrames=0;

  for (traj = trajList.begin(); traj != trajList.end(); traj++) {
    trajFrames = (*traj)->Total_Read_Frames();
    // If < 0 frames this indicates the number of frames could not be determined. 
    if (trajFrames < 0) return -1;
    (*traj)->TrajParm()->parmFrames += trajFrames;
    maxFrames+=trajFrames;
  }

  return maxFrames;
}

