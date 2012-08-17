// TrajoutList
#include <cstddef> // NULL
#include "TrajoutList.h"
#include "CpptrajStdio.h"
//#include "MpiRoutines.h" //worldsize

TrajoutList::TrajoutList() { }

TrajoutList::~TrajoutList() {
  for (std::vector<TrajectoryFile*>::iterator traj = trajout_.begin(); 
                                              traj != trajout_.end(); traj++) 
    delete *traj;
}

// TrajoutList::AddTrajout()
/** Add trajectory to the trajectory list as an output trajectory. 
  * Associate the trajectory with one of the parm files in the 
  * TopologyList. 
  * trajout <filename> <fileformat> [append] [nobox] [parm <parmfile> | parmindex <#>]
  *         [<range>]
  */
int TrajoutList::AddTrajout(ArgList *argIn, Topology *parmIn) {
  // Since we need to check if this filename is in use in order to prevent
  // overwrites, determine the filename here.
  ArgList::ConstArg filename = argIn->getNextString();
  if (filename==NULL) {
    mprinterr("Error: TrajoutList::Add: Called with NULL filename.\n");
    return 1;
  }
  // Check if filename is in use
  if (FindName(filename) != -1) {
    mprinterr("Error: trajout: Filename %s already in use.\n",filename);
    return 1;
  }

  TrajectoryFile *traj = new TrajectoryFile();
  if (traj==NULL) {
    mprinterr("Error: TrajoutList::Add: Could not allocate memory for traj.\n");
    return 1;
  }
  traj->SetDebug(debug_);
  // Default to AMBERTRAJ; format can be changed via args in the arg list
  if (traj->SetupTrajWrite(filename,argIn,parmIn,TrajectoryFile::UNKNOWN_TRAJ)) {
    mprinterr("Error: trajout: Could not set up trajectory.\n");
    delete traj;
    return 1;
  }

  // Add to trajectory file list
  trajout_.push_back(traj);
  // Add filename to filename list
  AddFilename( filename ); 

  return 0;
}

// TrajoutList::Write()
/** Go through each output traj, call write. The first time CurrentParm
  * matches the parm the trajectory was originally set up with it will
  * be opened, no need to call BeginTraj.
  */ 
int TrajoutList::Write(int set, Topology *CurrentParm, Frame *CurrentFrame) { 
  for (std::vector<TrajectoryFile*>::iterator traj = trajout_.begin(); 
                                              traj != trajout_.end(); traj++) 
  {
    if ( (*traj)->WriteFrame(set, CurrentParm, *CurrentFrame) ) {
      mprinterr("Error writing output trajectory.\n");
      return 1;
    }
  }

  return 0;
}

// TrajoutList::Close()
/** Close output trajectories. Called after input traj processing completed.
  */
void TrajoutList::Close() {
  for (std::vector<TrajectoryFile*>::iterator traj = trajout_.begin(); 
                                              traj != trajout_.end(); traj++) 
    (*traj)->EndTraj();
}

void TrajoutList::Info() {
  if (trajout_.empty()) {
    mprintf("  No files.\n");
  } else {
    for (std::vector<TrajectoryFile*>::iterator traj = trajout_.begin();
                                                traj != trajout_.end(); traj++)
      (*traj)->PrintInfo( 1 );
  }
}

