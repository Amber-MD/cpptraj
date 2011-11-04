// CoordFileList
#include "CoordFileList.h"
#include "CpptrajStdio.h"
#include <cstddef> // NULL

// CONSTRUCTOR
CoordFileList::CoordFileList() {
  debug=0;
}

// DESTRUCTOR
CoordFileList::~CoordFileList() {
  std::list<TrajectoryFile*>::iterator traj;
  //fprintf(stderr,"CoordFileList destructor\n");
  for (traj = trajList.begin(); traj != trajList.end(); traj++)
    delete *traj;
}

// CoordFileList::SetDebug()
/** \param debugIn debug level to set.
  */
void CoordFileList::SetDebug(int debugIn) {
  debug=debugIn;
  if (debug>0)
    mprintf("CoordFileList() DEBUG LEVEL SET TO %i\n",debug);
}

// CoordFileList::CheckFilename()
/** \param filenameIn filename to check
  * \return true if filenameIn is already in the list.
  */
bool CoordFileList::FilenameInUse(char *filenameIn) {
  std::list<TrajectoryFile*>::iterator traj;
  if (filenameIn==NULL) {
    mprinterr("Error: CoordFileList::CheckFilename: Called with NULL filename.\n");
    return 1;
  }
  for (traj = trajList.begin(); traj != trajList.end(); traj++)
    if ( (*traj)->TrajFilenameIs(filenameIn) ) return true;

  return false;
}

// CoordFileList::Info()
/** \param showExtended if true show extended information about the trajectory
  * \param indent number of spaces to indent before printing info. 
  */
void CoordFileList::Info(int showExtended, int indent) {
  std::list<TrajectoryFile*>::iterator traj;
  if (trajList.empty()) 
    mprintf("  No files.\n");
  for (traj = trajList.begin(); traj != trajList.end(); traj++) {
    if (indent>0) mprintf("%*s",indent,"");
    (*traj)->PrintInfo(showExtended);
  }
}

// CoordFileList::Begin()
void CoordFileList::Begin() {
  currentTraj = trajList.begin();
}

// CoordFileList::NextTraj()
/** /return Trajectory pointed to by iterator, or NULL if no more trajectories.
  */
TrajectoryFile *CoordFileList::NextTraj() {
  TrajectoryFile *trj;
  if (currentTraj == trajList.end()) return NULL;
  trj = *currentTraj;
  currentTraj++;
  return trj;
}

