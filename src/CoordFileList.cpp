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

/* CoordFileList::SetDebug()
 * Set trajectory list debug level.
 */
void CoordFileList::SetDebug(int debugIn) {
  debug=debugIn;
  if (debug>0)
    mprintf("CoordFileList() DEBUG LEVEL SET TO %i\n",debug);
}

/* CoordFileList::CheckFilename()
 * Check if filenameIn is already in use, return true if so.
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

/* CoordFileList::Info() 
 * Call PrintInfo for each traj in the list.
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

/* CoordFileList::Begin()
 * Set iterator to beginning of trajectory list
 */
void CoordFileList::Begin() {
  currentTraj = trajList.begin();
}

/* CoordFileList::NextTraj()
 * Return the current trajectory in the list and increment the iterator.
 * If at the end of the list return NULL.
 */
TrajectoryFile *CoordFileList::NextTraj() {
  TrajectoryFile *trj;
  if (currentTraj == trajList.end()) return NULL;
  trj = *currentTraj;
  currentTraj++;
  return trj;
}

