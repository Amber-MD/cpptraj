// TrajectoryFileBase
#include "TrajectoryFileBase.h"
#include "CpptrajStdio.h"
#include <cstddef> //NULL

// CONSTRUCTOR
TrajectoryFileBase::TrajectoryFileBase() {
  debug = 0;
  progress=NULL;
  trajName=NULL;
  trajParm=NULL;
  fileAccess=READ;
  start=1;
  stop=-1;
  offset=1;
  Frames=0;
  numFramesRead=0;
  total_read_frames=-1;
}

// DESTRUCTOR
TrajectoryFileBase::~TrajectoryFileBase() {
  if (progress!=NULL) delete progress;
}

/* TrajectoryFileBase::SetDebug()
 * Set debug level.
 */
void TrajectoryFileBase::SetDebug(int debugIn) {
  debug = debugIn;
  if (debug>0) mprintf("  TrajectoryFile debug level set to %i\n",debug);
}

