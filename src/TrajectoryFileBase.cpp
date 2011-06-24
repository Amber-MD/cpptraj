// TrajectoryFileBase
#include "TrajectoryFileBase.h"
#include "CpptrajStdio.h"
#include <cstdlib>
#include <cstring>

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
  total_frames=0;
  numFramesRead=0;
  total_read_frames=-1;
}

// DESTRUCTOR
TrajectoryFileBase::~TrajectoryFileBase() {
  if (progress!=NULL) delete progress;
  if (trajName!=NULL) free(trajName);
}

/* TrajectoryFileBase::SetDebug()
 * Set debug level.
 */
void TrajectoryFileBase::SetDebug(int debugIn) {
  debug = debugIn;
  if (debug>0) mprintf("  TrajectoryFile debug level set to %i\n",debug);
}

/*
 * TrajectoryFileBase::SetTrajName()
 * Set trajName to be a copy of input string.
 */
void TrajectoryFileBase::SetTrajName(char *nameIn) {
  trajName=(char*) realloc( trajName, (strlen(nameIn)+1) * sizeof(char) );
  strcpy(trajName,nameIn);
}

/* TrajectoryFileBase::SetArgs()
 * Called after initial trajectory setup, set the start, stop, and offset args.
 * Do some bounds checking.
 * For compatibility with ptraj frames start at 1. So for a traj with 10 frames:
 * cpptraj: 0 1 2 3 4 5 6 7 8 9
 *   ptraj: 1 2 3 4 5 6 7 8 9 10
 * Defaults: startArg=1, stopArg=-1, offsetArg=1
 */ 
void TrajectoryFileBase::SetArgs(int startArg, int stopArg, int offsetArg) {
  //mprintf("DEBUG: setArgs: Original start, stop: %i %i\n",startArg,stopArg);
  if (startArg!=1) {
    if (startArg<1) {
      mprintf("    Warning: %s start argument %i < 1, setting to 1.\n",trajName,startArg);
      start=0; // cpptraj = ptraj - 1
    } else if (total_frames>=0 && startArg>total_frames) {
      // If startArg==stopArg and is greater than # frames, assume we want
      // the last frame (useful when reading for reference structure).
      if (startArg==stopArg) {
        mprintf("    Warning: %s start %i > #Frames (%i), setting to last frame.\n",
                trajName,startArg,total_frames);
        start=total_frames - 1;
      } else {
        mprintf("    Warning: %s start %i > #Frames (%i), no frames will be processed.\n",
                trajName,startArg,total_frames);
        start=startArg - 1;
      }
    } else
      start=startArg - 1;
  }
  if (stopArg!=-1) {
    if ((stopArg - 1)<start) { // cpptraj = ptraj - 1
      mprintf("    Warning: %s stop %i < start, no frames will be processed.\n",
              trajName,stopArg);
      stop = start;
    } else if (total_frames>=0 && stopArg>total_frames) {
      mprintf("    Warning: %s stop %i >= #Frames (%i), setting to max.\n",
              trajName,stopArg,total_frames);
      stop=total_frames;
    } else
      stop=stopArg;
  }

  if (offsetArg!=1) {
    if (offset<1) {
      mprintf("    Warning: %s offset %i < 1, setting to 1.\n",
              trajName,offsetArg);
      offset=1;
    } else if (stop!=-1 && offsetArg > stop - start) {
      mprintf("    Warning: %s offset %i is so large that only 1 set will be processed.\n",
              trajName,offsetArg);
      offset=offsetArg;
    } else
      offset=offsetArg;
  }
  if (debug>0)
    mprintf("  [%s] Args: Start %i Stop %i  Offset %i\n",trajName,start,stop,offset);
}

/* TrajectoryFileBase::SingleFrame()
 * Tell the trajectory to set up stop and offset so that only start frame
 * will be processed.
 */
void TrajectoryFileBase::SingleFrame() {
  stop = start;
  offset = 1;
}

