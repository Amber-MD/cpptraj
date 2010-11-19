// TrajinList
#include "TrajinList.h"
#include "RemdTraj.h"

// CONSTRUCTOR
TrajinList::TrajinList() {
  Command="trajin";
  fileAccess=READ;
}

// DESTRUCTOR
TrajinList::~TrajinList() { }

/* 
 * TrajinList::Add()
 * Add trajectory to the trajectory list as an input trajectory. 
 * Associate the trajectory with one of the parm files in the 
 * ParmFileList. 
 * trajin <filename> [start] [stop] [offset] [parm <parmfile> | parmindex <#>]
 */
int TrajinList::Add(ArgList *A, ParmFileList *parmFileList, int worldsize) {
  TrajFile *T;
  RemdTraj *R; // Needed to access non-inherited functions
  int startArg, stopArg, offsetArg;

  // Set up common arguments from arglist
  if (this->ProcessArgList(A,parmFileList)) return coordErr;

  // Determine trajectory format and set up appropriate class
  if (A->hasKey("remdtraj")) {
    // Set up as replica
    R = new RemdTraj();
    R->SetReplicaName(trajfilename, A);
    T = (TrajFile*) R;
  } else
    // Set up basic file to determine type and format
    T = this->SetupTrajectory(trajfilename, fileAccess, UNKNOWN_FORMAT, UNKNOWN_TYPE);

  if (T==NULL) {
    fprintf(stdout,"ERROR: Setting up file for trajectory %s\n",trajfilename);
    return 1;
  }

  // Set parameter file
  T->P=P;
 
  // Set up trajectory. 
  if ( T->SetupRead() ) {
    fprintf(stdout,"ERROR: Setting up %s for read.\n",trajfilename);
    delete T;
    return 1;
  }
  // Get any user-specified start, stop, and offset args
  // NOTE: For compatibility with ptraj start from 1
  startArg=A->getNextInteger(1);
  stopArg=A->getNextInteger(-1);
  offsetArg=A->getNextInteger(1);
  T->SetArgs(startArg,stopArg,offsetArg);

  // Add to trajectory file list
  this->push_back(T);

  return 0;
}

/* 
 * TrajinList::SetupFrames()
 * Only called for input trajectories.
 * Loop over all trajectories and call their setup frames routine to calc
 * actual start and stop and how many frames total will be processed.
 * Return the number of frames to be processed.
 */
int TrajinList::SetupFrames() {
  int maxFrames, trajFrames;

  fprintf(stdout,"\nTRAJECTORIES:\n");

  maxFrames=0;

  for (it = this->begin(); it != this->end(); it++) {
    trajFrames = (*it)->setupFrameInfo(maxFrames);
    if (trajFrames==-1) {
      maxFrames=-1;
    }
    if (maxFrames>=0)
      maxFrames+=trajFrames;
    // Print Trajectory information
    (*it)->PrintInfo(1);
  }

  if (maxFrames==-1)
    fprintf(stdout,"  Coordinate processing will occur until EOF (unknown number of frames).\n");
  else
    fprintf(stdout,"  Coordinate processing will occur on %i frames.\n",maxFrames);
  return maxFrames;
}


