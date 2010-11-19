// ReferenceList
#include "ReferenceList.h"

// CONSTRUCTOR
ReferenceList::ReferenceList() {
  Command="reference";
  fileAccess=READ;
}

// DESTRUCTOR
ReferenceList::~ReferenceList() { }

/* 
 * ReferenceList::Add()
 * Add trajectory to the trajectory list as a reference trajectory. The list
 * will be converted to a list of reference frames by SetupRefFrames before
 * trajectories are processed. Associate the trajectory with one of the parm 
 * files in the ParmFileList. 
 * reference <filename> [start] [parm <parmfile> | parmindex <#>]
 */
int ReferenceList::Add(ArgList *A, ParmFileList *parmFileList, int worldsize) {
  TrajFile *T;
  int startArg;

  // Set up common arguments from arglist
  if (this->ProcessArgList(A,parmFileList)) return coordErr;

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
  // Get user-specified start arg
  // NOTE: For compatibility with ptraj start from 1
  startArg=A->getNextInteger(1);
  T->SetArgs(startArg,startArg,1);

  // Add to trajectory file list
  this->push_back(T); 

  return 0;
}

/* 
 * ReferenceList::SetupRefFrames()
 * Only called for reference trajectories. 
 * Special case of setup frames. Only want one frame from each trajectory,
 * place that frame in refFrames.
 */
int ReferenceList::SetupRefFrames(FrameList *refFrames) {
  int trajFrames;
  Frame *F;

  fprintf(stdout,"\nREFERENCE COORDS:\n");
  if (this->empty()) {
    fprintf(stdout,"  No reference coordinates.\n");
    return 1;
  }

  for (it = this->begin(); it != this->end(); it++) {
    // Setup the reference traj for reading. Should only be 1 frame.
    trajFrames=(*it)->setupFrameInfo(-1);
    if ((*it)->total_read_frames<1) {
      fprintf(stdout,"Error: No frames could be read for reference %s, skipping\n",
      (*it)->trajfilename);
      continue;
    }
    if ((*it)->P==NULL) {
      fprintf(stdout,"Error: No parm associated with ref coords from %s - ignoring.\n",
              (*it)->trajfilename);
      return 1;
    }
    // Reset skip flag if set since all threads need reference coords
    if ((*it)->skip)
      (*it)->skip=0;
    (*it)->Begin(&trajFrames, 0);
    // Get and copy the 1 frame from Traj, then close
    // NOTE: What happens when not seekable?
    (*it)->PrintInfo(1);
    (*it)->NextFrame(&trajFrames);
    F=(*it)->F->Copy();
    // DEBUG
    //fprintf(stdout,"DEBUG: Ref Coord Atom 0\n");
    //F->printAtomCoord(0);
    // Associate this frame with the correct parmfile
    //F->P=(*it)->P;
    // NOTE: Also use full file path??
    refFrames->Add(F,(*it)->File->basefilename,(*it)->P);
    (*it)->End();
  }
  return 0;
}


