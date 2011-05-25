// ReferenceList
#include "ReferenceList.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
ReferenceList::ReferenceList() {
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
int ReferenceList::Add(ArgList *A, ParmFileList *parmFileList) {
  TrajFile *T;
  int startArg,stopArg,offsetArg;
  bool average = false;

  // Set up common arguments from arglist
  if (this->ProcessArgList(A,parmFileList)) return 1;

  // Check if we want to obtain the average structure
  average = A->hasKey("average");

  // Set up basic file to determine type and format
  T = this->SetupTrajectory(trajfilename, fileAccess, UNKNOWN_FORMAT, UNKNOWN_TYPE);

  if (T==NULL) {
    rprintf("ERROR: Setting up file for trajectory %s\n",trajfilename);
    return 1;
  }

  // Set parameter file
  T->P=P;

  // Set up trajectory. 
  if ( T->SetupRead() ) {
    rprintf("ERROR: Setting up %s for read.\n",trajfilename);
    delete T;
    return 1;
  }
  // Get user-specified start arg
  // NOTE: For compatibility with ptraj start from 1
  startArg=A->getNextInteger(1);
  stopArg=startArg;
  offsetArg=1;
  // Get user-specified stop and offset only if getting avg structure
  if (average) {
    stopArg=A->getNextInteger(-1);
    offsetArg=A->getNextInteger(1);
  }
  T->SetArgs(startArg,stopArg,offsetArg);

  // Add to trajectory file list
  this->push_back(T);
  Average.push_back(average); 

  return 0;
}

/* 
 * ReferenceList::SetupRefFrames()
 * Only called for reference trajectories. 
 * Special case of setup frames. Only want one frame from each trajectory,
 * place that frame in refFrames.
 */
int ReferenceList::SetupRefFrames(FrameList *refFrames) {
  int trajFrames, global_set;
  double Nframes;
  Frame *F, *AvgFrame;
  int skipValue;
  int refTrajNum = 0;

  mprintf("\nREFERENCE COORDS:\n");
  if (this->empty()) {
    mprintf("  No reference coordinates.\n");
    return 1;
  }

  for (it = this->begin(); it != this->end(); it++) {
    // Setup the reference traj for reading. Should only be 1 frame.
    // NOTE: For MPI, calling setupFrameInfo with worldrank 0, worldsize 1 for 
    //       all ranks. This is to ensure each thread has a copy of the ref 
    //       struct.
    //       Calling setupFrameInfo with -1 to ensure the Parm frame count is
    //       not updated.

    trajFrames=(*it)->setupFrameInfo(-1,0,1);
    if ((*it)->total_read_frames<1) {
      rprintf("Error: No frames could be read for reference %s, skipping\n",
      (*it)->trajfilename);
      continue;
    }
    if ((*it)->P==NULL) {
      rprintf("Error: No parm associated with ref coords from %s - ignoring.\n",
              (*it)->trajfilename);
      return 1;
    }
    // Reset skip flag if set since all threads need reference coords
    skipValue=0;
    if ((*it)->skip) {
      skipValue=(*it)->skip;
      (*it)->skip=0;
    }
    // Start trajectory read
    global_set=0;
    (*it)->Begin(&global_set, 0);
    (*it)->PrintInfo(1);
    // If averaging requested, loop over specified frames and avg coords.
    if (Average[refTrajNum++]) {
      mprintf("    Averaging over %i frames.\n",trajFrames);
      AvgFrame = new Frame((*it)->P->natom, (*it)->P->mass);
      AvgFrame->ZeroCoords();
      global_set = 0;
      Nframes = 0.0;
      while ( (*it)->NextFrame(&global_set) ) {
        AvgFrame->AddCoord( (*it)->F );
        Nframes++;
      }
      if (Nframes < 1.0) { 
        mprintf("Error: reference average: # frames read is less than 1.\n");
        F=NULL;
      } else {
        AvgFrame->Divide( Nframes );
        F=AvgFrame->Copy();
      }
      delete AvgFrame; 
    // If no averaging, get and copy the 1 frame from Traj, then close
    } else {
      (*it)->NextFrame(&trajFrames);
      F=(*it)->F->Copy();
    }
    // DEBUG
    //fprintf(stdout,"DEBUG: Ref Coord Atom 0\n");
    //F->printAtomCoord(0);
    // Associate this frame with the correct parmfile
    //F->P=(*it)->P;
    // NOTE: Also use full file path??
    refFrames->Add(F,(*it)->File->basefilename,(*it)->P,(*it)->Start());
    (*it)->End();
    // Restore skip value
    (*it)->skip=skipValue;
  }
  return 0;
}


