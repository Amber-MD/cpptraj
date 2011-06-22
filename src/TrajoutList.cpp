// TrajoutList
#include "TrajoutList.h"
#include "CpptrajStdio.h"
#include "PtrajMpi.h" //worldsize

// CONSTRUCTOR
TrajoutList::TrajoutList() {
  fileAccess=WRITE;
}

// DESTRUCTOR
TrajoutList::~TrajoutList() { }

/* 
 * TrajoutList::Add()
 * Add trajectory to the trajectory list as an output trajectory. 
 * Associate the trajectory with one of the parm files in the 
 * ParmFileList. 
 * trajout <filename> <fileformat> [append] [nobox] [parm <parmfile> | parmindex <#>]
 *         [<range>]
 */
int TrajoutList::Add(ArgList *A, AmberParm *parmIn) {
  TrajFile *T;
  FileFormat writeFormat;
  FileType writeType;
  char *onlyframes, *trajfilename; 

  // Filename must be the first argument
  trajfilename = A->getNextString();

  // Must be called with a parm file.
  if (parmIn==NULL) {
    mprinterr("Error: trajout %s: Could not associate with a parm file.\n",trajfilename);
    return 1;
  }

  // Init variables
  writeFormat=AMBERTRAJ; 
  writeType=UNKNOWN_TYPE;
 
  // Check that this filename is not already in use
  if (this->CheckFilename(trajfilename)) {
    mprinterr("Error: trajout: Filename %s has already been used for trajout.\n",
            trajfilename);
    return 1;
  }

  // Set write or append access
  if ( A->hasKey("append") ) fileAccess=APPEND;

  // Set the write file format
  if      ( A->hasKey("pdb")      ) writeFormat=PDBFILE;
  else if ( A->hasKey("data")     ) writeFormat=DATAFILE;
  else if ( A->hasKey("netcdf")   ) writeFormat=AMBERNETCDF;
  else if ( A->hasKey("restart")  ) writeFormat=AMBERRESTART;
  else if ( A->hasKey("ncrestart")) writeFormat=AMBERRESTARTNC;
  else if ( A->hasKey("mol2")     ) writeFormat=MOL2FILE;

  // Set the write file type
  // Since Amber Restart files are written 1 at a time no MPI needed
  if (worldsize>1 && writeFormat!=AMBERRESTART) 
    writeType=MPIFILE;

  // Set up basic file for given type and format
  // If type is unknown it will be determined from extension or will be standard (default)
  T = this->SetupTrajectory(trajfilename, fileAccess, writeFormat, writeType);

  if (T==NULL) {
    rprintf("ERROR: Setting up file for trajectory %s\n",trajfilename);
    return 1;
  }

  // Get specified title if any - will not set if NULL
  T->SetTitle( A->getKeyString("title", NULL) );

  // Process any write arguments
  T->WriteArgs(A);

  // Get a frame range for trajout
  onlyframes = A->getKeyString("onlyframes",NULL);
  if (onlyframes!=NULL) {
    T->FrameRange = new Range();
    if ( T->FrameRange->SetRange(onlyframes) ) {
      mprintf("Warning: trajout: onlyframes: %s is not a valid range.\n");
      delete T->FrameRange;
    } else {
      T->FrameRange->PrintRange("      Saving frames",0);
    }
  }
    
  // Set parameter file
  T->P=parmIn;

  // Set box type from parm file unless "nobox" specified  
  T->BoxType=(int)parmIn->boxType;
  if (A->hasKey("nobox")) T->BoxType=0;

  // No more setup here; Write is set up when first frame written.
  // Add to trajectory file list
  this->push_back(T); 

  return 0;
}

/*
 * TrajoutList::Write()
 * Go through each output traj, set frame from input frame, call write.
 * Only set up an output traj in the TrajList for writing the first time Write
 * is called with the correct (matching) parmtop, based on pindex.
 * If called with outputSet==-2, end write.
 */ 
int TrajoutList::Write(int outputSet, Frame *Fin, AmberParm *CurrentParm) {
  AmberParm *ParmHolder;

  for (it = this->begin(); it != this->end(); it++) {
    // Skip if this input frame parm does not match output frame parm
    if ((*it)->P->pindex!=CurrentParm->pindex) continue;
    // Temporarily store the trajectory parm in ParmHolder in case currentParm is stripped etc
    ParmHolder = (*it)->P;
    (*it)->P = CurrentParm;
    // If there is a framerange defined, check if this frame matches. If so, pop
    if ((*it)->FrameRange!=NULL) {
      // If no more frames in the framerange skip
      if ( (*it)->FrameRange->empty() ) continue;
      // NOTE: For compatibility with ptraj user frame args start at 1
      if ( (*it)->FrameRange->front() - 1 != outputSet ) continue;
      (*it)->FrameRange->pop_front();
    }
    // Open if this is first call - skip flag not otherwise used for output
    if ((*it)->skip==0) {
      if (debug>0) rprintf("    Setting up %s for WRITE (%i atoms, orig. %i atoms)\n",
                           (*it)->trajfilename,(*it)->P->natom,ParmHolder->natom);
      if ((*it)->SetupWrite()) return 1;
      if ((*it)->Begin()) return 1;
      (*it)->skip=1;
    }
    // Set frame and write
    (*it)->F=Fin;
    //fprintf(stdout,"DEBUG: %20s: Writing %i\n",(*it)->File->filename,outputSet);
    if ((*it)->writeFrame(outputSet)) return 1;
    // Dont want this F deallocd since it belongs to input traj or action
    (*it)->F=NULL;
    // Reset the parm
    (*it)->P = ParmHolder;
  }
  return 0;
}

/*
 * TrajoutList::Close()
 * Close output trajectories. Called after input traj processing completed.
 */
void TrajoutList::Close() {
  for (it = this->begin(); it != this->end(); it++) {
    if ( (*it)->skip==1 ) 
      (*it)->End();
  }
}

