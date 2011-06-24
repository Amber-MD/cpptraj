// ReferenceList
#include "ReferenceList.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
ReferenceList::ReferenceList() {
  fileAccess=READ;
}

// DESTRUCTOR
ReferenceList::~ReferenceList() { }

/* ReferenceList::Add()
 * Add trajectory to the trajectory list as a reference trajectory. The list
 * will be converted to a list of reference frames by SetupRefFrames before
 * trajectories are processed. Associate the trajectory with one of the parm 
 * files in the ParmFileList. 
 * reference <filename> [start] [parm <parmfile> | parmindex <#>]
 */
int ReferenceList::Add(char *filename, ArgList *A, AmberParm *parmIn) {
  TrajectoryFile *traj;
  bool average = false;

  traj = new TrajectoryFile();
  if (traj==NULL) {
    mprinterr("Error: ReferenceList::Add: Could not allocate memory for traj.\n");
    return 1;
  }
  traj->SetDebug(debug);
  // Check if we want to obtain the average structure
  average = A->hasKey("average");
  // Set up trajectory
  if ( traj->SetupRead(filename,A,parmIn) ) {
    mprinterr("Error: reference: Could not set up trajectory.\n");
    delete traj;
    return 1;
  }
  // If not obtaining average structure, tell trajectory to only process the
  // start frame.
  traj->SingleFrame();

  // Add trajectory and average status 
  this->push_back(traj);
  Average.push_back(average); 

  return 0;
}

/* ReferenceList::SetupRefFrames()
 * Only called for reference trajectories. Convert each traj in the list
 * to either a single frame of that traj or an average over specified frames
 * of the traj. 
 */
int ReferenceList::SetupRefFrames(FrameList *refFrames) {
  std::list<TrajectoryFile*>::iterator traj;
  int trajFrames;
  double Nframes;
  Frame *CurrentFrame, *AvgFrame;
  AmberParm *CurrentParm;
  int refTrajNum = 0;

  mprintf("\nREFERENCE COORDS:\n");
  if (this->empty()) {
    mprintf("  No reference coordinates.\n");
    return 1;
  }

  for (traj = this->begin(); traj != this->end(); traj++) {
    // Setup the reference traj for reading. Should only be 1 frame
    // if not averaged.
    // NOTE: For MPI, calling setupFrameInfo with worldrank 0, worldsize 1 for 
    //       all ranks. This is to ensure each thread has a copy of the ref 
    //       struct.
    //       Calling setupFrameInfo with -1 to ensure the Parm frame count is
    //       not updated.

    trajFrames=(*traj)->SetupFrameInfo();
    if ((*traj)->Total_Read_Frames()<1) {
      mprinterr("Error: No frames could be read for reference %s, skipping\n",
                (*traj)->TrajName());
      continue;
    }
    // Start trajectory read
    if ( (*traj)->BeginTraj(false) ) {
      mprinterr("Error: Could not open reference %s\n.",(*traj)->TrajName());
      return 1;
    }
    (*traj)->PrintInfo(1);
    CurrentParm = (*traj)->TrajParm();
    CurrentFrame = new Frame(CurrentParm->natom, CurrentParm->mass);
    // If averaging requested, loop over specified frames and avg coords.
    if (Average[refTrajNum++]) {
      mprintf("    Averaging over %i frames.\n",trajFrames);
      AvgFrame = new Frame(CurrentParm->natom, CurrentParm->mass);
      AvgFrame->ZeroCoords();
      Nframes = 0.0;
      while ( (*traj)->GetNextFrame(CurrentFrame->X, CurrentFrame->box, &(CurrentFrame->T)) ) {
        AvgFrame->AddCoord( CurrentFrame );
        Nframes++;
      }
      if (Nframes < 1.0) { 
        mprintf("Error: reference average: # frames read is less than 1.\n");
        delete AvgFrame;
        AvgFrame=NULL;
      } else {
        AvgFrame->Divide( Nframes );
      }
      delete CurrentFrame;
      CurrentFrame = AvgFrame;
    // If no averaging, get and copy the 1 frame from Traj, then close
    } else {
      (*traj)->GetNextFrame(CurrentFrame->X, CurrentFrame->box, &(CurrentFrame->T));
    }
    (*traj)->EndTraj();
    // DEBUG
    //fprintf(stdout,"DEBUG: Ref Coord Atom 0\n");
    //F->printAtomCoord(0);
    // Associate this frame with the correct parmfile
    //F->P=(*it)->P;
    if (CurrentFrame==NULL) {
      mprinterr("Error getting frame for %s\n",(*traj)->TrajName());
      return 1;
    }
    // NOTE: Also use full file path??
    refFrames->Add(CurrentFrame,(*traj)->TrajName(),CurrentParm,(*traj)->Start());
  }
  return 0;
}

