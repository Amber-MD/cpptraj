// ReferenceList
#include "ReferenceList.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
ReferenceList::ReferenceList() {
  fileAccess=READ;
}

// DESTRUCTOR
ReferenceList::~ReferenceList() {
  // Free any stripped reference parms since they dont exist in the
  // main parm file list. 
  for (std::vector<AmberParm*>::iterator prm = StrippedRefParms.begin();
                                         prm != StrippedRefParms.end();
                                         prm++)
    delete *prm;
}

// ReferenceList::AddReference()
/** Add trajectory to the trajectory list as a reference trajectory. The list
  * will be converted to a list of reference frames by SetupRefFrames before
  * trajectories are processed. Associate the trajectory with one of the parm 
  * files in the ParmFileList. 
  * reference <filename> [start] [parm <parmfile> | parmindex <#>]
  */
// NOTE: Do not allocate Frames with new, should be static?
int ReferenceList::AddReference(char *filename, ArgList *A, AmberParm *parmIn) {
  TrajectoryFile *traj;
  bool average = false;
  char *maskexpr;
  std::string MaskExpr;

  traj = new TrajectoryFile();
  if (traj==NULL) {
    mprinterr("Error: ReferenceList::Add: Could not allocate memory for traj.\n");
    return 1;
  }
  traj->SetDebug(debug);
  // Check if we want to obtain the average structure
  average = A->hasKey("average");
  // Check for mask expression
  maskexpr = A->getNextMask();
  if (maskexpr!=NULL)
    MaskExpr.assign( maskexpr );
  
  // Set up trajectory
  if ( traj->SetupRead(filename,A,parmIn) ) {
    mprinterr("Error: reference: Could not set up trajectory.\n");
    delete traj;
    return 1;
  }

  // Check for tag after set up to allow SetupRead to look for a parm tag
  std::string reftag = A->getNextTag();

  // If not obtaining average structure, tell trajectory to only process the
  // start frame.
  if (!average)
    traj->SingleFrame();

  // Add trajectory, average status, and mask expression
  trajList.push_back(traj);
  Average.push_back(average); 
  MaskExpressions.push_back( MaskExpr );
  RefTags.push_back( reftag );

  return 0;
}

// ReferenceList::SetupRefFrames()
/** Only called for reference trajectories. Convert each traj in the list
  * to either a single frame of that traj or an average over specified frames
  * of the traj. 
  */
int ReferenceList::SetupRefFrames(FrameList *refFrames) {
  std::list<TrajectoryFile*>::iterator traj;
  int trajFrames;
  double Nframes;
  Frame *CurrentFrame, *AvgFrame;
  AmberParm *CurrentParm;
  int refTrajNum;
  AtomMask Mask;

  mprintf("\nREFERENCE COORDS:\n");
  if (trajList.empty()) {
    mprintf("  No reference coordinates.\n");
    return 1;
  }

  refTrajNum = -1; // Set to -1 since increment is at top of loop
  for (traj = trajList.begin(); traj != trajList.end(); traj++) {
    ++refTrajNum;
    // Setup the reference traj for reading. Should only be 1 frame
    // if not averaged.
    // NOTE: For MPI, calling setupFrameInfo with worldrank 0, worldsize 1 for 
    //       all ranks. This is to ensure each thread has a copy of the ref 
    //       struct.
    //       Calling setupFrameInfo with -1 to ensure the Parm frame count is
    //       not updated.

    trajFrames=(*traj)->Total_Read_Frames();
    if (trajFrames<1) {
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
    // NOTE: If ever need ref velocity change this alloc
    CurrentFrame = new Frame();
    CurrentFrame->SetupFrame(CurrentParm->natom, CurrentParm->mass);
    // If averaging requested, loop over specified frames and avg coords.
    if (Average[refTrajNum]) {
      mprintf("    Averaging over %i frames.\n",trajFrames);
      AvgFrame = new Frame();
      AvgFrame->SetupFrame(CurrentParm->natom, CurrentParm->mass);
      AvgFrame->ZeroCoords();
      Nframes = 0.0;
      while ( (*traj)->GetNextFrame(*CurrentFrame) ) {
        //AvgFrame->AddCoord( CurrentFrame );
        *AvgFrame += *CurrentFrame;
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
    // If no averaging, get and copy the 1 frame from Traj
    } else {
      (*traj)->GetNextFrame(*CurrentFrame);
    }
    (*traj)->EndTraj();
    // DEBUG
    //fprintf(stdout,"DEBUG: Ref Coord Atom 0\n");
    //CurrentFrame->printAtomCoord(0);
    if (CurrentFrame==NULL) {
      mprinterr("Error getting frame for %s\n",(*traj)->TrajName());
      return 1;
    }
    // If a mask expression was specified, strip to match the expression,
    // storing the new parm file here.
    if (!MaskExpressions[refTrajNum].empty()) {
      Mask.SetMaskString((char*)MaskExpressions[refTrajNum].c_str());
      mprintf("    reference: Keeping atoms in mask [%s]\n",Mask.MaskString());
      if (CurrentParm->SetupIntegerMask(Mask, CurrentFrame->X)) return 1;
      if (Mask.None()) {
        mprinterr("Error: No atoms kept for reference.\n");
        return 1;
      }
      // Create new stripped parm
      AmberParm *strippedRefParm = CurrentParm->modifyStateByMask(Mask.Selected,NULL);
      if (strippedRefParm==NULL) {
        mprinterr("Error: could not strip reference.\n");
        return 1;
      }
      strippedRefParm->Summary();
      // Store the new stripped parm in this class so it can be freed later
      StrippedRefParms.push_back(strippedRefParm);
      // Create new stripped frame and set up from current frame
      Frame *strippedRefFrame = new Frame();
      strippedRefFrame->SetupFrame(strippedRefParm->natom,strippedRefParm->mass);
      strippedRefFrame->SetFrameFromMask(CurrentFrame, &Mask);
      // Set the current frame and parm to be stripped frame and parm
      // No need to free CurrentParm since it exists in the parm file list.
      delete CurrentFrame;
      CurrentFrame = strippedRefFrame;
      CurrentParm = strippedRefParm;
    }
    
    // NOTE: Also use full file path??
    refFrames->AddRefFrame(CurrentFrame,(*traj)->TrajName(),CurrentParm,
                           (*traj)->Start(),RefTags[refTrajNum]);
  }
  return 0;
}

