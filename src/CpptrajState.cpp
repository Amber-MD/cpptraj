#include "CpptrajState.h"
#include "MpiRoutines.h"
#include "CpptrajStdio.h"

// Constructor
CpptrajState::CpptrajState() {
  debug=0;
  showProgress=true;
  activeRef = -1;
  exitOnError = true;
}

// Destructor
CpptrajState::~CpptrajState() {
}

/* CpptrajState::SetGlobalDebug()
 * Set the debug level for all components of CpptrajState.
 */
void CpptrajState::SetGlobalDebug(int debugIn) {
  debug = debugIn;
  rprintf("DEBUG LEVEL SET TO %i\n",debug);
  trajinList.SetDebug(debug);
  referenceList.SetDebug(debug);
  trajoutList.SetDebug(debug);
  parmFileList.SetDebug(debug);
  actionList.SetDebug(debug);
  DFL.SetDebug(debug);
}

// CpptrajState::Dispatch()
/** Send commands to their appropriate classes.
 * The command is tried on each class in turn. If the class rejects command
 * move onto the next one. If command is accepted return.
 * \param inputLine NULL-terminated string consisting of commands and arguments.
 */
// NOTE: Should differentiate between keyword rejection and outright error.
void CpptrajState::Dispatch(char *inputLine) {
  ArgList dispatchArg;
  AmberParm *tempParm; // For coordinate lists

  dispatchArg.SetList(inputLine," "); // Space delimited only?
  //printf("    *** %s ***\n",dispatchArg.ArgLine());
  // First argument is the command
  if (dispatchArg.Command()==NULL) {
    if (debug>0) mprintf("NULL Command.\n");
    return;
  }

  // General commands
  // noprogress: Turn off progress bar when processing trajectories.
  if (dispatchArg.CommandIs("noprogress")) {
    showProgress=false;
    mprintf("    noprogress: Progress bar will not be shown.\n");
    return;
  }
  if (dispatchArg.CommandIs("noexitonerror")) {
    mprintf("    noexitonerror: cpptraj will attempt to ignore errors if possible.\n");
    exitOnError=false;
    return;
  }
  // debug: Set global debug level
  if (dispatchArg.CommandIs("debug")) {
    SetGlobalDebug( dispatchArg.getNextInteger(0) );
    return ;
  }
  // actiondebug: Set actions debug level
  if (dispatchArg.CommandIs("actiondebug")) {
    actionList.SetDebug( dispatchArg.getNextInteger(0) );
    return;
  }
  // trajindebug: Set input trajectory debug level
  if (dispatchArg.CommandIs("trajindebug")) {
    trajinList.SetDebug( dispatchArg.getNextInteger(0) );
    return;
  }
  // trajoutdebug: Set output trajectory debug level
  if (dispatchArg.CommandIs("trajoutdebug")) {
    trajoutList.SetDebug( dispatchArg.getNextInteger(0) );
    return;
  }
  // referencedebug: Set reference trajectory debug level
  if (dispatchArg.CommandIs("referencedebug")) {
    referenceList.SetDebug( dispatchArg.getNextInteger(0) );
    return;
  }
  // parmdebug: Set parm debug level
  if (dispatchArg.CommandIs("parmdebug")) {
    parmFileList.SetDebug( dispatchArg.getNextInteger(0) );
    return;
  }
  // datafiledebug: Set master data file list debug
  if (dispatchArg.CommandIs("datafiledebug")) {
    DFL.SetDebug( dispatchArg.getNextInteger(0) );
    return;
  }

  // Mask Selection.
  if (dispatchArg.CommandIs("select")) {
    tempParm = parmFileList.GetParm(dispatchArg);
    AtomMask *tempMask = new AtomMask();
    tempMask->SetMaskString( dispatchArg.getNextMask() );
    // NOTE: No coords for now, Frame list not set up.
    tempParm->SetupIntegerMask( *tempMask, NULL );
    tempMask->PrintMaskAtoms("Selected");
    delete tempMask;
    return;
  }

  // Check if command pertains to coordinate lists
  // If it does, get a parm based on parm/parmindex keywords in arg list
  if (dispatchArg.CommandIs("trajin")) {
    tempParm = parmFileList.GetParm(dispatchArg);
    trajinList.AddTrajin(NULL, &dispatchArg, tempParm);
    return;
  }
  if (dispatchArg.CommandIs("reference")) {
    tempParm = parmFileList.GetParm(dispatchArg);
    referenceList.AddReference(NULL, &dispatchArg, tempParm);
    return;
  }
  if (dispatchArg.CommandIs("activeref")) {
    activeRef = dispatchArg.getNextInteger(0);
    return;
  }
  if (dispatchArg.CommandIs("trajout")) {
    tempParm = parmFileList.GetParm(dispatchArg);
    trajoutList.AddTrajout(NULL, &dispatchArg, tempParm);
    return;
  }

  // Check if command pertains to a parm file
  if (parmFileList.CheckCommand(&dispatchArg)==0) return;

  // Check if command pertains to datafiles
  if ( dispatchArg.CommandIs("datafile") ) {
    DFL.DF_Args.push_back(dispatchArg);
    return;
  }

  // Check if command pertains to an action
  if ( actionList.AddAction(dispatchArg)==0 ) return;

  // Check if command pertains to analysis
  if ( analysisList.AddAnalysis(dispatchArg)==0 ) return; 

  mprintf("Warning: Unknown Command %s.\n",dispatchArg.Command());
}

// CpptrajState::Run()
/** Process trajectories in trajinList. Each frame in trajinList is sent
 *  to the actions in actionList for processing.
 */
int CpptrajState::Run() {
  TrajectoryFile* traj;
  int maxFrames=0;        // Total # of frames that will be read
  int actionSet=0;        // Internal data frame
  int readSets=0;         // Number of frames actually read
  int lastPindex=-1;      // Index of the last loaded parm file
  AmberParm *CurrentParm=NULL; // Parm for actions; can be modified 
  Frame *CurrentFrame=NULL;    // Frame for actions; can be modified
  Frame TrajFrame;       // Original Frame read in from traj
  FrameList refFrames;         // List of reference frames from referenceList

  // ========== S E T U P   P H A S E ========== 
  // Calculate frame division among trajectories
  mprintf("\nINPUT TRAJECTORIES:\n");
  maxFrames=trajinList.SetupFrames();
  trajinList.Info(1,0);
  if (maxFrames<0)  
    mprintf("  Coordinate processing will occur on an unknown number of frames.\n");
  else
    mprintf("  Coordinate processing will occur on %i frames.\n",maxFrames);

  // Parameter file information
  parmFileList.Print();

  // Setup reference frames if reference files were specified 
  referenceList.SetupRefFrames(&refFrames);
  if (activeRef!=-1)
    refFrames.SetActiveRef(activeRef);
  refFrames.Info();

  // Output traj
  mprintf("\nOUTPUT TRAJECTORIES:\n");
  trajoutList.Info(1,0);
 
  // Set max frames in the data set list
  DSL.SetMax(maxFrames); 
  
  // Initialize actions and set up data set and data file list
  if (actionList.Init( &DSL, &refFrames, &DFL, &parmFileList, exitOnError)) 
    return 1;

  // Set up analysis - checks that datasets are present etc
  analysisList.Setup(&DSL, &parmFileList);

  // ========== R U N  P H A S E ==========
  // Loop over every trajectory in trajFileList
  rprintf("BEGIN TRAJECTORY PROCESSING:\n");
  trajinList.Begin();
  while ( (traj = trajinList.NextTraj()) != NULL ) {
    // Open up the trajectory file. If an error occurs, bail 
    if ( traj->BeginTraj(showProgress) ) {
      mprinterr("Error: Could not open trajectory %s.\n",traj->TrajName());
      break;
    }
    // Set current parm from current traj.
    CurrentParm = traj->TrajParm();

    // If Parm has changed, reset Frame and actions for new topology.
    if (lastPindex != CurrentParm->pindex) {
      // Set up the incoming trajectory frame for this parm
      TrajFrame.SetupFrameV(CurrentParm->natom, CurrentParm->mass, traj->HasVelocity());
      // Set up actions for this parm
      if (actionList.Setup( &CurrentParm )) {
        mprintf("WARNING: Could not set up actions for %s: skipping.\n",
                CurrentParm->parmName);
        continue;
      }
      //fprintf(stdout,"DEBUG: After setup of Actions in CpptrajState parm name is %s\n",
      //        CurrentParm->parmName);
      lastPindex = CurrentParm->pindex;
    }

    // Loop over every Frame in trajectory
    traj->PrintInfoLine();
    while ( traj->GetNextFrame(TrajFrame) ) {
      // Since Frame can be modified by actions, save original and use CurrentFrame
      CurrentFrame = &TrajFrame;
      // Perform Actions on Frame
      bool suppress_output = actionList.DoActions(&CurrentFrame, actionSet);
      // Do Output
      if (!suppress_output)
        trajoutList.Write(actionSet, CurrentParm, CurrentFrame);
//#ifdef DEBUG
//      dbgprintf("\tDEBUG: %30s: %4i\n",CurrentParm->parmName,CurrentParm->outFrame);
//#endif
      // Increment frame counters
      actionSet++; 
    }

    // Close the trajectory file
    traj->EndTraj();
    // Update how many frames across all threads have been written for parm
    // Do this for the original parm since it may have been modified.
    //traj->TrajParm()->outFrame += traj->Total_Read_Frames();
    readSets+=traj->NumFramesProcessed();
    mprintf("\n");
  } // End loop over trajin
  rprintf("Read %i frames and processed %i frames.\n",readSets,actionSet);

  // Close output traj
  trajoutList.Close();

  // Do action output
  actionList.Print( );

  // Print Dataset information
  DSL.Info();

  // DataSet Analysis
  analysisList.Analyze(&DFL);

  // Process any datafile commands
  DFL.ProcessDataFileArgs(&DSL);
  // Print Datafile information
  DFL.Info();

  // Do dataset output - first sync datasets
  DSL.Sync();
  // Only Master does DataFile output
  if (worldrank==0)
    DFL.Write();
 
  return 0;
}
