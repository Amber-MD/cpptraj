#include "CpptrajState.h"
#include "MpiRoutines.h"
#include "CpptrajStdio.h"

// Constructor
CpptrajState::CpptrajState() {
  TotalErrors=0; 
  debug=0;
  showProgress=true;
  activeRef = -1;
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

/* CpptrajState::Dispatch()
 * Send commands to their appropriate classes.
 * The command is tried on each class in turn. If the class rejects command
 * move onto the next one. If command is accepted return.
 * NOTE: Should differentiate between keyword rejection and outright error.
 */
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
    DF_Args.push_back(dispatchArg);
    return;
  }

  // Check if command pertains to an action
  if ( actionList.AddAction(dispatchArg)==0 ) return;

  // Check if command pertains to analysis
  if ( analysisList.AddAnalysis(dispatchArg)==0 ) return; 

  mprintf("Warning: Unknown Command %s.\n",dispatchArg.Command());
}

/* CpptrajState::ProcessDataFileCmd()
 * Process command relating to data files. All datafile commands have format:
 *   datafile <cmd> <datafile> ...
 */
void CpptrajState::ProcessDataFileCmd() {
  char *df_cmd = NULL;
  char *name1 = NULL;
  char *name2 = NULL;
  int width,precision;
  DataFile *df;

  if (DF_Args.empty()) return;
  mprintf("DATAFILE SETUP:\n");

  // Loop through all "datafile" arguments
  for (std::list<ArgList>::iterator dataArg=DF_Args.begin(); 
                                    dataArg!=DF_Args.end(); 
                                    dataArg++) 
  {
    // Next string will be the argument passed to datafile
    df_cmd = (*dataArg).getNextString();
    mprintf("  [%s]\n",(*dataArg).ArgLine());
    // Next string is datafile that command pertains to
    name1 = (*dataArg).getNextString();
    if (name1==NULL) {
      mprintf("Error: datafile %s: No filename given.\n",df_cmd);
      continue;
    }
    df = DFL.GetDataFile(name1);

    // datafile create
    // Usage: datafile create <filename> <dataset0> <dataset1> ...
    if ( (*dataArg).ArgIs(1,"create") ) {
      if (df==NULL)
        mprintf("    Creating file %s\n",name1);
      while ( (name2=(*dataArg).getNextString())!=NULL ) {
        if ( DFL.Add(name1, DSL.Get(name2))==NULL ) {
          mprintf("Warning: Dataset %s does not exist in main dataset list, skipping.\n",name2);
        }
      }

    // datafile xlabel
    // Usage: datafile xlabel <filename> <label>
    } else if ( (*dataArg).ArgIs(1,"xlabel") ) {
      if (df==NULL) {
        mprintf("Error: datafile xlabel: DataFile %s does not exist.\n",name1);
        continue;
      }
      df->SetXlabel((*dataArg).getNextString());

    // datafile ylabel
    // Usage: datafile ylabel <filename> <label>
    } else if ( (*dataArg).ArgIs(1,"ylabel") ) {
      if (df==NULL) {
        mprintf("Error: datafile ylabel: DataFile %s does not exist.\n",name1);
        continue;
      }
      df->SetYlabel((*dataArg).getNextString());

    // datafile invert
    // Usage: datafile invert <filename>
    } else if ( (*dataArg).ArgIs(1,"invert") ) {
      if (df==NULL) {
        mprintf("Error: datafile invert: DataFile %s does not exist.\n",name1);
        continue;
      }
      mprintf("    Inverting datafile %s\n",name1);
      df->SetInverted();

    // datafile noxcol
    // Usage: datafile noxcol <filename>
    } else if ( (*dataArg).ArgIs(1,"noxcol") ) {
      if (df==NULL) {
        mprintf("Error: datafile noxcol: DataFile %s does not exist.\n",name1);
        continue;
      }
      mprintf("    Not printing x column for datafile %s\n",name1);
      df->SetNoXcol();
    
    // datafile precision
    // Usage: datafile precision <filename> <dataset> [<width>] [<precision>]
    //        If width/precision not specified default to 12.4
    } else if ( (*dataArg).ArgIs(1,"precision") ) {
      if (df==NULL) {
        mprintf("Error: datafile precision: DataFile %s does not exist.\n",name1);
        continue;
      }
      name2 = (*dataArg).getNextString();
      width = (*dataArg).getNextInteger(12);
      precision = (*dataArg).getNextInteger(4);
      df->SetPrecision(name2,width,precision);
    }

  } // END loop over datafile args
}  

/* CpptrajState::Run()
 * Process trajectories in trajFileList. Each frame in trajFileList is sent
 * to the actions in actionList for processing.
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
  actionList.Init( &DSL, &refFrames, &DFL, &parmFileList);

  // Set up analysis - checks that datasets are present etc
  analysisList.Setup(&DSL);

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
    traj->PrintInfoLine();
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
    while ( traj->GetNextFrame(TrajFrame.X, TrajFrame.V, TrajFrame.box, &(TrajFrame.T)) ) {
      // Since Frame can be modified by actions, save original and use CurrentFrame
      CurrentFrame = &TrajFrame;
      // Perform Actions on Frame
      actionList.DoActions(&CurrentFrame, actionSet);
      // Do Output
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
  ProcessDataFileCmd();
  // Print Datafile information
  DFL.Info();

  // Do dataset output - first sync datasets
  DSL.Sync();
  // Only Master does DataFile output
  if (worldrank==0)
    DFL.Write();
 
  return 0;
}
