#include "PtrajState.h"
#include "PtrajMpi.h"
#include "CpptrajStdio.h"

// Constructor
PtrajState::PtrajState() {
  TotalErrors=0; 
  debug=0;
  showProgress=true;
}

// Destructor
PtrajState::~PtrajState() {
  for (std::list<ArgList*>::iterator it=DF_Args.begin(); it!=DF_Args.end(); it++) 
    delete (*it); 
}

/* PtrajState::SetGlobalDebug()
 * Set the debug level for all components of PtrajState.
 */
void PtrajState::SetGlobalDebug(int debugIn) {
  debug = debugIn;
  rprintf("DEBUG LEVEL SET TO %i\n",debug);
  trajinList.SetDebug(debug);
  referenceList.SetDebug(debug);
  trajoutList.SetDebug(debug);
  parmFileList.SetDebug(debug);
  ptrajActionList.SetDebug(debug);
  DFL.SetDebug(debug);
}

/* PtrajState::Dispatch()
 * Send commands to their appropriate classes.
 * The command is tried on each class in turn. If the class rejects command
 * move onto the next one. If command is accepted return.
 * NOTE: Should differentiate between keyword rejection and outright error.
 */
void PtrajState::Dispatch(char *inputLine) {
  AtomMask *tempMask;  // For ParmInfo
  AmberParm *tempParm; // For ParmInfo
  ArgList dispatchArg;

  dispatchArg.SetList(inputLine," "); // Space delimited only?
  //printf("    *** %s ***\n",dispatchArg.ArgLine());
  // First argument is the command
  if (dispatchArg.Command()==NULL) {
    if (debug>0) mprintf("NULL Command.\n");
    return;
  }

  // Check if command pertains to coordinate lists
  // If it does, get a parm based on parm/parmindex keywords in arg list
  if (dispatchArg.CommandIs("trajin")) {
    tempParm = parmFileList.GetParm(&dispatchArg);
    trajinList.Add(NULL, &dispatchArg, tempParm);
    return;
  }
  if (dispatchArg.CommandIs("reference")) {
    tempParm = parmFileList.GetParm(&dispatchArg);
    referenceList.Add(NULL, &dispatchArg, tempParm);
    return;
  }
  if (dispatchArg.CommandIs("trajout")) {
    tempParm = parmFileList.GetParm(&dispatchArg);
    trajoutList.Add(NULL, &dispatchArg, tempParm);
    return;
  }

  if (dispatchArg.CommandIs("parm")) {
    parmFileList.Add(dispatchArg.getNextString());
    return;
  }

  if (dispatchArg.CommandIs("noprogress")) {
    showProgress=false;
    mprintf("    noprogress: Progress bar will not be shown.\n");
    return;
  }

  if (dispatchArg.CommandIs("debug")) {
    SetGlobalDebug( dispatchArg.getNextInteger(0) );
    return ;
  }

  if (dispatchArg.CommandIs("parminfo")) {
    if ( (tempParm=parmFileList.GetParm(dispatchArg.getNextInteger(0)))!=NULL ) {
      tempMask = new AtomMask();
      char *maskarg = dispatchArg.getNextMask();
      if (maskarg!=NULL) {
        tempMask->SetMaskString( maskarg );
        tempMask->SetupCharMask( tempParm, debug);
        for (int atom=0; atom < tempParm->natom; atom++) {
          if (tempMask->AtomInCharMask(atom)) tempParm->AtomInfo(atom);
        }
      } else {
        tempParm->Summary();
      }
      delete tempMask;
    }
    return;
  }

  if (dispatchArg.CommandIs("parmbondinfo")) {
    if ( (tempParm=parmFileList.GetParm(dispatchArg.getNextInteger(0)))!=NULL ) {
      tempParm->PrintBondInfo();
    }
    return;
  }

  if (dispatchArg.CommandIs("parmmolinfo")) {
    if ( (tempParm=parmFileList.GetParm(dispatchArg.getNextInteger(0)))!=NULL ) {
      tempParm->PrintMoleculeInfo();
    }
    return;
  }

  if (dispatchArg.CommandIs("bondsearch")) {
    parmFileList.SetBondSearch();
    return;
  }

  if (dispatchArg.CommandIs("molsearch")) {
    parmFileList.SetMolSearch();
    return;
  }

  if (dispatchArg.CommandIs("nobondsearch")) {
    parmFileList.SetNoBondSearch();
    return;
  }

  if (dispatchArg.CommandIs("nomolsearch")) {
    parmFileList.SetNoMolSearch();
    return;
  }

  // Check if command pertains to datafiles
  if ( dispatchArg.CommandIs("datafile") ) {
    DF_Args.push_back(dispatchArg.Copy());
    return;
  }

  // Check if command pertains to an action
  if ( ptrajActionList.Add(&dispatchArg)==0 ) return;

  // Check if command pertains to analysis
  if ( analysisList.Add(&dispatchArg)==0 ) return; 

  mprintf("Warning: Unknown Command %s.\n",dispatchArg.Command());
}

/* PtrajState::ProcessDataFileCmd()
 * Process command relating to data files. All datafile commands have format:
 *   datafile <cmd> <datafile> ...
 */
void PtrajState::ProcessDataFileCmd() {
  char *df_cmd = NULL;
  char *name1 = NULL;
  char *name2 = NULL;
  int width,precision;
  DataFile *df;
  ArgList *dataArg;

  if (DF_Args.empty()) return;
  mprintf("DATAFILE SETUP:\n");

  // Loop through all "datafile" arguments
  for (std::list<ArgList*>::iterator it=DF_Args.begin(); it!=DF_Args.end(); it++) {
    dataArg = (*it);
    // Next string will be the argument passed to datafile
    df_cmd = dataArg->getNextString();
    mprintf("  [%s]\n",dataArg->ArgLine());
    // Next string is datafile that command pertains to
    name1 = dataArg->getNextString();
    if (name1==NULL) {
      mprintf("Error: datafile %s: No filename given.\n",df_cmd);
      continue;
    }
    df = DFL.GetDataFile(name1);

    // datafile create
    // Usage: datafile create <filename> <dataset0> <dataset1> ...
    if ( dataArg->ArgIs(1,"create") ) {
      if (df==NULL)
        mprintf("    Creating file %s\n",name1);
      while ( (name2=dataArg->getNextString())!=NULL ) {
        if ( DFL.Add(name1, DSL.Get(name2))==NULL ) {
          mprintf("Warning: Dataset %s does not exist in main dataset list, skipping.\n",name2);
        }
      }

    // datafile xlabel
    // Usage: datafile xlabel <filename> <label>
    } else if ( dataArg->ArgIs(1,"xlabel") ) {
      if (df==NULL) {
        mprintf("Error: datafile xlabel: DataFile %s does not exist.\n",name1);
        continue;
      }
      df->SetXlabel(dataArg->getNextString());

    // datafile ylabel
    // Usage: datafile ylabel <filename> <label>
    } else if ( dataArg->ArgIs(1,"ylabel") ) {
      if (df==NULL) {
        mprintf("Error: datafile ylabel: DataFile %s does not exist.\n",name1);
        continue;
      }
      df->SetYlabel(dataArg->getNextString());

    // datafile invert
    // Usage: datafile invert <filename>
    } else if ( dataArg->ArgIs(1,"invert") ) {
      if (df==NULL) {
        mprintf("Error: datafile invert: DataFile %s does not exist.\n",name1);
        continue;
      }
      mprintf("    Inverting datafile %s\n",name1);
      df->SetInverted();

    // datafile noxcol
    // Usage: datafile noxcol <filename>
    } else if ( dataArg->ArgIs(1,"noxcol") ) {
      if (df==NULL) {
        mprintf("Error: datafile noxcol: DataFile %s does not exist.\n",name1);
        continue;
      }
      mprintf("    Not printing x column for datafile %s\n",name1);
      df->SetNoXcol();
    
    // datafile precision
    // Usage: datafile precision <filename> <dataset> [<width>] [<precision>]
    //        If width/precision not specified default to 12.4
    } else if ( dataArg->ArgIs(1,"precision") ) {
      if (df==NULL) {
        mprintf("Error: datafile precision: DataFile %s does not exist.\n",name1);
        continue;
      }
      name2 = dataArg->getNextString();
      width = dataArg->getNextInteger(12);
      precision = dataArg->getNextInteger(4);
      df->SetPrecision(name2,width,precision);
    }

  } // END loop over datafile args
}  

/* PtrajState::Run()
 * Process trajectories in trajFileList. Each frame in trajFileList is sent
 * to the actions in ptrajActionList for processing.
 */
int PtrajState::Run() {
  TrajectoryFile* traj;
  int maxFrames=0;        // Total # of frames that will be read
  int actionSet=0;        // Internal data frame
  int readSets=0;         // Number of frames actually read
  int lastPindex=-1;      // Index of the last loaded parm file
  AmberParm *CurrentParm=NULL; 
  Frame *CurrentFrame=NULL;
  Frame *TrajFrame=NULL;
  FrameList refFrames;

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
  refFrames.Info();

  // Output traj
  mprintf("\nOUTPUT TRAJECTORIES:\n");
  trajoutList.Info(1,0);
 
  // Set max frames in the data set list
  DSL.SetMax(maxFrames); 
  
  // Initialize actions and set up data set and data file list
  ptrajActionList.Init( &DSL, &refFrames, &DFL, &parmFileList);

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
      if (TrajFrame!=NULL) delete TrajFrame;
      TrajFrame = new Frame(CurrentParm->natom, CurrentParm->mass, traj->HasVelocity());
      // Set up actions for this parm
      if (ptrajActionList.Setup( &CurrentParm )) {
        mprintf("WARNING: Could not set up actions for %s: skipping.\n",
                CurrentParm->parmName);
        continue;
      }
      //fprintf(stdout,"DEBUG: After setup of Actions in PtrajState parm name is %s\n",
      //        CurrentParm->parmName);
      lastPindex = CurrentParm->pindex;
    }

    // Loop over every Frame in trajectory
    while ( traj->GetNextFrame(TrajFrame->X, TrajFrame->V, TrajFrame->box, &(TrajFrame->T)) ) {
      // Since Frame can be modified by actions, save original and use CurrentFrame
      CurrentFrame = TrajFrame;
      // Perform Actions on Frame
      ptrajActionList.DoActions(&CurrentFrame, actionSet);
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
  if (TrajFrame!=NULL) delete TrajFrame;
  rprintf("Read %i frames and processed %i frames.\n",readSets,actionSet);

  // Close output traj
  trajoutList.Close();

  // Do action output
  ptrajActionList.Print( );

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
