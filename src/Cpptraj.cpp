#include <cstdlib> // system
#include "Cpptraj.h"
#include "MpiRoutines.h"
#include "CpptrajStdio.h"
#include "ReadLine.h"

void Cpptraj::Help_List() {
  mprintf("list <type> (<type> = actions,)\n");
}

void Cpptraj::Help_Help() {
  mprintf("help [<cmd>]\n");
}

void Cpptraj::Help_Debug() {
  mprintf("debug [<type>] <#> ((<type> = actions,trajin,trajout,ref,parm,analysis,datafile)\n");
}

enum GeneralCmds { LIST = 0, HELP, QUIT, RUN, DEBUG };

const DispatchObject::Token Cpptraj::DispatchArray[] = {
  { DispatchObject::GENERAL, "list" , 0,  Help_List, LIST },
  { DispatchObject::GENERAL, "help" , 0,  Help_Help, HELP },
  { DispatchObject::GENERAL, "quit" , 0,          0, QUIT },
  { DispatchObject::GENERAL, "go"   , 0,          0, RUN  },
  { DispatchObject::GENERAL, "debug", 0, Help_Debug, DEBUG},
  { DispatchObject::NONE,         0,  0,          0,     0}
};

// Constructor
Cpptraj::Cpptraj() {
  debug=0;
  showProgress=true;
  exitOnError = true;
}

void Cpptraj::List(ArgList& argIn) {
  if (argIn.hasKey("actions")) actionList.List();
  else {
    mprinterr("Error: list: unrecognized list type (%s)\n", argIn.ArgLine());
    Help_List();
  }
}

void Cpptraj::Help(ArgList& argIn) {
  bool listAllCommands = false;
  ArgList arg = argIn;
  arg.RemoveFirstArg();
  if (arg.empty()) {
    listAllCommands = true;
    mprintf("General Commands:\n");
    SearchTokenArray( DispatchArray, listAllCommands, arg );
    mprintf("Action Commands:\n");
    SearchTokenArray( ActionList::DispatchArray, listAllCommands, arg );
  } else {
    if (SearchToken( arg )==0 || dispatchToken_->Help == 0) 
      mprinterr("No help found for %s\n", arg.Command());
    else
      dispatchToken_->Help();
  }
}

void Cpptraj::Debug(ArgList& argIn) {
  debug = argIn.getNextInteger(0);
  if (argIn.hasKey("actions")) actionList.SetDebug( debug );
  else if (argIn.hasKey("trajin")) trajinList.SetDebug( debug );
  else if (argIn.hasKey("ref")) refFrames.SetDebug( debug );
  else if (argIn.hasKey("trajout")) trajoutList.SetDebug( debug );
  else if (argIn.hasKey("parm")) parmFileList.SetDebug( debug );
  else if (argIn.hasKey("analysis")) analysisList.SetDebug( debug );
  else if (argIn.hasKey("datafile")) DFL.SetDebug( debug );
  else SetGlobalDebug(debug);
}

// Cpptraj::SetGlobalDebug()
/** Set the debug level for all components of Cpptraj. */
void Cpptraj::SetGlobalDebug(int debugIn) {
  debug = debugIn;
  rprintf("DEBUG LEVEL SET TO %i\n",debug);
  trajinList.SetDebug(debug);
  refFrames.SetDebug(debug);
  trajoutList.SetDebug(debug);
  parmFileList.SetDebug(debug);
  actionList.SetDebug(debug);
  analysisList.SetDebug(debug);
  DFL.SetDebug(debug);
}

/// Used to add parm files from the command line.
void Cpptraj::AddParm(const char* parmfile) {
  if (parmfile==NULL) return;
  parmFileList.AddParmFile( parmfile );
}

// Cpptraj::SearchTokenArray()
/** Search the given array for command. If command is found set token
  * and return 1, otherwise return 0.
  */
int Cpptraj::SearchTokenArray(const DispatchObject::Token* DispatchArray,
                              bool listAllCommands, const ArgList& arg)
{
  // List/search Action Commands
  for (const DispatchObject::Token* token = DispatchArray;
                                    token->Type != DispatchObject::NONE; ++token)
  {
    //mprintf("DBG: CMD [%s] HELP=%i LISTALLCMD=%i\n", token->Cmd,(int)help,(int)listAllCommands);
    if (listAllCommands)
       mprintf("\t%s\n", token->Cmd);
    else if ( arg.CommandIs( token->Cmd ) ) {
      dispatchToken_ = token;
      return 1;
    }
  }
  return 0;
}

// Cpptraj::SearchToken()
/** Search each token list for the given command. If the command is found in
  * a list then dispatchToken is set by SearchTokenArray and 1 is returned.
  * \return 1 if the token is found, 0 if not.
  */
int Cpptraj::SearchToken(const ArgList& argIn) {
  dispatchToken_ = 0;
  if (SearchTokenArray( DispatchArray, false, argIn )) return 1;
  if (SearchTokenArray( ActionList::DispatchArray, false, argIn)) return 1;
  mprinterr("[%s]: Command not found.\n",argIn.Command());
  return 0;
}

// Cpptraj::Interactive()
void Cpptraj::Interactive() {
  ReadLine inputLine;
  bool readLoop = true;
  while ( readLoop ) {
    inputLine.GetInput(); 
    mprintf("\t[%s]\n", inputLine.c_str());
    ArgList command( inputLine.c_str() );
    command.MarkArg(0); // Always mark the command
    if ( SearchToken( command ) == 0)
      Dispatch( inputLine.c_str() ); // TODO: Remove this
    else {
      switch (dispatchToken_->Type) {
        case DispatchObject::ACTION : 
          actionList.AddAction( dispatchToken_->Alloc, command ); 
          break;
        case DispatchObject::GENERAL :
          switch ( dispatchToken_->Idx ) {
            case LIST : List(command); break;
            case HELP : Help(command); break;
            case DEBUG: Debug(command); break;
            case RUN  : Run(); // Fall through to quit
            case QUIT : readLoop = false; break;
          }
          break;
        default: mprintf("Dispatch type is currently not handled.\n");
      }
    }
  }
}

// Cpptraj::Dispatch()
/** Send commands to their appropriate classes.
 * The command is tried on each class in turn. If the class rejects command
 * move onto the next one. If command is accepted return.
 * \param inputLine NULL-terminated string consisting of commands and arguments.
 */
// NOTE: Should differentiate between keyword rejection and outright error.
void Cpptraj::Dispatch(const char* inputLine) {
  ArgList dispatchArg;
  Topology *tempParm; // For coordinate lists

  dispatchArg.SetList(inputLine," "); // Space delimited only?
  //printf("    *** %s ***\n",dispatchArg.ArgLine());
  // First argument is the command
  if (dispatchArg.Command()==NULL) {
    if (debug>0) mprintf("NULL Command.\n");
    return;
  }
  // Always mark the first argument.
  dispatchArg.MarkArg(0);

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
  if (dispatchArg.CommandIs("debug") || dispatchArg.CommandIs("prnlev")) {
    SetGlobalDebug( dispatchArg.getNextInteger(0) );
    return ;
  }
  // ls, pwd
  if (dispatchArg.CommandIs("ls") || dispatchArg.CommandIs("pwd")) {
    system( dispatchArg.ArgLine() );
    return;
  }
  // actiondebug: Set actions debug level
  if (dispatchArg.CommandIs("actiondebug")) {
    actionList.SetDebug( dispatchArg.getNextInteger(0) );
    return;
  }
  // analysisdebug: Set analyses debug level
  if (dispatchArg.CommandIs("analysisdebug")) {
    analysisList.SetDebug( dispatchArg.getNextInteger(0) );
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
    refFrames.SetDebug( dispatchArg.getNextInteger(0) );
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
    tempParm->SetupIntegerMask( *tempMask );
    tempMask->PrintMaskAtoms("Selected");
    delete tempMask;
    return;
  }

  // Check if command pertains to coordinate lists
  // If it does, get a parm based on parm/parmindex keywords in arg list
  if (dispatchArg.CommandIs("trajin")) {
    tempParm = parmFileList.GetParm(dispatchArg);
    trajinList.AddTrajin(&dispatchArg, tempParm);
    return;
  }
  if (dispatchArg.CommandIs("reference")) {
    tempParm = parmFileList.GetParm(dispatchArg);
    refFrames.AddReference(&dispatchArg, tempParm);
    return;
  }
  if (dispatchArg.CommandIs("activeref")) {
    refFrames.SetActiveRef( dispatchArg.getNextInteger(0) );
    return;
  }
  if (dispatchArg.CommandIs("trajout")) {
    tempParm = parmFileList.GetParm(dispatchArg);
    trajoutList.AddTrajout(&dispatchArg, tempParm);
    return;
  }

  // Check if command pertains to a parm file
  if (parmFileList.CheckCommand(dispatchArg)==0) return;

  // Check if command pertains to datafiles
  if ( dispatchArg.CommandIs("datafile") ) {
    DFL.AddDatafileArg(dispatchArg);
    return;
  }

  // Check if we are reading sets from a datafile
  if ( dispatchArg.CommandIs("readdata") ) {
    DataFile dataIn;
    if (dataIn.ReadData( dispatchArg, DSL )!=0) {
      mprinterr("Error: Could not read data file.\n");
    }
    return;
  }

  // Check if command pertains to an action
  if ( actionList.AddAction(dispatchArg)==0 ) return;

  // Check if command pertains to analysis
  if ( analysisList.AddAnalysis(dispatchArg)==0 ) return; 

  mprintf("Warning: Unknown Command %s.\n",dispatchArg.Command());
}

// Cpptraj::Run()
/** Process trajectories in trajinList. Each frame in trajinList is sent
 *  to the actions in actionList for processing.
 */
int Cpptraj::Run() {
  TrajectoryFile* traj;
  int maxFrames=0;        // Total # of frames that will be read
  int actionSet=0;        // Internal data frame
  int readSets=0;         // Number of frames actually read
  int lastPindex=-1;      // Index of the last loaded parm file
  Topology *CurrentParm=NULL; // Parm for actions; can be modified 
  Frame *CurrentFrame=NULL;    // Frame for actions; can be modified
  Frame TrajFrame;       // Original Frame read in from traj

  // ========== S E T U P   P H A S E ========== 
  // Calculate frame division among trajectories
  mprintf("\nINPUT TRAJECTORIES:\n");
  maxFrames = trajinList.SetupFrames();
  if (maxFrames<0)  
    mprintf("  Coordinate processing will occur on an unknown number of frames.\n");
  else
    mprintf("  Coordinate processing will occur on %i frames.\n",maxFrames);

  // Parameter file information
  parmFileList.Print();

  // Print reference information 
  mprintf("\nREFERENCE COORDS:\n");
  refFrames.Info();

  // Output traj
  mprintf("\nOUTPUT TRAJECTORIES:\n");
  trajoutList.Info();
 
  // Set max frames in the data set list
  DSL.SetMax(maxFrames); 
  
  // Initialize actions and set up data set and data file list
  if (actionList.Init( &DSL, &refFrames, &DFL, &parmFileList, exitOnError)) 
    return 1;

  // Set up analysis - checks that datasets are present etc
  //if (analysisList.Setup(&DSL, &parmFileList) > 0 && exitOnError)
  //  return 1;

  // ========== A C T I O N  P H A S E ==========
  // Loop over every trajectory in trajFileList
  rprintf("BEGIN TRAJECTORY PROCESSING:\n");
  trajinList.Begin();
  while ( (traj = trajinList.NextTraj()) != NULL ) {
    // Open up the trajectory file. If an error occurs, bail 
    if ( traj->BeginTraj(showProgress) ) {
      mprinterr("Error: Could not open trajectory %s.\n",traj->FullTrajStr());
      break;
    }
    // Set current parm from current traj.
    CurrentParm = traj->TrajParm();
    // Check if parm has changed
    bool parmHasChanged = (lastPindex != CurrentParm->Pindex());

    // If Parm has changed or trajectory velocity status has changed,
    // reset the frame.
    if (parmHasChanged || (TrajFrame.HasVelocity() != traj->HasVelocity()))
      TrajFrame.SetupFrameV(CurrentParm->Atoms(), traj->HasVelocity());

    // If Parm has changed, reset actions for new topology.
    if (parmHasChanged) {
      // Set active reference for this parm
      CurrentParm->SetReferenceCoords( refFrames.ActiveReference() );
      // Set up actions for this parm
      if (actionList.Setup( &CurrentParm )) {
        mprintf("WARNING: Could not set up actions for %s: skipping.\n",
                CurrentParm->c_str());
        continue;
      }
      //fprintf(stdout,"DEBUG: After setup of Actions in Cpptraj parm name is %s\n",
      //        CurrentParm->parmName);
      lastPindex = CurrentParm->Pindex();
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
//#     ifdef DEBUG
//      dbgprintf("\tDEBUG: %30s: %4i\n",CurrentParm->parmName,CurrentParm->outFrame);
//#     endif
      // Increment frame counter
      ++actionSet; 
    }

    // Close the trajectory file
    traj->EndTraj();
    // Update how many frames have been processed.
    readSets += traj->NumFramesProcessed();
    mprintf("\n");
  } // End loop over trajin
  rprintf("Read %i frames and processed %i frames.\n",readSets,actionSet);

  // Close output traj
  trajoutList.Close();

  // ========== A C T I O N  O U T P U T  P H A S E ==========
  actionList.Print( );

  // Sync DataSets and print DataSet information
  // TODO - Also have datafilelist call a sync??
  DSL.Sync();
  mprintf("\nDATASETS BEFORE ANALYSIS:\n");
  DSL.Info();

  // ========== A N A L Y S I S  P H A S E ==========
  analysisList.Setup(&DSL, &parmFileList);
  analysisList.Analyze(&DFL);

  // DEBUG: DataSets, post-Analysis
  mprintf("\nDATASETS AFTER ANALYSIS:\n");
  DSL.Info();

  // ========== D A T A  W R I T E  P H A S E ==========
  // Process any datafile commands
  DFL.ProcessDataFileArgs(&DSL);
  // Print Datafile information
  DFL.Info();
  // Only Master does DataFile output
  if (worldrank==0)
    DFL.Write();
 
  return 0;
}
