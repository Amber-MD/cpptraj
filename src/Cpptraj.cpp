#include <cstdlib> // system
#include "Cpptraj.h"
#include "MpiRoutines.h"
#include "CpptrajStdio.h"
#include "ReadLine.h"

void Cpptraj::Help_List() {
  mprintf("list <type> (<type> = actions,trajin,trajout,ref,parm,analysis,datafile,dataset)\n");
}

void Cpptraj::Help_Help() {
  mprintf("help [<cmd>]\n");
}

void Cpptraj::Help_Debug() {
  mprintf("debug [<type>] <#> ((<type> = actions,trajin,trajout,ref,parm,analysis,datafile)\n");
}

void Cpptraj::Help_ActiveRef() {
  mprintf("activeref <#>\n");
  mprintf("\tSet the reference structure to be used for coordinate-based mask parsing.\n");
}

void Cpptraj::Help_Create_DataFile() {
  mprintf("create <filename> <dataset0> <dataset1> ...\n");
}

void Cpptraj::Help_Precision() {
  mprintf("datafile precision <filename> <dataset> [<width>] [<precision>]\n");
  mprintf("\tIf width/precision not specified default to 12.4\n");
}

enum GeneralCmdTypes { LIST = 0, HELP, QUIT, RUN, DEBUG, NOPROG, NOEXITERR, SYSTEM,
                       ACTIVEREF, READDATA, CREATE, PRECISION, DATAFILE, REFERENCE,
                       TRAJIN, TRAJOUT };

const DispatchObject::Token Cpptraj::GeneralCmds[] = {
  { DispatchObject::GENERAL, "activeref",     0, Help_ActiveRef, ACTIVEREF },
  { DispatchObject::GENERAL, "create",        0, Help_Create_DataFile, CREATE },
  { DispatchObject::GENERAL, "datafile",        0, 0, DATAFILE },
  { DispatchObject::GENERAL, "debug",         0, Help_Debug, DEBUG    },
  { DispatchObject::GENERAL, "go"   ,         0,          0, RUN      },
  { DispatchObject::GENERAL, "help" ,         0,  Help_Help, HELP     },
  { DispatchObject::GENERAL, "list" ,         0,  Help_List, LIST     },
  { DispatchObject::GENERAL, "ls",            0,          0, SYSTEM   },
  { DispatchObject::GENERAL, "noexitonerror", 0,          0, NOEXITERR},
  { DispatchObject::GENERAL, "noprogress",    0,          0, NOPROG   },
  { DispatchObject::GENERAL, "precision",    0,   Help_Precision, PRECISION   },
  { DispatchObject::GENERAL, "prnlev",        0, Help_Debug, DEBUG    },
  { DispatchObject::GENERAL, "pwd",           0,          0, SYSTEM   },
  { DispatchObject::GENERAL, "quit" ,         0,          0, QUIT     },
  { DispatchObject::GENERAL, "readdata",      0,          0, READDATA },
  { DispatchObject::NONE,                  0, 0,          0,      0   }
};

const DispatchObject::Token Cpptraj::CoordCmds[] = {
  { DispatchObject::GENERAL, "reference",     0, FrameList::Help,   REFERENCE },
  { DispatchObject::GENERAL, "trajin",        0, TrajinList::Help,  TRAJIN },
  { DispatchObject::GENERAL, "trajout",       0, TrajoutList::Help, TRAJOUT },
  { DispatchObject::NONE,                  0, 0,                 0, 0 }
};

// Constructor
Cpptraj::Cpptraj() : 
  debug_(0),
  showProgress_(true),
  exitOnError_(true)
{}

void Cpptraj::Help(ArgList& argIn) {
  bool listAllCommands = false;
  ArgList arg = argIn;
  arg.RemoveFirstArg();
  if (arg.empty()) {
    listAllCommands = true;
    mprintf("General Commands:\n");
    SearchTokenArray( GeneralCmds, listAllCommands, arg );
    mprintf("Topology Commands:\n");
    SearchTokenArray( TopologyList::ParmCmds, listAllCommands, arg );
    mprintf("Coordinate Commands:\n");
    SearchTokenArray( CoordCmds, listAllCommands, arg );
    mprintf("Action Commands:\n");
    SearchTokenArray( ActionList::DispatchArray, listAllCommands, arg );
    mprintf("Analysis Commands:\n");
    SearchTokenArray( AnalysisList::DispatchArray, listAllCommands, arg );
  } else {
    if (SearchToken( arg )==0 || dispatchToken_->Help == 0) 
      mprinterr("No help found for %s\n", arg.Command());
    else
      dispatchToken_->Help();
  }
}

void Cpptraj::List(ArgList& argIn) {
  if      (argIn.hasKey("actions")) actionList.List();
  else if (argIn.hasKey("trajin")) trajinList.List();
  else if (argIn.hasKey("ref")) refFrames.List();
  else if (argIn.hasKey("trajout")) trajoutList.List();
  else if (argIn.hasKey("parm")) parmFileList.List();
  else if (argIn.hasKey("analysis")) analysisList.List();
  else if (argIn.hasKey("datafile")) DFL.List();
  else if (argIn.hasKey("dataset")) DSL.List();
  else {
    mprinterr("Error: list: unrecognized list type (%s)\n", argIn.ArgLine());
    Help_List();
  }
}

void Cpptraj::Debug(ArgList& argIn) {
  debug_ = argIn.getNextInteger(0);
  if (argIn.hasKey("actions")) actionList.SetDebug( debug_ );
  else if (argIn.hasKey("trajin")) trajinList.SetDebug( debug_ );
  else if (argIn.hasKey("ref")) refFrames.SetDebug( debug_ );
  else if (argIn.hasKey("trajout")) trajoutList.SetDebug( debug_ );
  else if (argIn.hasKey("parm")) parmFileList.SetDebug( debug_ );
  else if (argIn.hasKey("analysis")) analysisList.SetDebug( debug_ );
  else if (argIn.hasKey("datafile")) DFL.SetDebug( debug_ );
  else if (argIn.hasKey("dataset")) DSL.SetDebug( debug_ );
  else SetGlobalDebug(debug_);
}

// Cpptraj::SetGlobalDebug()
/** Set the debug level for all components of Cpptraj. */
void Cpptraj::SetGlobalDebug(int debugIn) {
  debug_ = debugIn;
  rprintf("GLOBAL DEBUG LEVEL SET TO %i\n",debug_);
  trajinList.SetDebug(debug_);
  refFrames.SetDebug(debug_);
  trajoutList.SetDebug(debug_);
  parmFileList.SetDebug(debug_);
  actionList.SetDebug(debug_);
  analysisList.SetDebug(debug_);
  DFL.SetDebug(debug_);
  DSL.SetDebug(debug_);
}

/// Used to add parm files from the command line.
void Cpptraj::AddParm(const char* parmfile) {
  if (parmfile==NULL) return;
  parmFileList.AddParmFile( parmfile );
}

int Cpptraj::Create_DataFile(ArgList& dataArg) {
  // Next string is datafile that command pertains to.
  std::string name1 = dataArg.GetStringNext();
  if (name1.empty()) {
    mprinterr("Error: create: No filename given.\nError: Usage: ");
    Help_Create_DataFile();
    return 1;
  }
  DataFile* df = DFL.GetDataFile(name1);
  if (df==NULL)
    mprintf("    Creating file %s:",name1.c_str());
  else
    mprintf("    Adding sets to file %s:",name1.c_str());
  int err = 0;
  while ( dataArg.ArgsRemain() ) {
    std::string name2 = dataArg.GetStringNext();
    DataSetList Sets = DSL.GetMultipleSets( name2 );
    if (Sets.empty())
      mprintf("Warning: %s does not correspond to any data sets.\n", name2.c_str());
    for (DataSetList::const_iterator set = Sets.begin(); set != Sets.end(); ++set) {
      mprintf(" %s", (*set)->Legend().c_str());
      if ( DFL.AddSetToFile(name1, *set)==NULL ) {
        mprinterr("Error: Could not add data set %s to file.\n", (*set)->Legend().c_str());
        ++err;
      }
    }
  }
  mprintf("\n");
  return err;
}

int Cpptraj::Precision(ArgList& dataArg) {
  // Next string is datafile that command pertains to.
  std::string name1 = dataArg.GetStringNext();
  if (name1.empty()) {
    mprinterr("Error: precision: No filename given.\nError: Usage: ");
    Help_Precision();
    return 1;
  }
  DataFile* df = DFL.GetDataFile(name1);
  if (df==NULL) {
    mprinterr("Error: precision: DataFile %s does not exist.\n",name1.c_str());
    return 1;
  }
  // This will break if dataset name starts with a digit...
  int width = dataArg.getNextInteger(12);
  int precision = dataArg.getNextInteger(4);
  std::string name2 = dataArg.GetStringNext();
  df->SetPrecision(name2, width, precision);
  return 0;
}

int Cpptraj::ReadData(ArgList& argIn) {
  DataFile dataIn;
  if (dataIn.ReadData( argIn, DSL )!=0) {
    mprinterr("Error: Could not read data file.\n");
    return 1;
  }
  return 0;
}

// Cpptraj::SearchTokenArray()
/** Search the given array for command. If command is found set token
  * and return 1, otherwise return 0.
  */
int Cpptraj::SearchTokenArray(const DispatchObject::Token* DispatchArray,
                              bool listAllCommands, const ArgList& arg)
{
  int col = 0;
  if (listAllCommands) mprintf("\t");
  // List/search Action Commands
  for (const DispatchObject::Token* token = DispatchArray;
                                    token->Type != DispatchObject::NONE; ++token)
  {
    //mprintf("DBG: CMD [%s] LISTALLCMD=%i  ARG=%s\n", token->Cmd,(int)listAllCommands,arg.Command());
    if (listAllCommands) {
      mprintf("%s  ", token->Cmd);
      ++col;
      if (col == 8) {
        mprintf("\n\t");
        col = 0;
      }
    } else if ( arg.CommandIs( token->Cmd ) ) {
      dispatchToken_ = token;
      return 1;
    }
  }
  if (listAllCommands && col > 0) mprintf("\n");
  return 0;
}

// Cpptraj::SearchToken()
/** Search each token list for the given command. If the command is found in
  * a list then dispatchToken is set by SearchTokenArray and 1 is returned.
  * \return 1 if the token is found, 0 if not.
  */
int Cpptraj::SearchToken(const ArgList& argIn) {
  dispatchToken_ = 0;
  // SPECIAL CASE: For backwards compat. remove analyze prefix
  if (argIn.CommandIs("analyze")) {
    ArgList analyzeArg = argIn;
    analyzeArg.RemoveFirstArg();
    if (SearchTokenArray( AnalysisList::DispatchArray, false, analyzeArg)) return 1;
  } else {
    if (SearchTokenArray( GeneralCmds, false, argIn )) return 1;
    if (SearchTokenArray( TopologyList::ParmCmds, false, argIn )) return 1;
    if (SearchTokenArray( CoordCmds, false, argIn )) return 1;
    if (SearchTokenArray( ActionList::DispatchArray, false, argIn)) return 1;
    if (SearchTokenArray( AnalysisList::DispatchArray, false, argIn)) return 1;
  }
  mprinterr("[%s]: Command not found.\n",argIn.Command());
  return 0;
}

// Cpptraj::Interactive()
void Cpptraj::Interactive() {
  ReadLine inputLine;
  // By default when interactive do not exit on errors
  exitOnError_ = false;
  bool readLoop = true;
  while ( readLoop ) {
    inputLine.GetInput(); 
    readLoop = Dispatch( inputLine.c_str() );
  }
}

// Cpptraj::Dispatch()
/** Send commands to their appropriate classes.
 * The command is tried on each class in turn. If the class rejects command
 * move onto the next one. If command is accepted return.
 * \param inputLine NULL-terminated string consisting of commands and arguments.
 */
// NOTE: Should differentiate between keyword rejection and outright error.
bool Cpptraj::Dispatch(const char* inputLine) {
  Topology* tempParm = 0; // For coordinate lists
  //mprintf("\t[%s]\n", inputLine);
  ArgList command( inputLine );
  command.MarkArg(0); // Always mark the command
  if ( SearchToken( command ) ) {
    //mprintf("TOKEN FOUND. CMD=%s  TYPE=%i\n", dispatchToken_->Cmd, (int)dispatchToken_->Type);
    switch (dispatchToken_->Type) {
      case DispatchObject::ACTION : 
        if (actionList.AddAction( dispatchToken_->Alloc, command, &parmFileList,
                                  &refFrames, &DSL, &DFL ) != 0 && exitOnError_ )
          return false; 
        break;
      case DispatchObject::ANALYSIS :
        if ( analysisList.AddAnalysis( dispatchToken_->Alloc, command, &parmFileList, &DSL ) 
             && exitOnError_)
          return false;
        break;
      case DispatchObject::GENERAL :
        switch ( dispatchToken_->Idx ) {
          case TRAJIN :
            tempParm = parmFileList.GetParm(command);
            trajinList.AddTrajin(&command, tempParm);
            break;
          case TRAJOUT :
            tempParm = parmFileList.GetParm(command);
            trajoutList.AddTrajout(&command, tempParm);
            break;
          case REFERENCE :
            tempParm = parmFileList.GetParm(command);
            refFrames.AddReference(&command, tempParm);
            break;
          case LIST  : List(command); break;
          case HELP  : Help(command); break;
          case DEBUG : Debug(command); break;
          case NOPROG: 
            showProgress_ = false; 
            mprintf("\tProgress bar will not be shown.\n");
            break;
          case NOEXITERR:
            exitOnError_ = false;
            mprintf("\tcpptraj will attempt to ignore errors if possible.\n");
            break;
          case ACTIVEREF:
            refFrames.SetActiveRef( command.getNextInteger(0) );
            break;
          case READDATA: 
            if (ReadData( command ) && exitOnError_) return false;
            break;
          case CREATE:
            if (Create_DataFile( command ) && exitOnError_) return false;
            break;
          case PRECISION:
            if (Precision( command ) && exitOnError_) return false;
            break;
          case DATAFILE:
            if (DFL.ProcessDataFileArgs( command ) && exitOnError_) return false;
            break;
          case SYSTEM  : system( command.ArgLine() ); break;
          case RUN     : Run(); // Fall through to quit
          case QUIT    : return false; break;
        }
        break;
      case DispatchObject::PARM :
        parmFileList.CheckCommand(dispatchToken_->Idx, command);
        break;
      default: mprintf("Dispatch type is currently not handled.\n");
    }
  }
  return true;
/*   // Check if command pertains to datafiles
  if ( dispatchArg.CommandIs("datafile") ) {
    DFL.AddDatafileArg(dispatchArg);
    return;
  }*/
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
  parmFileList.List();

  // Print reference information 
  mprintf("\nREFERENCE COORDS:\n");
  refFrames.List();

  // Output traj
  mprintf("\nOUTPUT TRAJECTORIES:\n");
  trajoutList.List();
 
  // Set max frames in the data set list
  DSL.SetMax(maxFrames); 
  
  // Initialize actions and set up data set and data file list
  //if (actionList.Init( &DSL, &refFrames, &DFL, &parmFileList, exitOnError_)) 
  //  return 1;

  // Set up analysis - checks that datasets are present etc
  //if (analysisList.Setup(&DSL, &parmFileList) > 0 && exitOnError_)
  //  return 1;

  // ========== A C T I O N  P H A S E ==========
  // Loop over every trajectory in trajFileList
  rprintf("BEGIN TRAJECTORY PROCESSING:\n");
  trajinList.Begin();
  while ( (traj = trajinList.NextTraj()) != NULL ) {
    // Open up the trajectory file. If an error occurs, bail 
    if ( traj->BeginTraj(showProgress_) ) {
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
      if (actionList.SetupActions( &CurrentParm )) {
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
  DSL.List();

  // ========== A N A L Y S I S  P H A S E ==========
  analysisList.DoAnalyses(&DFL);

  // DEBUG: DataSets, post-Analysis
  mprintf("\nDATASETS AFTER ANALYSIS:\n");
  DSL.List();

  // ========== D A T A  W R I T E  P H A S E ==========
  // Print Datafile information
  DFL.List();
  // Only Master does DataFile output
  if (worldrank==0)
    DFL.Write();
 
  return 0;
}
