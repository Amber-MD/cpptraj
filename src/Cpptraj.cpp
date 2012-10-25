#include <cstdio> 
#include <cstdlib> // system
#include "Cpptraj.h"
#include "MpiRoutines.h"
#include "CpptrajStdio.h"
#include "Trajin_Multi.h"
#include "FrameArray.h"
#include "ReadLine.h"

void Cpptraj::Usage(const char* programName) {
  mprinterr("\nUsage: %s [-p <Top1>, -p <Top2>, ...] [-i <Input1>, -i <Input2>, ...]\n",
            programName);
  mprinterr(  "       %s <Top> <Input>\n",programName);
  mprinterr(  "       Additional options:\n");
  mprinterr(  "         --help, -help : Print usage information and exit.\n");
  mprinterr(  "         -V, --version : Print version information and exit.\n");
  mprinterr(  "         --defines     : Print list of defines used in compilation.\n");
  mprinterr(  "         -debug <N>    : Set global debug level.\n");
  mprinterr(  "         --interactive : Enter interactive mode.\n");
}

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
  mprintf("precision <filename> [<dataset>] [<width>] [<precision>]\n");
  mprintf("\tSet precision for dataset(s) in given datafile to <width>.<precision>\n");
  mprintf("If width/precision not specified default to 12.4\n");
}

void Cpptraj::Help_SelectDS() {
  mprintf("selectds <dataset selection>\n");
  mprintf("\tShow results of data set selection. Data set selection format is:\n");
  mprintf("\t\t<name>[<aspect]:<idx range>\n");
  mprintf("Where '<name>' is the data set name, '[<aspect>]' is the data set aspect,\n");
  mprintf("and <idx range> is a numerical range specifying data set indices (i.e. 2-5,7 etc).\n");
  mprintf("The aspect and index portions may be optional. An asterisk '*' may be used as\n");
  mprintf("a wildcard. E.g. 'selectds R2', 'selectds RoG[Max]', 'selectds PR[res]:2-12'\n");
}

enum GeneralCmdTypes { LIST = 0, HELP, QUIT, RUN, DEBUG, NOPROG, NOEXITERR, 
                       SYSTEM, ACTIVEREF, READDATA, CREATE, PRECISION, DATAFILE,
                       SELECTDS, READINPUT, RUN_ANALYSIS, WRITEDATA };

const DispatchObject::Token Cpptraj::GeneralCmds[] = {
  { DispatchObject::GENERAL, "activeref",     0, Help_ActiveRef, ACTIVEREF },
  { DispatchObject::GENERAL, "runanalysis",     0, 0, RUN_ANALYSIS },
  { DispatchObject::GENERAL, "create",        0, Help_Create_DataFile, CREATE },
  { DispatchObject::GENERAL, "datafile",        0, 0, DATAFILE },
  { DispatchObject::GENERAL, "debug",         0, Help_Debug, DEBUG    },
  { DispatchObject::GENERAL, "go"   ,         0,          0, RUN      },
  { DispatchObject::GENERAL, "head" ,         0,  0, SYSTEM     },
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
  { DispatchObject::GENERAL, "readinput",      0,          0, READINPUT },
  { DispatchObject::GENERAL, "writedata",      0,          0, WRITEDATA },
  { DispatchObject::GENERAL, "selectds",      0, Help_SelectDS, SELECTDS },
  { DispatchObject::NONE,                  0, 0,          0,      0   }
};

enum CoordCmdTypes { REFERENCE, TRAJIN, TRAJOUT };

const DispatchObject::Token Cpptraj::CoordCmds[] = {
  { DispatchObject::COORD, "reference",     0, FrameList::Help,   REFERENCE },
  { DispatchObject::COORD, "trajin",        0, TrajinList::Help_Trajin,  TRAJIN },
  { DispatchObject::COORD, "ensemble",        0, TrajinList::Help_Ensemble,  TRAJIN },
  { DispatchObject::COORD, "trajout",       0, TrajoutList::Help, TRAJOUT },
  { DispatchObject::NONE,                  0, 0,                 0, 0 }
};

// -----------------------------------------------------------------------------
// Constructor
Cpptraj::Cpptraj() : 
  debug_(0),
  showProgress_(true),
  exitOnError_(true),
  nrun_(0)
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

void Cpptraj::SelectDS(ArgList& argIn) {
  std::string dsarg = argIn.GetStringNext();
  DataSetList dsets = DSL.GetMultipleSets( dsarg );
  mprintf("SelectDS: Arg [%s] selects %i data sets:\n", dsarg.c_str(), dsets.size());
  for (DataSetList::const_iterator set = dsets.begin(); set != dsets.end(); ++set)
    (*set)->Info();
}

// -----------------------------------------------------------------------------
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
int Cpptraj::SearchToken(ArgList& argIn) {
  dispatchToken_ = 0;
  // SPECIAL CASE: For backwards compat. remove analyze prefix
  if (argIn.CommandIs("analyze")) {
    argIn.RemoveFirstArg();
    argIn.MarkArg(0); // Mark new first arg as command
    if (SearchTokenArray( AnalysisList::DispatchArray, false, argIn)) return 1;
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

// -----------------------------------------------------------------------------
// Cpptraj::Interactive()
void Cpptraj::Interactive() {
  ReadLine inputLine;
  // By default when interactive do not exit on errors
  exitOnError_ = false;
  Mode readLoop = C_OK;
  while ( readLoop == C_OK ) {
    inputLine.GetInput(); 
    readLoop = Dispatch( inputLine.c_str() );
  }
}

static inline bool EndChar(char ptr) {
  if (ptr=='\n' || ptr=='\r' || ptr=='\0' || ptr==EOF) return true;
  return false;
}

/** Read commands from an input file. '#' indicates the beginning of a
  * comment, backslash at the end of a line indicates continuation
  * (otherwise indicates 'literal').
  * \return 0 if successfully read, 1 on error.
  */
int Cpptraj::ProcessInput(std::string const& inputFilename) {
  FILE *infile;
  if (inputFilename.empty()) return 1;
  mprintf("INPUT: Reading Input from file %s\n",inputFilename.c_str());
  if ( (infile=fopen(inputFilename.c_str(),"r"))==NULL ) {
    rprintf("Error: Could not open input file %s\n",inputFilename.c_str());
    return 1;
  }
  // Read in each line of input. Newline or NULL terminates. \ continues line.
  std::string inputLine;
  unsigned int idx = 0;
  char lastchar = '0';
  char ptr = 0;
  Mode cmode = C_OK;
  while ( ptr != EOF ) {
    ptr = (char)fgetc(infile);
    // Skip leading whitespace
    if (idx == 0 && isspace(ptr)) {
      while ( (ptr = (char)fgetc(infile))!=EOF )
        if ( !isspace(ptr) ) break;
    }
    // If '#' is encountered, skip the rest of the line
    if (ptr=='#') 
      while (!EndChar(ptr)) ptr=(char)fgetc(infile);
    // newline, NULL, or EOF terminates the line
    if (EndChar(ptr)) {
      // If no chars in string continue
      if (inputLine.empty()) continue;
      // Print the input line that will be sent to dispatch
      mprintf("  [%s]\n",inputLine.c_str());
      // Call Dispatch to convert input to arglist and process.
      cmode = Dispatch(inputLine.c_str());
      if (cmode != C_OK) break;
      // Reset Input line
      inputLine.clear();
      idx = 0;
      continue;
    }
    // Any consecutive whitespace is skipped
    if (idx > 0) lastchar = inputLine[idx-1];
    if (isspace(ptr) && isspace(lastchar)) continue;
    // Backslash followed by newline continues to next line. Otherwise backslash
    // followed by next char will be inserted. 
    if (ptr=='\\') {
      ptr = (char)fgetc(infile);
      if ( ptr == EOF ) break;
      if (ptr == '\n' || ptr == '\r') continue;
      inputLine += "\\";
      inputLine += ptr;
      idx += 2;
      continue;
    }
    // Add character to input line
    inputLine += ptr;
    ++idx;
  }
  fclose(infile);
  if (cmode == C_ERR) return 1;
  return 0;
} 

/** Read command line args. */
Cpptraj::Mode Cpptraj::ProcessCmdLineArgs(int argc, char** argv) {
  if (argc == 1) return C_INTERACTIVE;
  bool hasInput = false;
  bool interactive = false;
  for (int i = 1; i < argc; i++) {
    std::string arg(argv[i]); 
    if ( arg == "--help" || arg == "-help" ) {
      // --help, -help: Print usage and exit
      Usage( argv[0] );
      return C_QUIT;
    }
    if ( arg == "-V" || arg == "--version" ) 
      // -V, --version: Print version number and exit
      // Since version number should be printed before this is called, quit.
      return C_QUIT;
    if ( arg == "--defines" ) {
      // --defines: Print information on compiler defines used and exit
      mprintf("\nCompiled with:");
#     ifdef DEBUG
      mprintf(" -DDEBUG");
#     endif
#     ifdef HASBZ2
      mprintf(" -DHASBZ2");
#     endif
#     ifdef HASGZ
      mprintf(" -DHASGZ");
#     endif
#     ifdef BINTRAJ
      mprintf(" -DBINTRAJ");
#     endif
#     ifdef MPI
      mprintf(" -DMPI");
#     endif
#     ifdef _OPENMP
      mprintf(" -D_OPENMP");
#     endif
#     ifdef NO_MATHLIB
      mprintf(" -DNO_MATHLIB");
#     endif
      mprintf("\n");
      return C_QUIT;
    }
    if ( arg == "--interactive" )
      interactive = true;
    else if ( arg == "-debug" && i+1 != argc) 
      // -debug: Set overall debug level
      SetGlobalDebug( convertToInteger( argv[++i] ) ); 
    else if ( arg == "-p" && i+1 != argc) {
      // -p: Topology file
      if (parmFileList.AddParmFile( argv[++i] )) return C_ERR;
    } else if (arg == "-i" && i+1 != argc) {
      // -i: Input file(s)
      if (ProcessInput( argv[++i] )) return C_ERR;
      hasInput = true;
    } else if ( i == 1 ) {
      // For backwards compatibility with PTRAJ; Position 1 = TOP file
      if (parmFileList.AddParmFile( argv[i])) return C_ERR;
    } else if ( i == 2 ) {
      // For backwards compatibility with PTRAJ; Position 2 = INPUT file
      if (ProcessInput( argv[i])) return C_ERR;
      hasInput = true;
    } else {
      // Unrecognized
      mprintf("  Unrecognized input on command line: %i: %s\n", i,argv[i]);
      Usage(argv[0]);
      return C_QUIT;
    }
  }
  if (!hasInput || interactive) return C_INTERACTIVE;
  // If Run has already been called, just quit.
  if (nrun_ > 0) return C_QUIT;
  return C_OK;
}

// Cpptraj::Dispatch()
/** The input line is converted into a whitespace-delimited array of
  * arguments, the first of which is considered the command. This command
  * is searched for and if it is recognized it is sent to the appropriate
  * class. 
  * \param inputLine NULL-terminated string consisting of command and arguments.
  * \return true if command was accepted or no error occurred.
  * \return false if error occurred or exit requested.
  */
// NOTE: Should differentiate between keyword rejection and outright error.
Cpptraj::Mode Cpptraj::Dispatch(const char* inputLine) {
  //mprintf("\t[%s]\n", inputLine);
  ArgList command( inputLine );
  command.MarkArg(0); // Always mark the command
  if ( SearchToken( command ) ) {
    //mprintf("TOKEN FOUND. CMD=%s  TYPE=%i\n", dispatchToken_->Cmd, (int)dispatchToken_->Type);
    switch (dispatchToken_->Type) {
      case DispatchObject::PARM :
        if ( parmFileList.CheckCommand(dispatchToken_->Idx, command) 
             && exitOnError_ )
          return C_ERR;
        break;
      case DispatchObject::COORD :
        switch ( dispatchToken_->Idx ) {
          case TRAJIN :
            if (trajinList.AddTrajin(command, parmFileList) && exitOnError_)
              return C_ERR;
            break;
          case TRAJOUT :
            // For setting up ensemble, save trajout arg
            trajoutArgs_.push_back(command);
            if (trajoutList.AddTrajout(command, parmFileList) && exitOnError_)
              return C_ERR;
            break;
          case REFERENCE :
            if (refFrames.AddReference(command, parmFileList) && exitOnError_)
              return C_ERR;
            break;
        }
        break;
      case DispatchObject::ACTION : 
        // For setting up ensemble, save action arg
        actionArgs_.push_back(command);
        if (actionList.AddAction( dispatchToken_->Alloc, command, &parmFileList,
                                  &refFrames, &DSL, &DFL ) != 0 && exitOnError_ )
          return C_ERR;
        break;
      case DispatchObject::ANALYSIS :
        if ( analysisList.AddAnalysis( dispatchToken_->Alloc, command, &parmFileList, &DSL ) 
             && exitOnError_)
          return C_ERR;
        break;
      case DispatchObject::GENERAL :
        switch ( dispatchToken_->Idx ) {
          case LIST  : List(command); break;
          case HELP  : Help(command); break;
          case DEBUG : Debug(command); break;
          case SELECTDS: SelectDS(command); break;
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
            if (ReadData( command ) && exitOnError_) return C_ERR;
            break;
          case READINPUT:
            if (ProcessInput( command.GetStringNext() ) && exitOnError_) return C_ERR;
            break;
          case CREATE:
            if (Create_DataFile( command ) && exitOnError_) return C_ERR;
            break;
          case PRECISION:
            if (Precision( command ) && exitOnError_) return C_ERR;
            break;
          case DATAFILE:
            if (DFL.ProcessDataFileArgs( command ) && exitOnError_) return C_ERR;
            break;
          case SYSTEM      : system( command.ArgLine() ); break;
          case RUN         : Run(); break;
          case RUN_ANALYSIS: analysisList.DoAnalyses(&DFL); break;
          case WRITEDATA   : if (worldrank == 0) DFL.Write(); break;
          case QUIT        : return C_QUIT; break;
        }
        break;
      default: mprintf("Dispatch type is currently not handled.\n");
    }
  }
  return C_OK;
}

// -----------------------------------------------------------------------------
// Cpptraj::Run()
int Cpptraj::Run() {
  int err = 0;
  switch ( trajinList.Mode() ) {
    case TrajinList::NORMAL   : err = RunNormal(); break;
    case TrajinList::ENSEMBLE : err = RunEnsemble(); break;
    default: err = 1; break;
  }
  return err;
}

// Cpptraj::RunEnsemble()
int Cpptraj::RunEnsemble() {
  FrameArray FrameEnsemble;

  mprintf("\nINPUT ENSEMBLE:\n");
  // Ensure all ensembles are of the same size
  int ensembleSize = -1;
  for (TrajinList::const_iterator traj = trajinList.begin(); traj != trajinList.end(); ++traj) 
  {
    Trajin_Multi* mtraj = (Trajin_Multi*)*traj;
    if (ensembleSize == -1)
      ensembleSize = mtraj->EnsembleSize();
    else if (ensembleSize != mtraj->EnsembleSize()) {
      mprinterr("Error: Ensemble size (%i) does not match first ensemble size (%i).\n",
                mtraj->EnsembleSize(), ensembleSize);
      return 1;
    }
    // Perform ensemble setup - this also resizes FrameEnsemble
    if ( mtraj->EnsembleSetup( FrameEnsemble ) ) return 1;
  }
  mprintf("  Ensemble size is %i\n", ensembleSize);

  // Calculate frame division among trajectories
  int maxFrames = trajinList.MaxFrames();
  // Parameter file information
  parmFileList.List();
  // Print reference information 
  mprintf("\nREFERENCE COORDS:\n");
  refFrames.List();

  // Allocate an ActionList, TrajoutList, and DataSetList for each
  // member of the ensemble.
  std::vector<ActionList> ActionEnsemble( ensembleSize );
  std::vector<TrajoutList> TrajoutEnsemble( ensembleSize );
  std::vector<DataSetList> DataSetEnsemble( ensembleSize );

  // Set up output trajectories for each member of the ensemble
  for (ArgsArray::iterator targ = trajoutArgs_.begin(); targ != trajoutArgs_.end(); ++targ)
  {
    for (int member = 0; member < ensembleSize; ++member) 
      TrajoutEnsemble[member].AddEnsembleTrajout( *targ, parmFileList, member );
  }
  mprintf("\n");
  for (int member = 0; member < ensembleSize; ++member) {
    mprintf("OUTPUT TRAJECTORIES Member %i:\n", member);
    TrajoutEnsemble[member].List();
  }

  // TODO: One loop over member?
  for (int member = 0; member < ensembleSize; ++member) {
    mprintf("***** ENSEMBLE MEMBER %i: ", member);
    // Set max frames in the data set list
    DataSetEnsemble[member].AllocateSets(maxFrames);
    // Initialize actions 
    for (ArgsArray::iterator aarg = actionArgs_.begin(); aarg != actionArgs_.end(); ++aarg)
    {
      if ( SearchToken( *aarg ) ) {
        // Create copy of arg list so that args remain unmarked for next member
        ArgList command = *aarg;
        if (ActionEnsemble[member].AddAction( dispatchToken_->Alloc, command, &parmFileList,
                                          &refFrames, &(DataSetEnsemble[member]), &DFL ))
          return 1;
      }
    }
  }
      
  // ========== A C T I O N  P H A S E ==========
  int lastPindex=-1;          // Index of the last loaded parm file
  int readSets = 0;
  int actionSet = 0;
  bool hasVelocity = false;
  // Loop over every trajectory in trajFileList
  rprintf("BEGIN ENSEMBLE PROCESSING:\n");
  for ( TrajinList::const_iterator traj = trajinList.begin();
                                   traj != trajinList.end(); ++traj)
  {
    // Open up the trajectory file. If an error occurs, bail 
    if ( (*traj)->BeginTraj(showProgress_) ) {
      mprinterr("Error: Could not open trajectory %s.\n",(*traj)->FullTrajStr());
      break;
    }
    // Set current parm from current traj.
    Topology* CurrentParm = (*traj)->TrajParm();
    // Check if parm has changed
    bool parmHasChanged = (lastPindex != CurrentParm->Pindex());

    // If Parm has changed or trajectory velocity status has changed,
    // reset the frame.
    if (parmHasChanged || (hasVelocity != (*traj)->HasVelocity()))
      FrameEnsemble.SetupFrames(CurrentParm->Atoms(), (*traj)->HasVelocity());
    hasVelocity = (*traj)->HasVelocity();

    // If Parm has changed, reset actions for new topology.
    if (parmHasChanged) {
      // Set active reference for this parm
      CurrentParm->SetReferenceCoords( refFrames.ActiveReference() );
      // Set up actions for this parm
      bool setupOK = true;
      for (int member = 0; member < ensembleSize; ++member) {
        if (ActionEnsemble[member].SetupActions( &CurrentParm )) {
          mprintf("Warning: Ensemble member %i: Could not set up actions for %s: skipping.\n",
                  member, CurrentParm->c_str());
          setupOK = false;
        }
      }
      if (!setupOK) continue;
      //fprintf(stdout,"DEBUG: After setup of Actions in Cpptraj parm name is %s\n",
      //        CurrentParm->parmName);
      lastPindex = CurrentParm->Pindex();
    }

    // Loop over every collection of frames in the ensemble
    (*traj)->PrintInfoLine();
    Trajin_Multi* mtraj = (Trajin_Multi*)*traj;
    while ( mtraj->GetNextEnsemble(FrameEnsemble) ) {
      // Loop over all members of the ensemble
      for (int member = 0; member < ensembleSize; ++member) {
        // Get this members current position
        int pos = mtraj->EnsemblePosition( member );
        // Since Frame can be modified by actions, save original and use CurrentFrame
        Frame* CurrentFrame = &(FrameEnsemble[member]);
        // Perform Actions on Frame
        bool suppress_output = ActionEnsemble[pos].DoActions(&CurrentFrame, actionSet);
        // Do Output
        if (!suppress_output) 
          TrajoutEnsemble[pos].Write(actionSet, CurrentParm, CurrentFrame);
      } // END loop over ensemble
      // Increment frame counter
      ++actionSet;
    }

    // Close the trajectory file
    (*traj)->EndTraj();
    // Update how many frames have been processed.
    readSets += (*traj)->NumFramesProcessed();
    mprintf("\n");
  } // End loop over trajin
  rprintf("Read %i frames and processed %i frames.\n",readSets,actionSet);

  // Close output trajectories
  for (int member = 0; member < ensembleSize; ++member)
    TrajoutEnsemble[member].Close();

  // ========== A C T I O N  O U T P U T  P H A S E ==========
  for (int member = 0; member < ensembleSize; ++member)
    actionList.Print( );

  // Sync DataSets and print DataSet information
  // TODO - Also have datafilelist call a sync??
  for (int member = 0; member < ensembleSize; ++member) {
    DataSetEnsemble[member].Sync();
    DataSetEnsemble[member].sort();
    mprintf("\nENSEMBLE MEMBER %i DATASETS:\n",member);
    DataSetEnsemble[member].List();
  }

  // ========== A N A L Y S I S  P H A S E ==========
//  analysisList.Setup(&DSL, &parmFileList);
//  analysisList.Analyze(&DFL);

  // DEBUG: DataSets, post-Analysis
//  mprintf("\nDATASETS AFTER ANALYSIS:\n");
//  DSL.Info();

  // ========== D A T A  W R I T E  P H A S E ==========
  // Process any datafile commands
//  DFL.ProcessDataFileArgs(&DSL);
  // Print Datafile information
  DFL.List();
  // Only Master does DataFile output
  if (worldrank==0)
    DFL.Write();

  return 0;
}

// Cpptraj::RunNormal()
/** Process trajectories in trajinList. Each frame in trajinList is sent
 *  to the actions in actionList for processing.
 */
int Cpptraj::RunNormal() {
  int actionSet=0;            // Internal data frame
  int readSets=0;             // Number of frames actually read
  int lastPindex=-1;          // Index of the last loaded parm file
  Frame TrajFrame;            // Original Frame read in from traj
  ++nrun_;

  // ========== S E T U P   P H A S E ========== 
  // Parameter file information
  parmFileList.List();
  // Input coordinate file information
  trajinList.List();
  // Print reference information 
  mprintf("\nREFERENCE COORDS:\n");
  refFrames.List();
  // Output traj
  mprintf("\nOUTPUT TRAJECTORIES:\n");
  trajoutList.List();
  // Allocate DataSets in the master DataSetList based on # frames to be read
  DSL.AllocateSets(trajinList.MaxFrames()); 
  
  // ========== A C T I O N  P H A S E ==========
  // Loop over every trajectory in trajFileList
  rprintf("BEGIN TRAJECTORY PROCESSING:\n");
  for ( TrajinList::const_iterator traj = trajinList.begin();
                                   traj != trajinList.end(); ++traj)
  {
    // Open up the trajectory file. If an error occurs, bail 
    if ( (*traj)->BeginTraj(showProgress_) ) {
      mprinterr("Error: Could not open trajectory %s.\n",(*traj)->FullTrajStr());
      break;
    }
    // Set current parm from current traj.
    Topology* CurrentParm = (*traj)->TrajParm();
    // Check if parm has changed
    bool parmHasChanged = (lastPindex != CurrentParm->Pindex());

    // If Parm has changed or trajectory velocity status has changed,
    // reset the frame.
    if (parmHasChanged || (TrajFrame.HasVelocity() != (*traj)->HasVelocity()))
      TrajFrame.SetupFrameV(CurrentParm->Atoms(), (*traj)->HasVelocity());

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
    (*traj)->PrintInfoLine();
    while ( (*traj)->GetNextFrame(TrajFrame) ) {
      // Since Frame can be modified by actions, save original and use CurrentFrame
      Frame* CurrentFrame = &TrajFrame;
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
    (*traj)->EndTraj();
    // Update how many frames have been processed.
    readSets += (*traj)->NumFramesProcessed();
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
