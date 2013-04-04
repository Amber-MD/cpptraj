#include <cstdio> 
#include <cstdlib> // system
#include "Cpptraj.h"
#include "MpiRoutines.h"
#include "CpptrajStdio.h"
#include "Trajin_Multi.h"
#include "Trajin_Single.h"
#include "FrameArray.h"
#include "ReadLine.h"
#include "ParmFile.h"
#include "DataSet_Coords.h" // CrdAction
#include "Command.h"

void Cpptraj::Usage() {
  mprinterr("\nUsage: cpptraj [-p <Top0>] [-i <Input0>] [-y <trajin>] [-x <trajout>]\n");
  mprinterr("               [-h | --help] [-V | --version] [--defines] [-debug <#>]\n");
  mprinterr("               [--interactive] [--log <logfile>]\n");
  mprinterr("       cpptraj <Top> <Input>\n");
  mprinterr("\t-p <Top0>      : Load <Top0> as a topology file. May be specified more than once.\n");
  mprinterr("\t-i <Input0>    : Read input from <Input0>. May be specified more than once.\n");
  mprinterr("\t-y <trajin>    : Read from trajectory file <trajin>; same as input 'trajin <trajin>'.\n");
  mprinterr("\t-x <trajout>   : Write trajectory file <trajout>; same as input 'trajout <trajout>'.\n");
  mprinterr("\t-h | --help    : Print command line help and exit.\n");
  mprinterr("\t-V | --version : Print version and exit.\n");
  mprinterr("\t--defines      : Print compiler defines and exit.\n");
  mprinterr("\t-debug <#>     : Set global debug level to <#>; same as input 'debug <#>'.\n");
  mprinterr("\t--interactive  : Force interactive mode.\n");
  mprinterr("\t--log <logfile>: Record commands to <logfile> (interactive mode only). Default is 'cpptraj.log'.\n");
}

// -----------------------------------------------------------------------------
// Constructor
Cpptraj::Cpptraj() : 
  debug_(0),
  showProgress_(true),
  exitOnError_(true),
  nrun_(0)
{}

/** List all commands, or call help function of specific command. */
void Cpptraj::Help(ArgList& argIn) {
  ArgList arg = argIn;
  arg.RemoveFirstArg();
  if (arg.empty()) 
    // NONE in this context means list all commands
    Command::List(DispatchObject::NONE);
  else if (arg.CommandIs("General"))
    Command::List(DispatchObject::GENERAL);
  else if (arg.CommandIs("Topology"))
    Command::List(DispatchObject::PARM);
  else if (arg.CommandIs("Action"))
    Command::List(DispatchObject::ACTION);
  else if (arg.CommandIs("Analysis"))
    Command::List(DispatchObject::ANALYSIS);
  else if (arg.CommandIs("Trajectory"))
    Command::List(DispatchObject::TRAJ);
  else {
    DispatchObject::TokenPtr dispatchToken = Command::SearchToken( arg );
    if (dispatchToken == 0 || dispatchToken->Help == 0) 
      mprinterr("No help found for %s\n", arg.Command());
    else
      dispatchToken->Help();
  }
}

/// Types of lists
enum ListType { L_ACTION = 0, L_TRAJIN, L_REF, L_TRAJOUT, L_PARM, L_ANALYSIS,
                L_DATAFILE, L_DATASET, N_LISTS };
/// Select lists from ArgList
static std::vector<bool> ListsFromArg( ArgList& argIn, bool allowEnableAll ) {
  std::vector<bool> enabled( (int)N_LISTS );
  enabled[L_ACTION]   = argIn.hasKey("actions");
  enabled[L_TRAJIN]   = argIn.hasKey("trajin");
  enabled[L_REF]      = argIn.hasKey("ref");
  enabled[L_TRAJOUT]  = argIn.hasKey("trajout");
  enabled[L_PARM]     = argIn.hasKey("parm");
  enabled[L_ANALYSIS] = argIn.hasKey("analysis");
  enabled[L_DATAFILE] = argIn.hasKey("datafile");
  enabled[L_DATASET]  = argIn.hasKey("dataset");
  if (!allowEnableAll) return enabled;
  // If nothing is enabled, set all enabled
  bool nothing_enabled = true;
  for (std::vector<bool>::iterator en = enabled.begin(); en != enabled.end(); ++en)
    if (*en) {
      nothing_enabled = false;
      break;
    }
  if (nothing_enabled) enabled.assign( (int)N_LISTS, true );
  return enabled;
}
// Cpptraj::ListAction()
void Cpptraj::ListAction(ArgList& argIn, int cmdidx) {
  std::vector<bool> enabled;
  switch (cmdidx) {
    case Command::LIST: /** List all members of specified lists. */
      enabled = ListsFromArg( argIn, true );
      if ( enabled[L_ACTION]   ) actionList_.List();
      if ( enabled[L_TRAJIN]   ) trajinList_.List();
      if ( enabled[L_REF]      ) {mprintf("\nREFERENCE COORDS:\n");refFrames_.List();}
      if ( enabled[L_TRAJOUT]  ) {mprintf("\nOUTPUT TRAJECTORIES:\n");trajoutList_.List();}
      if ( enabled[L_PARM]     ) parmFileList_.List();
      if ( enabled[L_ANALYSIS] ) analysisList_.List();
      if ( enabled[L_DATAFILE] ) DFL_.List();
      if ( enabled[L_DATASET]  ) {mprintf("\nDATASETS:\n");DSL_.List();}
      break;
    case Command::DEBUG: /** Set debug level of specified lists */
      enabled = ListsFromArg( argIn, true );
      debug_ = argIn.getNextInteger(0);
      if ( enabled[L_ACTION]   ) actionList_.SetDebug( debug_ );
      if ( enabled[L_TRAJIN]   ) trajinList_.SetDebug( debug_ );
      if ( enabled[L_REF]      ) refFrames_.SetDebug( debug_ );
      if ( enabled[L_TRAJOUT]  ) trajoutList_.SetDebug( debug_ );
      if ( enabled[L_PARM]     ) parmFileList_.SetDebug( debug_ );
      if ( enabled[L_ANALYSIS] ) analysisList_.SetDebug( debug_ );
      if ( enabled[L_DATAFILE] ) DFL_.SetDebug( debug_ );
      if ( enabled[L_DATASET]  ) DSL_.SetDebug( debug_ );
      break;
    case Command::CLEAR: /** Clear specified lists */
      enabled = ListsFromArg( argIn, argIn.hasKey("all") );
      if ( enabled[L_ACTION]   ) actionList_.Clear();
      if ( enabled[L_TRAJIN]   ) trajinList_.Clear();
      if ( enabled[L_REF]      ) refFrames_.Clear();
      if ( enabled[L_TRAJOUT]  ) trajoutList_.Clear();
      if ( enabled[L_PARM]     ) parmFileList_.Clear();
      if ( enabled[L_ANALYSIS] ) analysisList_.Clear();
      if ( enabled[L_DATAFILE] ) DFL_.Clear();
      if ( enabled[L_DATASET]  ) DSL_.Clear();
      break;
  }
}

/** Add DataFile to DataFileList using specified sets. */
int Cpptraj::Create_DataFile(ArgList& dataArg, int cmdidxIn) {
  // Two modes:
  //   1) create: Add a new DataFile to DFL with specified DataSets.
  //   2) write:  Immediately write DataFile with specified DataSets.
  Command::Type cmdidx = (Command::Type)cmdidxIn;
  // Next string is datafile that command pertains to.
  std::string name1 = dataArg.GetStringNext();
  if (name1.empty()) {
    mprinterr("Error: No filename given.\n");
    return 1;
  }
  DataFile* df = 0;
  if ( cmdidx == Command::CREATE )
    df = DFL_.AddDataFile(name1, dataArg);
  else { // WRITE
    df = new DataFile();
    if (df->SetupDatafile( name1, dataArg, debug_ )) {
      delete df;
      df = 0;
    }
  }
  if (df==0) {
    mprinterr("Error: Could not create file %s:",name1.c_str());
    return 1;
  }
  // Treat all remaining args as dataset names
  int err = 0;
  ArgList dsetArgs = dataArg.RemainingArgs();
  for (ArgList::const_iterator dsa = dsetArgs.begin(); dsa != dsetArgs.end(); ++dsa) {
    DataSetList Sets = DSL_.GetMultipleSets( *dsa );
    if (Sets.empty())
      mprintf("Warning: %s does not correspond to any data sets.\n", (*dsa).c_str());
    for (DataSetList::const_iterator set = Sets.begin(); set != Sets.end(); ++set) {
      mprintf(" %s", (*set)->Legend().c_str());
      if ( df->AddSet(*set) ) {
        mprinterr("Error: Could not add data set %s to file.\n", (*set)->Legend().c_str());
        ++err;
      }
    }
  }
  mprintf("\n");
  if ( cmdidx == Command::WRITE && err == 0 ) {
    df->Write();
    delete df;
  }
  return err;
}

/** Set precision for specific set or all sets in specified DataFile */
int Cpptraj::Precision(ArgList& dataArg) {
  // Next string is DataSet(s)/DataFile that command pertains to.
  std::string name1 = dataArg.GetStringNext();
  if (name1.empty()) {
    mprinterr("Error: precision: No filename/setname given.\n");
    return 1;
  }
  // This will break if dataset name starts with a digit...
  int width = dataArg.getNextInteger(12);
  if (width < 1) {
    mprintf("Error: precision: Cannot set width < 1 (%i).\n", width);
    return 1;
  }
  int precision = dataArg.getNextInteger(4);
  if (precision < 0) precision = 0;
  DataFile* df = DFL_.GetDataFile(name1);
  if (df != 0) {
    mprintf("\tSetting precision for all sets in %s to %i.%i\n", df->DataFilename().base(),
            width, precision);
    df->SetPrecision(width, precision);
  } else {
    DataSetList dsets = DSL_.GetMultipleSets( name1 );
    mprintf("\tSetting precision for %i sets to %i.%i\n", dsets.size(),
            width, precision);
    for (DataSetList::const_iterator set = dsets.begin(); set != dsets.end(); ++set)
      (*set)->SetPrecision(width, precision);
  }
  return 0;
}

/** Read data from file into master DataSetList */
int Cpptraj::ReadData(ArgList& argIn) {
  DataFile dataIn;
  if (dataIn.ReadData( argIn, DSL_ )!=0) {
    mprinterr("Error: Could not read data file.\n");
    return 1;
  }
  return 0;
}

/** Show results of DataSet selection */
void Cpptraj::SelectDS(ArgList& argIn) {
  std::string dsarg = argIn.GetStringNext();
  DataSetList dsets = DSL_.GetMultipleSets( dsarg );
  mprintf("SelectDS: Arg [%s]:", dsarg.c_str());
  dsets.List();
}

// -----------------------------------------------------------------------------
/** Load file into TopologyList */
int Cpptraj::LoadParm(ArgList& argIn) {
  std::string parmtag = argIn.getNextTag();
  bool bondsearch = !argIn.hasKey("nobondsearch");
  double offset = argIn.getKeyDouble("bondsearch", -1.0);
  return parmFileList_.AddParmFile(argIn.GetStringNext(), parmtag, bondsearch, offset);
}

/** Print information for specified parm */
int Cpptraj::ParmInfo(ArgList& argIn, int cmdidxIn) {
  Command::Type cmdidx = (Command::Type) cmdidxIn;
  int pindex = argIn.getNextInteger(0);
  Topology* parm = parmFileList_.GetParm( pindex );
  if ( parm == 0 ) {
    mprinterr("Error: parm index %i not loaded.\n",pindex);
    return 1;
  }
  std::string maskarg = argIn.GetMaskNext();
  switch (cmdidx) {
    case Command::PARMINFO:
      if (!maskarg.empty())
        parm->PrintAtomInfo( maskarg );
      else
        parm->Summary();
      break;
    case Command::BONDINFO  : parm->PrintBondInfo( maskarg ); break;
    case Command::RESINFO   : parm->PrintResidueInfo( maskarg ); break;
    case Command::MOLINFO   : parm->PrintMoleculeInfo( maskarg ); break;
    case Command::CHARGEINFO: parm->PrintChargeInfo( maskarg ); break;
    default: return 1; // Should never get here
  }
  return 0;
}

/** Write parm to Amber Topology file. */
int Cpptraj::ParmWrite(ArgList& argIn) {
  std::string outfilename = argIn.GetStringKey("out");
  if (outfilename.empty()) {
    mprinterr("Error: parmwrite: No output filename specified (use 'out <filename>').\n");
    return 1;
  }
  int pindex = argIn.getNextInteger(0);
  Topology* parm = parmFileList_.GetParm( pindex );
  if ( parm == 0 ) {
    mprinterr("Error: parmwrite: parm index %i not loaded.\n",pindex);
    return 1;
  }
  mprintf("\tWriting parm %i (%s) to Amber parm %s\n",pindex,
          parm->c_str(), outfilename.c_str());
  ParmFile pfile;
  pfile.Write( *parm, outfilename, ParmFile::AMBERPARM, debug_ );
  return 0;
}

// Cpptraj::ParmStrip()
int Cpptraj::ParmStrip(ArgList& argIn) {
  int pindex = argIn.getNextInteger(0);
  Topology* parm = parmFileList_.GetParm( pindex );
  if ( parm == 0 ) {
    mprinterr("Error: parmstrip: parm index %i not loaded.\n",pindex);
    return 1;
  }
  AtomMask tempMask( argIn.GetMaskNext() );
  // Since want to keep atoms outside mask, invert selection
  tempMask.InvertMask();
  parm->SetupIntegerMask( tempMask );
  mprintf("\tStripping atoms in mask [%s] (%i) from %s\n",tempMask.MaskString(),
           parm->Natom() - tempMask.Nselected(), parm->c_str());
  Topology* tempParm = parm->modifyStateByMask(tempMask);
  if (tempParm==0) {
    mprinterr("Error: parmstrip: Could not strip parm.\n");
    return 1;
  } else {
    // Replace parm with stripped version
    tempParm->ParmInfo();
    mprintf("\n");
    parmFileList_.ReplaceParm(pindex, tempParm);
  }
  return 0;
}

/** Modify parm box information. */
int Cpptraj::ParmBox(ArgList& argIn) {
  int pindex = argIn.getNextInteger(0);
  Topology* parm = parmFileList_.GetParm( pindex );
  if ( parm == 0 ) {
    mprinterr("Error: parmbox: parm index %i not loaded.\n",pindex);
    return 1;
  }
  if ( argIn.hasKey("nobox") )
    parm->SetBox( Box() );
  else {
    Box pbox;
    pbox.SetX( argIn.getKeyDouble("x",0) );
    pbox.SetY( argIn.getKeyDouble("y",0) );
    pbox.SetZ( argIn.getKeyDouble("z",0) );
    pbox.SetAlpha( argIn.getKeyDouble("alpha",0) );
    pbox.SetBeta(  argIn.getKeyDouble("beta",0)  );
    pbox.SetGamma( argIn.getKeyDouble("gamma",0) );
    // Fill in missing parm box information from specified parm
    pbox.SetMissingInfo( parm->ParmBox() );
    parm->SetBox( pbox );
  }
  return 0;
}

/** Modify parm solvent information */
int Cpptraj::ParmSolvent(ArgList& argIn) {
  std::string maskexpr = argIn.GetMaskNext();
  if ( maskexpr.empty() ) {
    mprinterr("Error: solvent: No mask specified.\n");
    return 1;
  }
  // Get parm index
  int pindex = argIn.getNextInteger(0);
  Topology* parm = parmFileList_.GetParm( pindex );
  if ( parm == 0 ) {
    mprinterr("Error: solvent: parm index %i not loaded.\n",pindex);
    return 1;
  } 
  parm->SetSolvent( maskexpr );
  return 0;
}

/** Show results of mask expression */
int Cpptraj::Select(ArgList& argIn) {
  AtomMask tempMask( argIn.GetMaskNext() );
  int pindex = argIn.getNextInteger(0);
  Topology* parm = parmFileList_.GetParm( pindex );
  if ( parm == 0 ) {
    mprinterr("Error: solvent: parm index %i not loaded.\n",pindex);
    return 1;
  }
  parm->SetupIntegerMask( tempMask );
  mprintf("Selected %i atoms.\n", tempMask.Nselected());
  if (!argIn.hasKey("total"))
    tempMask.PrintMaskAtoms("Selected");
  return 0;
}

// -----------------------------------------------------------------------------
// Cpptraj::LoadCrd()
int Cpptraj::LoadCrd(ArgList& argIn) {
  // Get parm
  Topology* parm = parmFileList_.GetParm( argIn );
  if (parm == 0) {
    mprinterr("Error: loadcrd: No parm files loaded.\n");
    return 1;
  }
  // Load trajectory
  Trajin_Single trajin;
  trajin.SetDebug( debug_ );
  if (trajin.SetupTrajRead(argIn.GetStringNext(), &argIn, parm)) {
    mprinterr("Error: loadcrd: Could not set up input trajectory.\n");
    return 1;
  }
  // Create input frame
  Frame frameIn;
  frameIn.SetupFrameV(parm->Atoms(), trajin.HasVelocity());
  // Create DataSet, use base file name as set name if none specified. 
  // NOTE: Default name should NEVER get used.
  std::string setname = argIn.GetStringNext();
  if (setname.empty())
    setname = trajin.TrajFilename().Base();
  DataSet_Coords* coords = (DataSet_Coords*)DSL_.AddSet(DataSet::COORDS, setname, "__DCRD__");
  if (coords == 0) {
    mprinterr("Error: loadcrd: Could not set up COORDS data set.\n");
    return 1;
  }
  coords->SetTopology( *parm );
  // Read trajectory
  mprintf("\tLoading trajectory %s as \"%s\"\n", trajin.TrajFilename().full(), setname.c_str());
  trajin.BeginTraj(true);
  trajin.PrintInfoLine();
  while (trajin.GetNextFrame( frameIn ))
    coords->AddFrame( frameIn );
  trajin.EndTraj();
  return 0;
}
  
// Cpptraj::CrdAction()
/** Perform action on given COORDS dataset */
int Cpptraj::CrdAction(ArgList& argIn) {
  std::string setname = argIn.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: crdaction: Specify COORDS dataset name.\n");
    return 1;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)DSL_.FindSetOfType( setname, DataSet::COORDS );
  if (CRD == 0) {
    mprinterr("Error: crdaction: No COORDS set with name %s found.\n", setname.c_str());
    return 1;
  }
  // Start, stop, offset
  ArgList crdarg( argIn.GetStringKey("crdframes"), "," );
  int start = crdarg.getNextInteger(1) - 1;
  int stop;
  if (crdarg.hasKey("last"))
    stop = CRD->Size();
  else
    stop = crdarg.getNextInteger(CRD->Size());
  int offset = crdarg.getNextInteger(1);
  if (debug_ > 0) mprintf("\tDBG: Frames %i to %i, offset %i\n", start+1, stop, offset);
  ArgList actionargs = argIn.RemainingArgs();
  actionargs.MarkArg(0);
  DispatchObject::TokenPtr tkn = Command::SearchTokenType( DispatchObject::ACTION, actionargs);
  if ( tkn == 0 ) return 1;
  Action* act = (Action*)tkn->Alloc();
  if (act == 0) return 1;
  DataFileList dfl;
  if ( act->Init( actionargs, &parmFileList_, &refFrames_, &DSL_, &dfl, debug_ ) != Action::OK ) {
    delete act;
    return 1;
  }
  actionargs.CheckForMoreArgs();
  // Set up frame and parm for COORDS.
  Topology* originalParm = new Topology();
  *originalParm = CRD->Top();
  Frame* originalFrame = new Frame( CRD->Top().Atoms() );
  // Set up for this topology
  Topology* currentParm = originalParm;
  if ( act->Setup( currentParm, &currentParm ) == Action::ERR ) {
    delete act;
    return 1;
  }
  // Check if parm was modified. If so, update COORDS.
  if ( currentParm != originalParm ) {
    mprintf("Info: crdaction: Parm for %s was modified by action %s\n", 
            CRD->Legend().c_str(), actionargs.Command());
    CRD->SetTopology( *currentParm );
  }
  // Loop over all frames in COORDS.
  ProgressBar progress( stop );
  for (int frame = start; frame < stop; frame += offset) {
    progress.Update( frame );
    CRD->GetFrame( frame, *originalFrame );
    Frame* currentFrame = originalFrame;
    if (act->DoAction( frame, currentFrame, &currentFrame ) == Action::ERR) {
      mprinterr("Error: crdaction: Frame %i\n", frame + 1);
      break;
    }
    // Check if frame was modified. If so, update COORDS.
    // TODO: Have actions indicate whether they will modify coords
    //if ( currentFrame != originalFrame ) 
      CRD->SetCRD( frame, *currentFrame );
  }
  act->Print();
  if (worldrank == 0) dfl.Write();
  delete originalFrame;
  delete originalParm;
  delete act;
  return 0;
}

// Cpptraj::CrdOut()
/** Write out COORDS dataset */
int Cpptraj::CrdOut(ArgList& argIn) {
  std::string setname = argIn.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: crdout: Specify COORDS dataset name.\n");
    return 1;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)DSL_.FindSetOfType( setname, DataSet::COORDS );
  if (CRD == 0) {
    mprinterr("Error: crdout: No COORDS set with name %s found.\n", setname.c_str());
    return 1;
  }
  setname = argIn.GetStringNext();
  // Start, stop, offset
  ArgList crdarg( argIn.GetStringKey("crdframes"), "," );
  int start = crdarg.getNextInteger(1) - 1;
  int stop;
  if (crdarg.hasKey("last"))
    stop = CRD->Size();
  else
    stop = crdarg.getNextInteger(CRD->Size());
  int offset = crdarg.getNextInteger(1);
  if (debug_ > 0) mprintf("\tDBG: Frames %i to %i, offset %i\n", start+1, stop, offset);
  Trajout outtraj;
  Topology* currentParm = (Topology*)&(CRD->Top()); // TODO: Fix cast
  if (outtraj.SetupTrajWrite( setname, &argIn, currentParm, TrajectoryFile::UNKNOWN_TRAJ)) {
    mprinterr("Error: crdout: Could not set up output trajectory.\n");
    return 1;
  }
  outtraj.PrintInfo( 1 );
  Frame currentFrame( CRD->Top().Natom() );
  ProgressBar progress( stop );
  for (int frame = start; frame < stop; frame += offset) {
    progress.Update( frame );
    CRD->GetFrame( frame, currentFrame );
    if ( outtraj.WriteFrame( frame, currentParm, currentFrame ) ) {
      mprinterr("Error writing %s to output trajectory, frame %i.\n", 
                CRD->Legend().c_str(), frame + 1);
      break;
    }
  }
  return 0;
}

// Cpptraj::CrdAnalyze()
/** Run a single analysis. */
int Cpptraj::CrdAnalyze(ArgList& argIn) {
  ArgList analyzeargs = argIn.RemainingArgs();
  analyzeargs.MarkArg(0);
  DispatchObject::TokenPtr tkn = Command::SearchTokenType( DispatchObject::ANALYSIS, analyzeargs);
  if ( tkn == 0 ) return 1;
  Analysis* ana = (Analysis*)tkn->Alloc();
  if (ana == 0) return 1;
  DataFileList dfl;
  if ( ana->Setup( analyzeargs, &DSL_, &parmFileList_, &dfl, debug_ ) != Analysis::OK ) {
    delete ana;
    return 1;
  }
  int err = 0;
  if (ana->Analyze() == Analysis::ERR) 
    err = 1;
  else if (worldrank == 0) dfl.Write();
  delete ana;
  return err;
}

// -----------------------------------------------------------------------------
// Cpptraj::Interactive()
Cpptraj::Mode Cpptraj::Interactive() {
  ReadLine inputLine;
  // By default when interactive do not exit on errors
  exitOnError_ = false;
  // Open log file. If no name has been set, use default.
  if (!logfile_.Filename().empty())
    logfile_.OpenFile();
  else
    logfile_.OpenAppend("cpptraj.log");
  Mode readLoop = C_OK;
  while ( readLoop == C_OK ) {
    if (inputLine.GetInput()) break; 
    if (!inputLine.empty()) {
      readLoop = Dispatch( *inputLine );
      if (logfile_.IsOpen())
        logfile_.Printf("%s\n", inputLine.c_str());
    }
  }
  logfile_.CloseFile();
  // If we broke out of loop because of EOF and Run has been called at least
  // once, indicate that we can safely quit.
  if (readLoop == C_OK && nrun_ > 0) return C_QUIT;
  return readLoop;
}

static inline bool EndChar(char ptr) {
  if (ptr=='\n' || ptr=='\r' || ptr=='\0' || ptr==EOF) return true;
  return false;
}

/** Read commands from an input file, or from STDIN if given filename
  * is empty. '#' indicates the beginning of a comment, backslash at the 
  * end of a line indicates continuation (otherwise indicates 'literal').
  * \return 0 if successfully read, 1 on error.
  */
Cpptraj::Mode Cpptraj::ProcessInput(std::string const& inputFilename) {
  FILE* infile;
  if (inputFilename.empty()) {
    mprintf("INPUT: Reading Input from STDIN\n");
    infile = stdin;
  } else {
    mprintf("INPUT: Reading Input from file %s\n",inputFilename.c_str());
    if ( (infile=fopen(inputFilename.c_str(),"r"))==0 ) {
      rprintf("Error: Could not open input file %s\n",inputFilename.c_str());
      return C_ERR;
    }
  }
  // Read in each line of input. Newline or null terminates. \ continues line.
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
    // newline, null, or EOF terminates the line
    if (EndChar(ptr)) {
      // If no chars in string continue
      if (inputLine.empty()) continue;
      // Print the input line that will be sent to dispatch
      mprintf("  [%s]\n",inputLine.c_str());
      // Call Dispatch to convert input to arglist and process.
      cmode = Dispatch(inputLine);
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
  if (!inputFilename.empty())
    fclose(infile);
  // If everything OK and Run has already been called, just quit.
  if (cmode == C_OK && nrun_ > 0) return C_QUIT;
  return cmode;
} 

/** Read command line args. */
Cpptraj::Mode Cpptraj::ProcessCmdLineArgs(int argc, char** argv) {
  if (argc == 1) return C_INTERACTIVE;
  bool hasInput = false;
  bool interactive = false;
  for (int i = 1; i < argc; i++) {
    std::string arg(argv[i]); 
    if ( arg == "--help" || arg == "-h" ) {
      // --help, -help: Print usage and exit
      Usage();
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
    else if ( arg == "-debug" && i+1 != argc) { 
      // -debug: Set overall debug level
      ArgList dbgarg( argv[++i] );
      ListAction( dbgarg, Command::DEBUG );
    } else if ( arg == "--log" && i+1 != argc)
      // --log: Set up log file for interactive mode
      logfile_.SetupWrite( argv[++i], debug_ );
    else if ( arg == "-p" && i+1 != argc) {
      // -p: Topology file
      if (parmFileList_.AddParmFile( argv[++i] )) return C_ERR;
    } else if ( arg == "-y" && i+1 != argc) {
      // -y: Trajectory file in
      if (Dispatch("trajin " + std::string(argv[++i])) == C_ERR) return C_ERR;
    } else if ( arg == "-x" && i+1 != argc) {
      // -x: Trajectory file out 
      if (Dispatch("trajout " + std::string(argv[++i])) == C_ERR) return C_ERR;
      hasInput = true;
    } else if ( arg == "-c" && i+1 != argc) {
      // -c: Reference file
      if (Dispatch("reference " + std::string(argv[++i])) == C_ERR) return C_ERR;
    } else if (arg == "-i" && i+1 != argc) {
      // -i: Input file(s)
      Cpptraj::Mode cmode = ProcessInput( argv[++i] );
      if (cmode == C_ERR) return C_ERR;
      if (cmode == C_QUIT) return C_QUIT;
      hasInput = true;
    } else if (arg == "-ms" && i+1 != argc) {
      // -ms: Mask string
      ArgList maskArg( argv[++i] );
      ParmInfo( maskArg, Command::PARMINFO );
      return C_QUIT; 
    } else if ( i == 1 ) {
      // For backwards compatibility with PTRAJ; Position 1 = TOP file
      if (parmFileList_.AddParmFile( argv[i])) return C_ERR;
    } else if ( i == 2 ) {
      // For backwards compatibility with PTRAJ; Position 2 = INPUT file
      Cpptraj::Mode cmode = ProcessInput( argv[i] );
      if (cmode == C_ERR) return C_ERR;
      if (cmode == C_QUIT) return C_QUIT;
      hasInput = true;
    } else {
      // Unrecognized
      mprintf("  Unrecognized input on command line: %i: %s\n", i,argv[i]);
      Usage();
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
  * \param inputLine null-terminated string consisting of command and arguments.
  * \return C_OK if command was accepted or no error occurred.
  * \return C_ERR if error occurred.
  * \return C_QUIT if quit requested.
  */
Cpptraj::Mode Cpptraj::Dispatch(std::string const& inputLine) {
  int err = 0;
  //mprintf("\t[%s]\n", inputLine);
  ArgList command( inputLine );
  if ( command.empty() ) return C_OK;
  command.MarkArg(0); // Always mark the command
  DispatchObject::TokenPtr dispatchToken = Command::SearchToken( command );
  if ( dispatchToken != 0 ) {
    //mprintf("TOKEN FOUND. CMD=%s  TYPE=%i\n", dispatchToken->Cmd, (int)dispatchToken->Type);
    switch (dispatchToken->Type) {
      case DispatchObject::PARM : /** PARM COMMANDS */
        switch ( dispatchToken->Idx ) {
          case Command::LOADPARM : err = LoadParm( command ); break;
          case Command::BONDINFO :
          case Command::RESINFO  :
          case Command::MOLINFO  :
          case Command::CHARGEINFO:
          case Command::PARMINFO : err = ParmInfo( command, dispatchToken->Idx ); break;
          case Command::PARMWRITE: err = ParmWrite( command ); break;
          case Command::PARMSTRIP: err = ParmStrip( command ); break;
          case Command::PARMBOX  : err = ParmBox( command ); break;
          case Command::SOLVENT  : err = ParmSolvent(command); break;
        } 
        break;
      case DispatchObject::TRAJ : /** TRAJECTORY COMMANDS */
        switch ( dispatchToken->Idx ) {
          case Command::TRAJIN :
            // Update # of sets to be read in for master DSL
            err = trajinList_.AddTrajin(command, parmFileList_); 
            if (err == 0) DSL_.SetMax( trajinList_.MaxFrames() );
            break;
          case Command::TRAJOUT :
            // For setting up ensemble, save trajout arg
            trajoutArgs_.push_back(command);
            err = trajoutList_.AddTrajout(command, parmFileList_);
            break;
          case Command::REFERENCE : err = refFrames_.AddReference(command, parmFileList_); break;
        }
        break;
      case DispatchObject::ACTION : /** ACTION COMMANDS */ 
        // For setting up ensemble, save action arg
        actionArgs_.push_back(command);
        err = actionList_.AddAction( dispatchToken->Alloc, command, &parmFileList_, 
                                    &refFrames_, &DSL_, &DFL_ );
        break;
      case DispatchObject::ANALYSIS : /** ANALYSIS COMMANDS */
        err = analysisList_.AddAnalysis( dispatchToken->Alloc, command, &parmFileList_, &DSL_, &DFL_ );
        break;
      case DispatchObject::GENERAL : /** GENERAL COMMANDS */
        switch ( dispatchToken->Idx ) {
          case Command::HELP      : Help(command); break;
          case Command::LIST      : 
          case Command::DEBUG     : 
          case Command::CLEAR     : ListAction(command, dispatchToken->Idx); break;
          case Command::LOADCRD   : err = LoadCrd(command); break; 
          case Command::CRDACTION : err = CrdAction(command); break;
          case Command::CRDOUT    : err = CrdOut(command); break;
          case Command::SELECT    : err = Select(command); break; 
          case Command::SELECTDS  : SelectDS(command); break;
          case Command::NOPROG    : 
            showProgress_ = false; 
            mprintf("\tProgress bar will not be shown.\n");
            break;
          case Command::NOEXITERR:
            exitOnError_ = false;
            mprintf("\tcpptraj will attempt to ignore errors if possible.\n");
            break;
          case Command::ACTIVEREF : refFrames_.SetActiveRef( command.getNextInteger(0) ); break;
          case Command::READDATA  : err = ReadData( command ); break;
          case Command::READINPUT :
            switch (ProcessInput( command.GetStringNext() )) {
              case C_ERR  : if ( exitOnError_ ) return C_ERR; break;
              case C_QUIT : return C_QUIT; break;
              default     : break;
            } 
            break;
          case Command::WRITE       :
          case Command::CREATE      : err = Create_DataFile( command, dispatchToken->Idx ); break;
          case Command::PRECISION   : err = Precision( command ); break;
          case Command::DATAFILE    : err = DFL_.ProcessDataFileArgs( command ); break;
          case Command::SYSTEM      : system( command.ArgLine() ); break;
          case Command::RUN         : Run(); break;
          case Command::RUN_ANALYSIS:
            // If only 1 arg (the command) run all analyses in list
            if (command.Nargs() == 1) { 
              analysisList_.DoAnalyses();
              mprintf("Analysis complete. Use 'writedata' to write datafiles to disk.\n");
            } else
              err = CrdAnalyze(command);
            break;
          case Command::WRITEDATA   : if (worldrank == 0) DFL_.Write(); break;
          case Command::QUIT        : return C_QUIT; break;
        }
        break;
      case DispatchObject::DEPRECATED: 
        mprintf("Warning: %s is deprecated.\n", command.Command()); 
        break;
      default: mprintf("Dispatch type is currently not handled.\n");
    }
    if (err != 0 && exitOnError_) return C_ERR;
  } else { // Command not recognized
    if (exitOnError_) return C_ERR;
  }
  return C_OK;
}

// -----------------------------------------------------------------------------
// Cpptraj::Run()
int Cpptraj::Run() {
  int err = 0;
  ++nrun_;
  // Special case: check if _DEFAULTCRD_ COORDS DataSet is defined. If so,
  // this means 1 or more actions has requested that a default COORDS DataSet
  // be created.
  DataSet* default_crd = DSL_.FindSetOfType("_DEFAULTCRD_", DataSet::COORDS);
  if (default_crd != 0) {
    mprintf("Warning: One or more analyses requested creation of default COORDS DataSet.\n");
    // If the DataSet has already been written to do not create again.
    if (default_crd->Size() > 0)
      mprintf("Warning: Default COORDS DataSet has already been written to.\n");
    else { 
      if (Dispatch("createcrd _DEFAULTCRD_") == C_ERR) return 1;
    }
  }
  switch ( trajinList_.Mode() ) {
    case TrajinList::NORMAL   : 
      err = RunNormal();
      // Analysis has been run at this point. Clean up Analyses
      analysisList_.Clear();
      break;
    case TrajinList::ENSEMBLE : err = RunEnsemble(); break;
    default:
      // No trajectories loaded; If analyses are defined, try to run them.
      if (!analysisList_.Empty())
        analysisList_.DoAnalyses(); 
      else {
        mprinterr("No trajectories loaded. Exiting.\n");
        err = 1;
      }
  }
  // Clean up Actions.
  actionList_.Clear();
  return err;
}

// Cpptraj::RunEnsemble()
int Cpptraj::RunEnsemble() {
  FrameArray FrameEnsemble;

  mprintf("\nINPUT ENSEMBLE:\n");
  // Ensure all ensembles are of the same size
  int ensembleSize = -1;
  for (TrajinList::const_iterator traj = trajinList_.begin(); traj != trajinList_.end(); ++traj) 
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
  // At this point all ensembles should match (i.e. same map etc.)
  ((Trajin_Multi*)(trajinList_.front()))->EnsembleInfo();

  // Calculate frame division among trajectories
  trajinList_.List();
  int maxFrames = trajinList_.MaxFrames();
  // Parameter file information
  parmFileList_.List();
  // Print reference information 
  mprintf("\nREFERENCE COORDS:\n");
  refFrames_.List();

  // Allocate an ActionList, TrajoutList, and DataSetList for each
  // member of the ensemble. Use separate DataFileList.
  std::vector<ActionList> ActionEnsemble( ensembleSize );
  std::vector<TrajoutList> TrajoutEnsemble( ensembleSize );
  std::vector<DataSetList> DataSetEnsemble( ensembleSize );
  DataFileList DataFileEnsemble;

  // Set up output trajectories for each member of the ensemble
  for (ArgsArray::iterator targ = trajoutArgs_.begin(); targ != trajoutArgs_.end(); ++targ)
  {
    for (int member = 0; member < ensembleSize; ++member) 
      TrajoutEnsemble[member].AddEnsembleTrajout( *targ, parmFileList_, member );
  }
  mprintf("\nENSEMBLE OUTPUT TRAJECTORIES (Numerical filename suffix corresponds to above map):\n");
  TrajoutEnsemble[0].List();
  if (debug_ > 0) {
    for (int member = 1; member < ensembleSize; ++member) {
      mprintf("OUTPUT TRAJECTORIES Member %i:\n", member);
      TrajoutEnsemble[member].List();
    }
  }

  // TODO: One loop over member?
  for (int member = 0; member < ensembleSize; ++member) {
    // Set max frames in the data set list and allocate
    DataSetEnsemble[member].SetMax( maxFrames );
    DataSetEnsemble[member].AllocateSets();
    // Initialize actions 
    if (!actionArgs_.empty())
      mprintf("***** ACTIONS FOR ENSEMBLE MEMBER %i:\n", member);
    for (ArgsArray::iterator aarg = actionArgs_.begin(); aarg != actionArgs_.end(); ++aarg)
    {
      DispatchObject::TokenPtr dispatchToken = Command::SearchToken( *aarg );
      if ( dispatchToken != 0 ) {
        // Create copy of arg list so that args remain unmarked for next member
        ArgList command = *aarg;
        if (ActionEnsemble[member].AddAction( dispatchToken->Alloc, command, &parmFileList_,
                                              &refFrames_, &(DataSetEnsemble[member]), 
                                              &DataFileEnsemble ))
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
  rprintf("\nBEGIN ENSEMBLE PROCESSING:\n");
  for ( TrajinList::const_iterator traj = trajinList_.begin();
                                   traj != trajinList_.end(); ++traj)
  {
    // Open up the trajectory file. If an error occurs, bail 
    if ( (*traj)->BeginTraj(showProgress_) ) {
      mprinterr("Error: Could not open trajectory %s.\n",(*traj)->TrajFilename().full());
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
      CurrentParm->SetReferenceCoords( refFrames_.ActiveReference() );
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
      lastPindex = CurrentParm->Pindex();
    }

    // Loop over every collection of frames in the ensemble
    (*traj)->PrintInfoLine();
    Trajin_Multi* mtraj = (Trajin_Multi*)*traj;
    while ( mtraj->GetNextEnsemble(FrameEnsemble) ) {
      if (!mtraj->BadEnsemble()) {
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
      } else {
        mprinterr("Error: Could not read frame %i for ensemble.\n", actionSet + 1);
      }
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
  mprintf("\nENSEMBLE ACTION OUTPUT:\n");
  for (int member = 0; member < ensembleSize; ++member)
    actionList_.Print( );

  // Sync DataSets and print DataSet information
  // TODO - Also have datafilelist call a sync??
  int total_data_sets = DataSetEnsemble[0].size();
  mprintf("\nENSEMBLE DATASETS: Each member has %i sets total.\n", total_data_sets);
  for (int member = 0; member < ensembleSize; ++member) {
    DataSetEnsemble[member].Sync();
    DataSetEnsemble[member].sort();
    if (total_data_sets != DataSetEnsemble[member].size())
      mprintf("Warning: Ensemble member %i # data sets (%i) does not match member 0 (%i)\n",
              member, DataSetEnsemble[member].size(), total_data_sets);
    if (debug_ > 0)
      DataSetEnsemble[member].List();
  }

  // Print Datafile information
  DataFileEnsemble.List();
  // Only Master does DataFile output
  if (worldrank==0)
    DataFileEnsemble.Write();

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

  // ========== S E T U P   P H A S E ========== 
  // Parameter file information
  parmFileList_.List();
  // Input coordinate file information
  trajinList_.List();
  // Print reference information 
  mprintf("\nREFERENCE COORDS:\n");
  refFrames_.List();
  // Output traj
  mprintf("\nOUTPUT TRAJECTORIES:\n");
  trajoutList_.List();
  // Allocate DataSets in the master DataSetList based on # frames to be read
  DSL_.AllocateSets(); 
  
  // ========== A C T I O N  P H A S E ==========
  // Loop over every trajectory in trajFileList
  rprintf("\nBEGIN TRAJECTORY PROCESSING:\n");
  for ( TrajinList::const_iterator traj = trajinList_.begin();
                                   traj != trajinList_.end(); ++traj)
  {
    // Open up the trajectory file. If an error occurs, bail 
    if ( (*traj)->BeginTraj(showProgress_) ) {
      mprinterr("Error: Could not open trajectory %s.\n",(*traj)->TrajFilename().full());
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
      CurrentParm->SetReferenceCoords( refFrames_.ActiveReference() );
      // Set up actions for this parm
      if (actionList_.SetupActions( &CurrentParm )) {
        mprintf("WARNING: Could not set up actions for %s: skipping.\n",
                CurrentParm->c_str());
        continue;
      }
      lastPindex = CurrentParm->Pindex();
    }

    // Loop over every Frame in trajectory
    (*traj)->PrintInfoLine();
    while ( (*traj)->GetNextFrame(TrajFrame) ) {
      // Since Frame can be modified by actions, save original and use CurrentFrame
      Frame* CurrentFrame = &TrajFrame;
      // Perform Actions on Frame
      bool suppress_output = actionList_.DoActions(&CurrentFrame, actionSet);
      // Do Output
      if (!suppress_output)
        trajoutList_.Write(actionSet, CurrentParm, CurrentFrame);
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
  trajoutList_.Close();

  // ========== A C T I O N  O U T P U T  P H A S E ==========
  mprintf("\nACTION OUTPUT:\n");
  actionList_.Print( );

  // Sync DataSets and print DataSet information
  DSL_.Sync();

  // ========== A N A L Y S I S  P H A S E ==========
  mprintf("\nDATASETS:\n");
  if (!analysisList_.Empty()) {
    DSL_.List();
    analysisList_.DoAnalyses();
    // DEBUG: DataSets, post-Analysis
    mprintf("\nDATASETS AFTER ANALYSIS:\n");
  }
  DSL_.List();

  // ========== D A T A  W R I T E  P H A S E ==========
  // Print Datafile information
  DFL_.List();
  // Only Master does DataFile output
  if (worldrank==0)
    DFL_.Write();
 
  return 0;
}
