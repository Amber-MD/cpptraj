#include <cstdio> // for ReadInput
#include <cstdlib> // system
#include "Command.h"
#include "CpptrajStdio.h"
// INC_ACTION==================== ALL ACTION CLASSES GO HERE ===================
#include "Action_Distance.h"
#include "Action_Rmsd.h"
#include "Action_Dihedral.h"
#include "Action_Angle.h"
#include "Action_AtomMap.h"
#include "Action_Strip.h"
#include "Action_DSSP.h"
#include "Action_Center.h"
#include "Action_Hbond.h"
#include "Action_Image.h"
#include "Action_Surf.h"
#include "Action_Radgyr.h"
#include "Action_Mask.h"
#include "Action_Closest.h"
#include "Action_NAstruct.h"
#include "Action_Pucker.h"
#include "Action_Outtraj.h"
#include "Action_Average.h"
#include "Action_Radial.h"
#include "Action_DistRmsd.h"
#include "Action_Jcoupling.h"
#include "Action_Pairwise.h"
#include "Action_Molsurf.h"
#include "Action_CheckStructure.h"
#include "Action_DihedralScan.h"
#include "Action_Rotdif.h"
#include "Action_RunningAvg.h"
#include "Action_AtomicFluct.h"
#include "Action_Watershell.h"
#include "Action_Contacts.h"
#include "Action_Vector.h"
#include "Action_Principal.h"
#include "Action_Matrix.h"
#include "Action_LIE.h"
#include "Action_Grid.h"
#include "Action_GridFreeEnergy.h"
#include "Action_Dipole.h"
#include "Action_Projection.h"
#include "Action_ClusterDihedral.h"
#include "Action_Unwrap.h"
#include "Action_Diffusion.h"
#include "Action_DNAionTracker.h"
#include "Action_Scale.h"
#include "Action_RandomizeIons.h"
#include "Action_AutoImage.h"
#include "Action_STFC_Diffusion.h"
#include "Action_AtomicCorr.h"
#include "Action_Bounds.h"
#include "Action_Rotate.h"
#include "Action_Translate.h"
#include "Action_Box.h"
#include "Action_CreateCrd.h"
#include "Action_MultiDihedral.h"
#include "Action_MakeStructure.h"
#include "Action_SymmetricRmsd.h"
#include "Action_Volmap.h"
#include "Action_Spam.h"
#include "Action_Temperature.h"
#include "Action_CreateReservoir.h"
#include "Action_Density.h"
#include "Action_PairDist.h"
#include "Action_OrderParameter.h"
#include "Action_MinDist.h"
#include "Action_FixAtomOrder.h"

// INC_ANALYSIS================= ALL ANALYSIS CLASSES GO HERE ==================
#include "Analysis_Hist.h"
#include "Analysis_Corr.h"
#include "Analysis_Matrix.h"
#include "Analysis_Timecorr.h"
#include "Analysis_IRED.h"
#include "Analysis_Modes.h"
#include "Analysis_CrankShaft.h"
#include "Analysis_Statistics.h"
#include "Analysis_CrossCorr.h"
#include "Analysis_AutoCorr.h"
#include "Analysis_Lifetime.h"
#include "Analysis_FFT.h"
#include "Analysis_CrdFluct.h"
#include "Analysis_RmsAvgCorr.h"
#include "Analysis_Rms2d.h"
#include "Analysis_Clustering.h"
#include "Analysis_RunningAvg.h"
#include "Analysis_MeltCurve.h"
#include "Analysis_Overlap.h"
#include "Analysis_AmdBias.h"
#include "Analysis_RemLog.h"
// ---- Command Functions ------------------------------------------------------
/** Search Commands list for a specific type of command. */
Command::TokenPtr Command::SearchTokenType(CommandType dtype, 
                                                   ArgList const& argIn)
{
  for (TokenPtr token = Commands; token->Type != NONE; ++token)
  {
    if (dtype != token->Type) continue;
    if (argIn.CommandIs( token->Cmd )) return token;
  }
  mprintf("[%s]: Command not found.\n", argIn.Command());
  return 0;
}

/// Strings that correspond to CommandType
const char* Command::CommandTitle[] = { 0, "Topology", "Trajectory", "Action",
  "Analysis", "General", "Deprecated" };

/** List all commands of the given type, or all commands if type
  * is NONE.
  */
void Command::ListCommands(CommandType dtype) {
  CommandType lastType = NONE;
  int col = 0;
  for (TokenPtr token = Commands; token->Type != DEPRECATED; ++token)
  {
    CommandType currentType = token->Type;
    if (dtype != NONE && dtype != currentType) continue;
    // Command type title
    if (currentType != lastType) {
      if (col != 0) mprintf("\n");
      mprintf("%s Commands:\n", CommandTitle[currentType]);
      lastType = currentType;
      col = 0;
    }
    if (col == 0) mprintf("\t");
    mprintf("%s  ", token->Cmd);
    ++col;
    if (col == 8) {
      mprintf("\n");
      col = 0;
    }
  }
  mprintf("\n");
}

/** Search the Commands list for given command.
  * \return the token if found, 0 if not.
  */
Command::TokenPtr Command::SearchToken(ArgList& argIn) {
  // SPECIAL CASE: For backwards compat. remove analyze prefix
  if (argIn.CommandIs("analyze")) {
    argIn.RemoveFirstArg();
    argIn.MarkArg(0); // Mark new first arg as command
    return (SearchTokenType(ANALYSIS, argIn));
  }
  // Search for command.
  for (TokenPtr token = Commands; token->Type != NONE; ++token)
    if (argIn.CommandIs( token->Cmd )) return token;
  mprintf("[%s]: Command not found.\n", argIn.Command());
  return 0;
}

/** Search for the given command and execute it. */
Command::RetType Command::Dispatch(CpptrajState& State, std::string const& commandIn) {
  ArgList cmdArg( commandIn );
  cmdArg.MarkArg(0); // Always mark the first arg as the command 
  TokenPtr cmdToken = SearchToken( cmdArg );
  if (cmdToken == 0) return C_ERR;
  return ( cmdToken->Fxn( State, cmdArg, cmdToken->Alloc ) );
}

/// Used by ProcessInput to determine when line ends.
static inline bool EndChar(char ptr) {
  if (ptr=='\n' || ptr=='\r' || ptr=='\0' || ptr==EOF) return true;
  return false;
}

/** Read commands from an input file, or from STDIN if given filename
  * is empty. '#' indicates the beginning of a comment, backslash at the 
  * end of a line indicates continuation (otherwise indicates 'literal').
  * \return 0 if successfully read, 1 on error.
  */
Command::RetType Command::ProcessInput(CpptrajState& State, 
                                               std::string const& inputFilename)
{
  FILE* infile; // TODO: CpptrajFile
  if (inputFilename.empty()) {
    mprintf("INPUT: Reading Input from STDIN\n");
    infile = stdin;
  } else {
    mprintf("INPUT: Reading Input from file %s\n",inputFilename.c_str());
    if ( (infile=fopen(inputFilename.c_str(),"r"))==0 ) {
      rprinterr("Error: Could not open input file %s\n",inputFilename.c_str());
      return C_ERR;
    }
  }
  // Read in each line of input. Newline or null terminates. \ continues line.
  std::string inputLine;
  unsigned int idx = 0;
  char lastchar = '0';
  char ptr = 0;
  int nInputErrors = 0;
  RetType cmode = C_OK;
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
      cmode = Command::Dispatch(State, inputLine);
      if (cmode == C_ERR)
        nInputErrors++;
      else if (cmode == C_QUIT)
        break;
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
  if (nInputErrors > 0) {
    mprinterr("\t%i errors encountered reading input.\n", nInputErrors);
    return C_ERR;
  }
  return cmode;
}

// ====================== CPPTRAJ COMMANDS HELP ================================
static void Help_Help() {
  mprintf("\t{[<cmd>] | General | Action | Analysis | Topology | Trajectory}\n");
  mprintf("\tWith no arguments list all known commands, otherwise display help for\n");
  mprintf("\tcommand <cmd>. If General/Action/Analysis/Topology/Trajectory specified\n");
  mprintf("\tlist commands only in that category.\n");
}

static void Help_System() {
  mprintf("\tCall command from system.\n");
}

static void Help_NoProgress() {
  mprintf("\tDo not print progress while reading in trajectories.\n");
}

static void Help_NoExitOnError() {
  mprintf("\tDo not exit when errors are encountered. This is the default\n");
  mprintf("\tin interactive mode.\n");
}

static void Help_Run() {
  mprintf("\tProcess all trajectories currently in input trajectory list.\n");
  mprintf("\tAll actions in action list will be run on each frame.\n");
  mprintf("\tIf not processing ensemble input, all analyses in analysis\n");
  mprintf("\tlist will be run after trajectory processing.\n");
}

static void Help_Quit() {
  mprintf("\tExit CPPTRAJ\n");
}

static const char TypeList[] =
  "(<type> = actions,trajin,trajout,ref,parm,analysis,datafile,dataset)";

static void Help_List() {
  mprintf("\t[<type>] %s\n", TypeList);
  mprintf("\tList currently loaded objects of the specified type. If no type is given\n");
  mprintf("\tthen list all loaded objects.\n");
}

static void Help_Debug() {
  mprintf("\t[<type>] <#> %s\n", TypeList);
  mprintf("\tSet debug level for new objects of the specified type. If no type is given\n");
  mprintf("\tthen set debug level for all new objects. Does not affect current objects.\n");
}

static void Help_Clear() {
  mprintf("\t[ {all | <type>} ] %s\n", TypeList);
  mprintf("\tClear currently loaded objects of the specified type. If 'all' is specified\n");
  mprintf("\tthen clear all loaded objects.\n");
}

static void Help_ActiveRef() {
  mprintf("\t<#>\n");
  mprintf("\tSet the reference structure to be used for coordinate-based mask parsing.\n");
  mprintf("\t<#> starts from 0 (first loaded reference).\n");
}

static void Help_Create_DataFile() {
  mprintf("\t<filename> <dataset0> [<dataset1> ...]\n");
  mprintf("\tAdd a file with specified data sets to the data file list. Does not\n");
  mprintf("\timmediately write the data.\n");
}

static void Help_DataFile() {
  mprintf("\t<data filename> <datafile cmd>\n");
  mprintf("\tPass <datafile cmd> to specified data file currently in data file list.\n");
}

static void Help_ReadData() {
  mprintf("\t<filename>\n");
  mprintf("\tRead data from <filename> into data sets.\n");
}

static void Help_ReadInput() {
  mprintf("\t<filename>\n");
  mprintf("\tRead commands from <filename>\n");
}

static void Help_Write_DataFile() {
  mprintf("\t<filename> <dataset0> [<dataset1> ...]\n");
  mprintf("\tWrite specified data sets to <filename> immediately.\n");
}

static void Help_WriteData() {
  mprintf("\tWrite all files currently in the data file list.\n");
}

static void Help_Precision() {
  mprintf("\t{<filename> | <dataset arg>} [<width>] [<precision>]\n");
  mprintf("\tSet precision for all datasets in datafile <filename> or dataset(s)\n");
  mprintf("\tspecified by <dataset arg> to <width>.<precision>. If width/precision\n");
  mprintf("\tis not specified then default to 12.4\n");
}

static void Help_Select() {
  mprintf("\t[<parmindex>] <mask>\n");
  mprintf("\tShow atom numbers selected by <mask> for parm <parmindex>\n");
  mprintf("\t(default first parm)\n");
}

static void Help_SelectDS() {
  mprintf("\t<dataset selection>\n");
  mprintf("\tShow results of data set selection. Data set selection format is:\n");
  mprintf("\t\t<name>[<aspect]:<idx range>\n");
  mprintf("\tWhere '<name>' is the data set name, '[<aspect>]' is the data set aspect,\n");
  mprintf("\tand <idx range> is a numerical range specifying data set indices (i.e. 2-5,7 etc).\n");
  mprintf("\tThe aspect and index portions may be optional. An asterisk '*' may be used as\n");
  mprintf("\ta wildcard. E.g. 'selectds R2', 'selectds RoG[Max]', 'selectds PR[res]:2-12'\n");
}

static void Help_Trajin() {
  mprintf("\t<filename> {[<start>] [<stop> | last] [offset]} | lastframe\n");
  mprintf("\t           [parm <parmfile> | parmindex <#>]\n");
  mprintf("\t           [ remdtraj [remdtrajtemp <T> | remdtrajidx <#>]\n");
  mprintf("\t           [trajnames <rep1>,<rep2>,...,<repN> ] ]\n");
  mprintf("\tLoad trajectory specified by <filename> to the input trajectory list.\n");
}

static void Help_Ensemble() {
  mprintf("\t<file0> {[<start>] [<stop> | last] [offset]} | lastframe\n");
  mprintf("\t        [parm <parmfile> | parmindex <#>]\n");
  mprintf("\t        [trajnames <file1>,<file2>,...,<fileN>\n");
  mprintf("\tLoad an ensemble of trajectories starting with <file0> that will be processed together.\n");
}

static void Help_Trajout() {
  mprintf("\t<filename> [<fileformat>] [append] [nobox]\n");
  mprintf("\t           [parm <parmfile> | parmindex <#>] [onlyframes <range>] [title <title>]\n");
  mprintf("\t           %s\n", ActionFrameCounter::HelpText);
  mprintf("\t           [ <Format Options> ]\n");
  mprintf("\tSpecify output trajectory.\n");
}

static void Help_Reference() {
  mprintf("\t<filename> [<frame#>] [<mask>] [TAG] [lastframe]\n");
  mprintf("\t           [average [<stop>] [<offset>]]\n");
  mprintf("\tLoad trajectory <filename> as a reference frame.\n");
}

static void Help_Parm() {
  mprintf("\t<filename> [<tag>] [nobondsearch | bondsearch [<offset>]]\n");
  mprintf("\tAdd <filename> to the topology list.\n");
}
static void Help_ParmInfo() {
  mprintf("\t[<parmindex>] [<mask>]\n");
  mprintf("\tPrint information on topology <parmindex> (0 by default). If <mask> is given\n");
  mprintf("\tprint info on atoms in mask. If no mask given print overall information.\n");
}

static void Help_ParmWrite() {
  mprintf("\tout <filename> [<parmindex>]\n");
  mprintf("\tWrite topology <parmindex> to <filename> as an Amber Topology file.\n");
}

static void Help_ParmStrip() {
  mprintf("\t<mask> [<parmindex>]\n");
  mprintf("\tStrip atoms in mask from topology <parmindex>.\n");
}

static void Help_ParmBox() {
  mprintf("\t[<parmindex>] [x <xval>] [y <yval>] [z <zval>]\n");
  mprintf("\t              [alpha <a>] [beta <b>] [gamma <g>] [nobox]\n");
  mprintf("\tSet the specified topology box info to what is specified. If nobox, remove box info.\n");
}

static void Help_Solvent() {
  mprintf("\t[<parmindex>] { <mask> | none }\n");
  mprintf("\tSet solvent for the specified topology (default 0) based on <mask>.\n");
  mprintf("\tIf 'none' specified, remove all solvent information.\n");
}

static void Help_BondInfo() {
  mprintf("\t[<parmindex>] [<mask>]\n");
  mprintf("\tPrint bond information of atoms in <mask> for topology <parmindex> (0 by default).\n");
}

static void Help_ChargeInfo() {
  mprintf("\t[<parmindex>] <mask>\n");
  mprintf("\tPrint the total charge of atoms in <mask> for topology <parmindex> (0 by default).\n");
}

static void Help_ResInfo() {
  mprintf("\t[<parmindex>] [<mask>]\n");
  mprintf("\tPrint information for residues in <mask> for topology <parmindex> (0 by default).\n");
}
static void Help_MolInfo() {
  mprintf("\t[<parmindex>] [<mask>]\n");
  mprintf("\tPrint information for molecules in <mask> for topology <parmindex> (0 by default).\n");
}

static void Help_LoadCrd() {
  mprintf("\t<filename> [parm <parm> | parmindex<#>] [<trajin args>] [<name>]\n");
  mprintf("\tLoad trajectory <filename> as a COORDS data set named <name> (default <filename>).\n");
}
static void Help_CrdAction() {
  mprintf("\t<crd set> <actioncmd> [<action args>] [crdframes <start>,<stop>,<offset>]\n");
  mprintf("\tPerform action <actioncmd> on COORDS data set <crd set>.\n");
}

static void Help_CrdOut() {
  mprintf("\t<crd set> <filename> [<trajout args>] [crdframes <start>,<stop>,<offset>]\n");
  mprintf("\tWrite COORDS data set <crd set> to trajectory file <filename>\n");
}

static void Help_RunAnalysis() {
  mprintf("\t[<analysis> [<analysis args>]]\n");
  mprintf("\tIf specified alone, run all analyses in the analysis list.\n");
  mprintf("\tOtherwise run the specified analysis immediately.\n");
}

// ---------- GENERAL COMMANDS -------------------------------------------------
/// Set active reference for distance-based masks etc.
Command::RetType ActiveRef(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)State.FL()->SetActiveRef( argIn.getNextInteger(0) );
}

/// Clear data in specified lists
Command::RetType ClearList(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)State.ClearList( argIn );
}

/// Set debug value for specified list(s)
Command::RetType SetListDebug(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)State.SetListDebug( argIn );
}

/// List all members of specified list(s)
Command::RetType ListAll(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)State.ListAll( argIn );
}

/// Perform action on given COORDS dataset
Command::RetType CrdAction(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  std::string setname = argIn.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: %s: Specify COORDS dataset name.\n", argIn.Command());
    return Command::C_ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL()->FindSetOfType( setname, DataSet::COORDS );
  if (CRD == 0) {
    mprinterr("Error: %s: No COORDS set with name %s found.\n", argIn.Command(), setname.c_str());
    return Command::C_ERR;
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
  if (State.Debug() > 0) mprintf("\tDBG: Frames %i to %i, offset %i\n", start+1, stop, offset);
  ArgList actionargs = argIn.RemainingArgs();
  actionargs.MarkArg(0);
  Command::TokenPtr tkn = Command::SearchTokenType( Command::ACTION, actionargs);
  if ( tkn == 0 ) return Command::C_ERR;
  Action* act = (Action*)tkn->Alloc();
  if (act == 0) return Command::C_ERR;
  if ( act->Init( actionargs, State.PFL(), State.FL(), State.DSL(), State.DFL(), State.Debug() ) != Action::OK ) {
    delete act;
    return Command::C_ERR;
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
    return Command::C_ERR;
  }
  // Check if parm was modified. If so, update COORDS.
  if ( currentParm != originalParm ) {
    mprintf("Info: crdaction: Parm for %s was modified by action %s\n",
            CRD->Legend().c_str(), actionargs.Command());
    CRD->SetTopology( *currentParm );
  }
  // Loop over all frames in COORDS.
  ProgressBar progress( stop - start );
  int set = 0;
  for (int frame = start; frame < stop; frame += offset) {
    progress.Update( set );
    CRD->GetFrame( frame, *originalFrame );
    Frame* currentFrame = originalFrame;
    if (act->DoAction( set, currentFrame, &currentFrame ) == Action::ERR) {
      mprinterr("Error: crdaction: Frame %i, set %i\n", frame + 1, set + 1);
      break;
    }
    // Check if frame was modified. If so, update COORDS.
    // TODO: Have actions indicate whether they will modify coords
    //if ( currentFrame != originalFrame ) 
      CRD->SetCRD( frame, *currentFrame );
    set++;
  }
  act->Print();
  State.MasterDataFileWrite();
  delete originalFrame;
  delete originalParm;
  delete act;
  return Command::C_OK;
}

/// Write out COORDS dataset
Command::RetType CrdOut(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  std::string setname = argIn.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: crdout: Specify COORDS dataset name.\n");
    return Command::C_ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL()->FindSetOfType( setname, DataSet::COORDS );
  if (CRD == 0) {
    mprinterr("Error: crdout: No COORDS set with name %s found.\n", setname.c_str());
    return Command::C_ERR;
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
  if (State.Debug() > 0) mprintf("\tDBG: Frames %i to %i, offset %i\n", start+1, stop, offset);
  Trajout outtraj;
  Topology* currentParm = (Topology*)&(CRD->Top()); // TODO: Fix cast
  if (outtraj.InitTrajWrite( setname, &argIn, currentParm, TrajectoryFile::UNKNOWN_TRAJ)) {
    mprinterr("Error: crdout: Could not set up output trajectory.\n");
    return Command::C_ERR;
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
  return Command::C_OK;
}

/// Load single trajectory as DataSet_Coords
Command::RetType LoadCrd(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  // Get parm
  Topology* parm = State.PFL()->GetParm( argIn );
  if (parm == 0) {
    mprinterr("Error: loadcrd: No parm files loaded.\n");
    return Command::C_ERR;
  }
  // Load trajectory
  Trajin_Single trajin;
  trajin.SetDebug( State.Debug() );
  if (trajin.SetupTrajRead(argIn.GetStringNext(), &argIn, parm)) {
    mprinterr("Error: loadcrd: Could not set up input trajectory.\n");
    return Command::C_ERR;
  }
  // Create input frame
  Frame frameIn;
  frameIn.SetupFrameV(parm->Atoms(), trajin.HasVelocity(), trajin.NreplicaDimension());
  // Create DataSet, use base file name as set name if none specified. 
  // NOTE: Default name should NEVER get used.
  std::string setname = argIn.GetStringNext();
  if (setname.empty())
    setname = trajin.TrajFilename().Base();
  DataSet_Coords* coords = (DataSet_Coords*)State.DSL()->AddSet(DataSet::COORDS, 
                                                                setname, "__DCRD__");
  if (coords == 0) {
    mprinterr("Error: loadcrd: Could not set up COORDS data set.\n");
    return Command::C_ERR;
  }
  coords->SetTopology( *parm );
  // Read trajectory
  mprintf("\tLoading trajectory %s as \"%s\"\n", trajin.TrajFilename().full(), setname.c_str());
  trajin.BeginTraj(true);
  trajin.PrintInfoLine();
  while (trajin.GetNextFrame( frameIn ))
    coords->AddFrame( frameIn );
  trajin.EndTraj();
  return Command::C_OK;
}

/// Add DataSets specified by arguments to given DataFile.
// NOTE: Used byt Create_DataFile and Write_DataFile
// TODO: Put in DataFile?
static int AddSetsToDataFile(DataFile& df, ArgList const& dsetArgs, DataSetList& DSL)
{
  int err = 0;
  for (ArgList::const_iterator dsa = dsetArgs.begin(); dsa != dsetArgs.end(); ++dsa) {
    DataSetList Sets = DSL.GetMultipleSets( *dsa );
    if (Sets.empty())
      mprintf("Warning: %s does not correspond to any data sets.\n", (*dsa).c_str());
    for (DataSetList::const_iterator set = Sets.begin(); set != Sets.end(); ++set) {
      mprintf(" %s", (*set)->Legend().c_str());
      if ( df.AddSet(*set) ) {
        mprinterr("Error: Could not add data set %s to file.\n", (*set)->Legend().c_str());
        ++err;
      }
    }
  }
  mprintf("\n");
  return err;
}

/// Add a new DataFile to DFL with specified DataSets, to be written later.
Command::RetType Create_DataFile(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  // Next string is datafile that command pertains to.
  std::string name1 = argIn.GetStringNext();
  if (name1.empty()) {
    mprinterr("Error: No filename given.\n");
    return Command::C_ERR;
  }
  DataFile* df = State.DFL()->AddDataFile(name1, argIn);
  if (df == 0) return Command::C_ERR;
  return (Command::RetType)( AddSetsToDataFile(*df, argIn.RemainingArgs(), *(State.DSL())) );
}

/// Write DataFile with specified DataSets immediately.
Command::RetType Write_DataFile(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  // Next string is datafile that command pertains to.
  std::string name1 = argIn.GetStringNext();
  if (name1.empty()) {
    mprinterr("Error: No filename given.\n");
    return Command::C_ERR;
  }
  DataFile* df = new DataFile();
  if (df == 0) return Command::C_ERR;
  if (df->SetupDatafile( name1, argIn, State.Debug() )) {
    delete df;
    return Command::C_ERR;
  }
  int err = AddSetsToDataFile(*df, argIn.RemainingArgs(), *(State.DSL()));
  if (err == 0) df->WriteData();
  delete df;
  return (Command::RetType)err;
}

/// Process DataFile-specific command
Command::RetType DataFileCmd(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)( State.DFL()->ProcessDataFileArgs( argIn ) );
}

/// Read data from file into master DataSetList
Command::RetType ReadData(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  DataFile dataIn;
  if (dataIn.ReadData( argIn, *State.DSL() )!=0) {
    mprinterr("Error: Could not read data file.\n");
    return Command::C_ERR;
  }
  return Command::C_OK;
}

/// Exit
Command::RetType Quit(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return Command::C_QUIT;
}

/// Run a system command
Command::RetType SystemCmd(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  system( argIn.ArgLine() );
  return Command::C_OK;
}

/// Find help for command/topic
Command::RetType Help(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  ArgList arg = argIn;
  arg.RemoveFirstArg();
  if (arg.empty())
    // NONE in this context means list all commands
    Command::ListCommands(Command::NONE);
  else if (arg.CommandIs("General"))
    Command::ListCommands(Command::GENERAL);
  else if (arg.CommandIs("Topology"))
    Command::ListCommands(Command::PARM);
  else if (arg.CommandIs("Action"))
    Command::ListCommands(Command::ACTION);
  else if (arg.CommandIs("Analysis"))
    Command::ListCommands(Command::ANALYSIS);
  else if (arg.CommandIs("Trajectory"))
    Command::ListCommands(Command::TRAJ);
  else {
    Command::TokenPtr dispatchToken = Command::SearchToken( arg );
    if (dispatchToken == 0 || dispatchToken->Help == 0)
      mprinterr("No help found for %s\n", arg.Command());
    else
      dispatchToken->Help();
  }
  return Command::C_OK;
}

/// Run the current State
Command::RetType RunState(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)State.Run();
}

/// Read input from a file.
Command::RetType ReadInput(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  // Next arg should be a filename. Not allowed to be blank in command.
  std::string inputFilename = argIn.GetStringNext();
  if (inputFilename.empty()) {
    mprinterr("Error: No input filename given.\n");
    return Command::C_ERR;
  }
  return Command::ProcessInput(State, inputFilename);
}

/// Tell CpptrajState to ignore errors if possible
Command::RetType NoExitOnError(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  State.SetNoExitOnError();
  mprintf("\tAttempting to ignore errors if possible.\n");
  return Command::C_OK;
}

/// Tell CpptrajState not to use a progress bar during Run.
Command::RetType NoProgress(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  State.SetNoProgress();
  mprintf("\tProgress bar will not be used during Run.\n");
  return Command::C_OK;
}

///  Set precision for specific set or all sets in specified DataFile
Command::RetType Precision(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  // Next string is DataSet(s)/DataFile that command pertains to.
  std::string name1 = argIn.GetStringNext();
  if (name1.empty()) {
    mprinterr("Error: precision: No filename/setname given.\n");
    return Command::C_ERR;
  }
  // This will break if dataset name starts with a digit...
  int width = argIn.getNextInteger(12);
  if (width < 1) {
    mprintf("Error: precision: Cannot set width < 1 (%i).\n", width);
    return Command::C_ERR;
  }
  int precision = argIn.getNextInteger(4);
  if (precision < 0) precision = 0;
  DataFile* df = State.DFL()->GetDataFile(name1);
  if (df != 0) {
    mprintf("\tSetting precision for all sets in %s to %i.%i\n", df->DataFilename().base(),
            width, precision);
    df->SetPrecision(width, precision);
  } else {
    DataSetList dsets = State.DSL()->GetMultipleSets( name1 );
    mprintf("\tSetting precision for %i sets to %i.%i\n", dsets.size(),
            width, precision);
    for (DataSetList::const_iterator set = dsets.begin(); set != dsets.end(); ++set)
      (*set)->SetPrecision(width, precision);
  }
  return Command::C_OK;
}

/// Run specified analysis or all analyses in State.
Command::RetType RunAnalysis(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  // If only 1 arg (the command) run all analyses in list
  if (argIn.Nargs() == 1) {
    State.RunAnalyses();
    State.MasterDataFileWrite();
    return Command::C_OK;
  }
  // Run specified analysis
  // FIXME: Use RemoveFirstArg
  ArgList analyzeargs = argIn.RemainingArgs();
  analyzeargs.MarkArg(0);
  Command::TokenPtr tkn = Command::SearchTokenType( Command::ANALYSIS, analyzeargs );
  if ( tkn == 0 ) return Command::C_ERR;
  Analysis* ana = (Analysis*)tkn->Alloc();
  if (ana == 0) return Command::C_ERR;
  Command::RetType err = Command::C_ERR;
  if ( ana->Setup( analyzeargs, State.DSL(), State.PFL(), State.DFL(), State.Debug() ) == 
                   Analysis::OK )
  {
    if (ana->Analyze() != Analysis::ERR) {
      err = Command::C_OK;
      State.MasterDataFileWrite();
    }
  }
  delete ana;
  return err;
}

/// Show results of mask expression
Command::RetType SelectAtoms(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  AtomMask tempMask( argIn.GetMaskNext() );
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  parm->SetupIntegerMask( tempMask );
  mprintf("Selected %i atoms.\n", tempMask.Nselected());
  if (!argIn.hasKey("total"))
    tempMask.PrintMaskAtoms("Selected");
  return Command::C_OK;
}

/// Show results of DataSet expression
Command::RetType SelectDataSets(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  std::string dsarg = argIn.GetStringNext();
  DataSetList dsets = State.DSL()->GetMultipleSets( dsarg );
  mprintf("SelectDS: Arg [%s]:", dsarg.c_str());
  dsets.List();
  return Command::C_OK;
}

/// Write all DataFiles in State
Command::RetType WriteAllData(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  State.MasterDataFileWrite();
  return Command::C_OK;
}

// ---------- TRAJECTORY COMMANDS ----------------------------------------------
/// Add output trajectory to State
Command::RetType Trajout(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)State.AddTrajout( argIn );
}

/// Add input trajectory to State
Command::RetType Trajin(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)State.AddTrajin( argIn );
}

/// Add reference trajectory to State
Command::RetType Reference(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)State.AddReference( argIn );
}

// ---------- TOPOLOGY COMMANDS ------------------------------------------------
/// Load topology to State
Command::RetType LoadParm(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  std::string parmtag = argIn.getNextTag();
  bool bondsearch = !argIn.hasKey("nobondsearch");
  double offset = argIn.getKeyDouble("bondsearch", -1.0);
  return (Command::RetType)
    State.PFL()->AddParmFile(argIn.GetStringNext(), parmtag, bondsearch, offset);
}

/// Print info for specified parm or atoms in specified parm.
Command::RetType ParmInfo(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  std::string maskarg = argIn.GetMaskNext();
  if (!maskarg.empty())
    parm->PrintAtomInfo( maskarg );
  else
    parm->Summary();
  return Command::C_OK;
}

/// Print bond info for atoms in mask.
Command::RetType BondInfo(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  parm->PrintBondInfo( argIn.GetMaskNext() );
  return Command::C_OK;
}

/// Print residue info for atoms in mask.
Command::RetType ResInfo(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  parm->PrintResidueInfo( argIn.GetMaskNext() );
  return Command::C_OK;
}

/// Print molecule info for atoms in mask.
Command::RetType MolInfo(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  parm->PrintMoleculeInfo( argIn.GetMaskNext() );
  return Command::C_OK;
}

/// Print the total charge of atoms in mask
Command::RetType ChargeInfo(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  parm->PrintChargeInfo( argIn.GetMaskNext() );
  return Command::C_OK;
}

/// Modify specified parm box info
Command::RetType ParmBox(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Box pbox;
  bool nobox = false;
  if ( argIn.hasKey("nobox") )
    nobox = true;
  else {
    pbox.SetX( argIn.getKeyDouble("x",0) );
    pbox.SetY( argIn.getKeyDouble("y",0) );
    pbox.SetZ( argIn.getKeyDouble("z",0) );
    pbox.SetAlpha( argIn.getKeyDouble("alpha",0) );
    pbox.SetBeta(  argIn.getKeyDouble("beta",0)  );
    pbox.SetGamma( argIn.getKeyDouble("gamma",0) );
  }
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  if (nobox)
    mprintf("\tRemoving box information from parm %i:%s\n", parm->Pindex(), parm->c_str());
  else
    // Fill in missing parm box information from specified parm
    pbox.SetMissingInfo( parm->ParmBox() );
  parm->SetBox( pbox );
  return Command::C_OK;
}

/// Strip atoms from specified parm
Command::RetType ParmStrip(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  AtomMask tempMask( argIn.GetMaskNext() );
  // Since want to keep atoms outside mask, invert selection
  tempMask.InvertMask();
  parm->SetupIntegerMask( tempMask );
  mprintf("\tStripping atoms in mask [%s] (%i) from %s\n",tempMask.MaskString(),
           parm->Natom() - tempMask.Nselected(), parm->c_str());
  Topology* tempParm = parm->modifyStateByMask(tempMask);
  if (tempParm==0) {
    mprinterr("Error: %s: Could not strip parm.\n", argIn.Command());
    return Command::C_ERR;
  } else {
    // Replace parm with stripped version
    // TODO: Implement proper assignment op for Topology
    tempParm->ParmInfo();
    mprintf("\n");
    State.PFL()->ReplaceParm(parm->Pindex(), tempParm);
  }
  return Command::C_OK;
}

/// Write parm to Amber Topology file.
Command::RetType ParmWrite(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  if (State.PFL()->WriteParm( argIn )) return Command::C_ERR;
  return Command::C_OK;
}

/// Modify parm solvent information
Command::RetType ParmSolvent(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  std::string maskexpr;
  if (!argIn.hasKey("none")) {
    maskexpr = argIn.GetMaskNext();
    if ( maskexpr.empty() ) {
      mprinterr("Error: solvent: No mask specified.\n");
      return Command::C_ERR;
    }
  }
  // Get parm index
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  parm->SetSolvent( maskexpr );
  return Command::C_OK;
}

/// Scale dihedral force constants in specfied parm by factor.
Command::RetType ScaleDihedralK(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Topology* parm = State.PFL()->GetParm( argIn );
  if (parm == 0) return Command::C_ERR;
  double scale_factor = argIn.getNextDouble(1.0);
  mprintf("\tScaling dihedral force constants in %s by %f\n", parm->c_str(), scale_factor);
  parm->ScaleDihedralK( scale_factor );
  return Command::C_OK;
}

// ---------- DISPATCHABLE COMMANDS --------------------------------------------
/// Add an action to the State ActionList
Command::RetType AddAction(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return ( (Command::RetType)State.AddAction( Alloc, argIn ) );
}

/// Add an action to the State AnalysisList
Command::RetType AddAnalysis(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return ( (Command::RetType)State.AddAnalysis( Alloc, argIn ) );
}

// -----------------------------------------------------------------------------
/// Warn about deprecated commands.
Command::RetType Deprecated(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  mprintf("Warning: %s is deprecated.\n", argIn.Command());
  return Command::C_OK;
}

// ================ LIST OF ALL COMMANDS =======================================
/** Ideally keep this array first sorted by type (1st field), then 
  * alphabetically by command string (2nd field).
  */
const Command::Token Command::Commands[] = {
  // GENERAL COMMANDS
  { GENERAL, "activeref",     0, Help_ActiveRef,       ActiveRef       },
  { GENERAL, "clear",         0, Help_Clear,           ClearList       },
  { GENERAL, "crdaction",     0, Help_CrdAction,       CrdAction       },
  { GENERAL, "crdout",        0, Help_CrdOut,          CrdOut          },
  { GENERAL, "create",        0, Help_Create_DataFile, Create_DataFile },
  { GENERAL, "datafile",      0, Help_DataFile,        DataFileCmd     },
  { GENERAL, "debug",         0, Help_Debug,           SetListDebug    },
  { GENERAL, "exit" ,         0, Help_Quit,            Quit            },
  { GENERAL, "gnuplot",       0, Help_System,          SystemCmd       },
  { GENERAL, "go",            0, Help_Run,             RunState        },
  { GENERAL, "head",          0, Help_System,          SystemCmd       },
  { GENERAL, "help",          0, Help_Help,            Help            },
  { GENERAL, "list",          0, Help_List,            ListAll         },
  { GENERAL, "loadcrd",       0, Help_LoadCrd,         LoadCrd         },
  { GENERAL, "ls",            0, Help_System,          SystemCmd       },
  { GENERAL, "noexitonerror", 0, Help_NoExitOnError,   NoExitOnError   },
  { GENERAL, "noprogress",    0, Help_NoProgress,      NoProgress      },
  { GENERAL, "precision",     0, Help_Precision,       Precision       },
  { GENERAL, "prnlev",        0, Help_Debug,           SetListDebug    },
  { GENERAL, "pwd",           0, Help_System,          SystemCmd       },
  { GENERAL, "quit" ,         0, Help_Quit,            Quit            },
  { GENERAL, "readdata",      0, Help_ReadData,        ReadData        },
  { GENERAL, "readinput",     0, Help_ReadInput,       ReadInput       },
  { GENERAL, "run"   ,        0, Help_Run,             RunState        },
  { GENERAL, "runanalysis",   0, Help_RunAnalysis,     RunAnalysis     },
  { GENERAL, "select",        0, Help_Select,          SelectAtoms     },
  { GENERAL, "selectds",      0, Help_SelectDS,        SelectDataSets  },
  { GENERAL, "write",         0, Help_Write_DataFile,  Write_DataFile  },
  { GENERAL, "writedata",     0, Help_WriteData,       WriteAllData    },
  { GENERAL, "xmgrace",       0, Help_System,          SystemCmd       },
  // TRAJECTORY COMMANDS
  { TRAJ,    "ensemble",      0, Help_Ensemble,        Trajin          },
  { TRAJ,    "reference",     0, Help_Reference,       Reference       },
  { TRAJ,    "trajin",        0, Help_Trajin,          Trajin          },
  { TRAJ,    "trajout",       0, Help_Trajout,         Trajout         },
  // TOPOLOGY COMMANDS
  { PARM,    "bondinfo",      0, Help_BondInfo,        BondInfo        },
  { PARM,    "charge",        0, Help_ChargeInfo,      ChargeInfo      },
  { PARM,    "molinfo",       0, Help_MolInfo,         MolInfo         },
  { PARM,    "parm",          0, Help_Parm,            LoadParm        },
  { PARM,    "parmbondinfo",  0, Help_BondInfo,        BondInfo        },
  { PARM,    "parmbox",       0, Help_ParmBox,         ParmBox         },
  { PARM,    "parminfo",      0, Help_ParmInfo,        ParmInfo        },
  { PARM,    "parmmolinfo",   0, Help_MolInfo,         MolInfo         },
  { PARM,    "parmresinfo",   0, Help_ResInfo,         ResInfo         },
  { PARM,    "parmstrip",     0, Help_ParmStrip,       ParmStrip       },
  { PARM,    "parmwrite",     0, Help_ParmWrite,       ParmWrite       },
  { PARM,    "resinfo",       0, Help_ResInfo,         ResInfo         },
  { PARM,    "scaledihedralk",0, 0,                    ScaleDihedralK  },
  { PARM,    "solvent",       0, Help_Solvent,         ParmSolvent     },
  // INC_ACTION: ACTION COMMANDS
  { ACTION, "angle", Action_Angle::Alloc, Action_Angle::Help, AddAction },
  { ACTION, "atomiccorr", Action_AtomicCorr::Alloc, Action_AtomicCorr::Help, AddAction },
  { ACTION, "atomicfluct", Action_AtomicFluct::Alloc, Action_AtomicFluct::Help, AddAction },
  { ACTION, "atommap", Action_AtomMap::Alloc, Action_AtomMap::Help, AddAction },
  { ACTION, "autoimage", Action_AutoImage::Alloc, Action_AutoImage::Help, AddAction },
  { ACTION, "average", Action_Average::Alloc, Action_Average::Help, AddAction },
  { ACTION, "bounds", Action_Bounds::Alloc, Action_Bounds::Help, AddAction },
  { ACTION, "box", Action_Box::Alloc, Action_Box::Help, AddAction },
  { ACTION, "center", Action_Center::Alloc, Action_Center::Help, AddAction },
  { ACTION, "check", Action_CheckStructure::Alloc, Action_CheckStructure::Help, AddAction },
  { ACTION, "checkoverlap", Action_CheckStructure::Alloc, Action_CheckStructure::Help, AddAction },
  { ACTION, "checkstructure", Action_CheckStructure::Alloc, Action_CheckStructure::Help, AddAction },
  { ACTION, "closest", Action_Closest::Alloc, Action_Closest::Help, AddAction },
  { ACTION, "closestwaters", Action_Closest::Alloc, Action_Closest::Help, AddAction },
  { ACTION, "clusterdihedral", Action_ClusterDihedral::Alloc, Action_ClusterDihedral::Help, AddAction },
  { ACTION, "contacts", Action_Contacts::Alloc, Action_Contacts::Help, AddAction },
  { ACTION, "createcrd", Action_CreateCrd::Alloc, Action_CreateCrd::Help, AddAction },
  { ACTION, "createreservoir", Action_CreateReservoir::Alloc, Action_CreateReservoir::Help, AddAction },
  { ACTION, "density", Action_Density::Alloc, Action_Density::Help, AddAction },
  { ACTION, "diffusion", Action_Diffusion::Alloc, Action_Diffusion::Help, AddAction },
  { ACTION, "dihedral", Action_Dihedral::Alloc, Action_Dihedral::Help, AddAction },
  { ACTION, "dihedralscan", Action_DihedralScan::Alloc, Action_DihedralScan::Help, AddAction },
  { ACTION, "dipole", Action_Dipole::Alloc, Action_Dipole::Help, AddAction },
  { ACTION, "distance", Action_Distance::Alloc, Action_Distance::Help, AddAction },
//  { ACTION, "dnaiontracker", Action_DNAionTracker::Alloc, Action_DNAionTracker::Help, AddAction },
  { ACTION, "drms", Action_DistRmsd::Alloc, Action_DistRmsd::Help, AddAction },
  { ACTION, "drmsd", Action_DistRmsd::Alloc, Action_DistRmsd::Help, AddAction },
  { ACTION, "fixatomorder", Action_FixAtomOrder::Alloc, Action_FixAtomOrder::Help, AddAction },
//  { ACTION, "gfe", Action_GridFreeEnergy::Alloc, Action_GridFreeEnergy::Help, AddAction },
  { ACTION, "grid", Action_Grid::Alloc, Action_Grid::Help, AddAction },
  { ACTION, "hbond", Action_Hbond::Alloc, Action_Hbond::Help, AddAction },
  { ACTION, "image", Action_Image::Alloc, Action_Image::Help, AddAction },
  { ACTION, "jcoupling", Action_Jcoupling::Alloc, Action_Jcoupling::Help, AddAction },
  { ACTION, "lie", Action_LIE::Alloc, Action_LIE::Help, AddAction },
  { ACTION, "makestructure", Action_MakeStructure::Alloc, Action_MakeStructure::Help, AddAction },
  { ACTION, "mask", Action_Mask::Alloc, Action_Mask::Help, AddAction },
  { ACTION, "matrix", Action_Matrix::Alloc, Action_Matrix::Help, AddAction },
  { ACTION, "mindist", Action_MinDist::Alloc, Action_MinDist::Help, AddAction },
  { ACTION, "molsurf", Action_Molsurf::Alloc, Action_Molsurf::Help, AddAction },
  { ACTION, "multidihedral", Action_MultiDihedral::Alloc, Action_MultiDihedral::Help, AddAction },
  { ACTION, "nastruct", Action_NAstruct::Alloc, Action_NAstruct::Help, AddAction },
  { ACTION, "lipidorder", Action_OrderParameter::Alloc, Action_OrderParameter::Help, AddAction },
  { ACTION, "outtraj", Action_Outtraj::Alloc, Action_Outtraj::Help, AddAction },
  { ACTION, "pairdist", Action_PairDist::Alloc, Action_PairDist::Help, AddAction },
  { ACTION, "pairwise", Action_Pairwise::Alloc, Action_Pairwise::Help, AddAction },
  { ACTION, "principal", Action_Principal::Alloc, Action_Principal::Help, AddAction },
  { ACTION, "projection", Action_Projection::Alloc, Action_Projection::Help, AddAction },
  { ACTION, "pucker", Action_Pucker::Alloc, Action_Pucker::Help, AddAction },
  { ACTION, "radgyr", Action_Radgyr::Alloc, Action_Radgyr::Help, AddAction },
  { ACTION, "radial", Action_Radial::Alloc, Action_Radial::Help, AddAction },
  { ACTION, "randomizeions", Action_RandomizeIons::Alloc, Action_RandomizeIons::Help, AddAction },
  { ACTION, "rms", Action_Rmsd::Alloc, Action_Rmsd::Help, AddAction },
  { ACTION, "rmsd", Action_Rmsd::Alloc, Action_Rmsd::Help, AddAction },
  { ACTION, "rog", Action_Radgyr::Alloc, Action_Radgyr::Help, AddAction },
  { ACTION, "rotate", Action_Rotate::Alloc, Action_Rotate::Help, AddAction },
  { ACTION, "rotdif", Action_Rotdif::Alloc, Action_Rotdif::Help, AddAction },
  { ACTION, "runavg", Action_RunningAvg::Alloc, Action_RunningAvg::Help, AddAction },
  { ACTION, "runningaverage", Action_RunningAvg::Alloc, Action_RunningAvg::Help, AddAction },
  { ACTION, "scale", Action_Scale::Alloc, Action_Scale::Help, AddAction },
  { ACTION, "secstruct", Action_DSSP::Alloc, Action_DSSP::Help, AddAction },
  { ACTION, "spam", Action_Spam::Alloc, Action_Spam::Help, AddAction },
  { ACTION, "stfcdiffusion", Action_STFC_Diffusion::Alloc, Action_STFC_Diffusion::Help, AddAction },
  { ACTION, "strip", Action_Strip::Alloc, Action_Strip::Help, AddAction },
  { ACTION, "surf", Action_Surf::Alloc, Action_Surf::Help, AddAction },
  { ACTION, "symmrmsd", Action_SymmetricRmsd::Alloc, Action_SymmetricRmsd::Help, AddAction },
  { ACTION, "temperature", Action_Temperature::Alloc, Action_Temperature::Help, AddAction },
  { ACTION, "trans", Action_Translate::Alloc, Action_Translate::Help, AddAction },
  { ACTION, "translate", Action_Translate::Alloc, Action_Translate::Help, AddAction },
  { ACTION, "unstrip", Action_Unstrip::Alloc, Action_Unstrip::Help, AddAction },
  { ACTION, "unwrap", Action_Unwrap::Alloc, Action_Unwrap::Help, AddAction },
  { ACTION, "vector", Action_Vector::Alloc, Action_Vector::Help, AddAction },
  { ACTION, "watershell", Action_Watershell::Alloc, Action_Watershell::Help, AddAction },
  { ACTION, "volmap", Action_Volmap::Alloc, Action_Volmap::Help, AddAction},
  // INC_ANALYSIS: ANALYSIS COMMANDS
  { ANALYSIS, "2drms", Analysis_Rms2d::Alloc, Analysis_Rms2d::Help, AddAnalysis },
  { ANALYSIS, "amdbias", Analysis_AmdBias::Alloc, Analysis_AmdBias::Help, AddAnalysis },
  { ANALYSIS, "autocorr", Analysis_AutoCorr::Alloc, Analysis_AutoCorr::Help, AddAnalysis },
  { ANALYSIS, "cluster", Analysis_Clustering::Alloc, Analysis_Clustering::Help, AddAnalysis },
  { ANALYSIS, "corr", Analysis_Corr::Alloc, Analysis_Corr::Help, AddAnalysis },
  { ANALYSIS, "correlationcoe", Analysis_Corr::Alloc, Analysis_Corr::Help, AddAnalysis },
  { ANALYSIS, "crank", Analysis_CrankShaft::Alloc, Analysis_CrankShaft::Help, AddAnalysis },
  { ANALYSIS, "crankshaft", Analysis_CrankShaft::Alloc, Analysis_CrankShaft::Help, AddAnalysis },
  { ANALYSIS, "crdfluct", Analysis_CrdFluct::Alloc, Analysis_CrdFluct::Help, AddAnalysis },
  { ANALYSIS, "crosscorr", Analysis_CrossCorr::Alloc, Analysis_CrossCorr::Help, AddAnalysis },
  { ANALYSIS, "diagmatrix", Analysis_Matrix::Alloc, Analysis_Matrix::Help, AddAnalysis },
  { ANALYSIS, "fft", Analysis_FFT::Alloc, Analysis_FFT::Help, AddAnalysis },
  { ANALYSIS, "hist", Analysis_Hist::Alloc, Analysis_Hist::Help, AddAnalysis },
  { ANALYSIS, "histogram", Analysis_Hist::Alloc, Analysis_Hist::Help, AddAnalysis },
  { ANALYSIS, "ired", Analysis_IRED::Alloc, Analysis_IRED::Help, AddAnalysis },
  { ANALYSIS, "lifetime", Analysis_Lifetime::Alloc, Analysis_Lifetime::Help, AddAnalysis },
  { ANALYSIS, "matrix", Analysis_Matrix::Alloc, Analysis_Matrix::Help, AddAnalysis },
  { ANALYSIS, "meltcurve", Analysis_MeltCurve::Alloc, Analysis_MeltCurve::Help, AddAnalysis },
  { ANALYSIS, "modes", Analysis_Modes::Alloc, Analysis_Modes::Help, AddAnalysis },
  { ANALYSIS, "overlap", Analysis_Overlap::Alloc, Analysis_Overlap::Help, AddAnalysis },
  { ANALYSIS, "remlog", Analysis_RemLog::Alloc, Analysis_RemLog::Help, AddAnalysis },
  { ANALYSIS, "rms2d", Analysis_Rms2d::Alloc, Analysis_Rms2d::Help, AddAnalysis },
  { ANALYSIS, "rmsavgcorr", Analysis_RmsAvgCorr::Alloc, Analysis_RmsAvgCorr::Help, AddAnalysis },
  { ANALYSIS, "stat", Analysis_Statistics::Alloc, Analysis_Statistics::Help, AddAnalysis },
  { ANALYSIS, "statistics", Analysis_Statistics::Alloc, Analysis_Statistics::Help, AddAnalysis },
  { ANALYSIS, "timecorr", Analysis_Timecorr::Alloc, Analysis_Timecorr::Help, AddAnalysis },
  { ANALYSIS, "runningavg", Analysis_RunningAvg::Alloc, Analysis_RunningAvg::Help, AddAnalysis },
  // DEPRECATED COMMANDS
  { DEPRECATED, "acceptor",     0, 0, Deprecated },
  { DEPRECATED, "bondsearch",   0, 0, Deprecated },
  { DEPRECATED, "donor",        0, 0, Deprecated },
  { DEPRECATED, "molsearch",    0, 0, Deprecated },
  { DEPRECATED, "nobondsearch", 0, 0, Deprecated },
  { DEPRECATED, "nomolsearch",  0, 0, Deprecated },
  { NONE      , 0,              0, 0, 0          }
};
