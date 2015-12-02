#include <cstdlib> // system
#include "CmdInput.h" // ProcessInput()
#include "Command.h"
#include "CpptrajStdio.h"
#include "ParmFile.h" // ReadOptions, WriteOptions
#include "Timer.h"
#include "RPNcalc.h" // Calc
#include "Cmd_CompareTop.h"
#include "Cmd_Coords.h"
#include "Cmd_GenerateAmberRst.h"
#include "Cmd_SequenceAlign.h"
#include "ProgressBar.h"
// INC_ACTION==================== ALL ACTION CLASSES GO HERE ===================
#include "Action_Angle.h"
#include "Action_Distance.h"
#include "Action_Rmsd.h"
#include "Action_Dihedral.h"
#include "Action_AtomMap.h"
#include "Action_Strip.h"
#include "Action_Unstrip.h"
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
#include "Action_Gist.h"
#include "Action_CreateReservoir.h"
#include "Action_Density.h"
#include "Action_PairDist.h"
#include "Action_OrderParameter.h"
#include "Action_FixAtomOrder.h"
#include "Action_NMRrst.h"
#include "Action_FilterByData.h"
#include "Action_LESsplit.h"
#include "Action_NativeContacts.h"
#include "Action_VelocityAutoCorr.h"
#include "Action_SetVelocity.h"
#include "Action_MultiVector.h"
#include "Action_MinImage.h"
#include "Action_ReplicateCell.h"
#include "Action_AreaPerMol.h"
#include "Action_Energy.h"
#include "Action_CheckChirality.h"
#include "Action_Channel.h" // EXPERIMENTAL
#include "Action_Volume.h"

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
#include "Analysis_Integrate.h"
#include "Analysis_Spline.h"
#include "Analysis_Average.h"
#include "Analysis_KDE.h"
#include "Analysis_MultiHist.h"
#include "Analysis_Divergence.h"
#include "Analysis_VectorMath.h"
#include "Analysis_Regression.h"
#include "Analysis_LowestCurve.h"
#include "Analysis_CurveFit.h"
#include "Analysis_PhiPsi.h"
#include "Analysis_Rotdif.h"
#include "Analysis_Wavelet.h"
#include "Analysis_State.h"
#include "Analysis_Multicurve.h"
#include "Analysis_TI.h"
// ---- Command Functions ------------------------------------------------------
/// Warn about deprecated commands.
void Command::WarnDeprecated(Cmd::TokenPtr token)
{
  mprinterr("Error: '%s' is deprecated.\n", token->Cmd);
  if (token->Help != 0)
    token->Help();
}

/** Search Commands list for a specific type of command. */
Cmd::TokenPtr Command::SearchTokenType(Cmd::Ctype dtype,
                                           ArgList const& argIn)
{
  for (Cmd::TokenPtr token = Commands; token->Type != Cmd::NONE; ++token)
  {
    if (token->Type == Cmd::DEPRECATED && argIn.CommandIs( token->Cmd )) {
      WarnDeprecated( token );
      return 0;
    }
    if (dtype != token->Type) continue;
    if (argIn.CommandIs( token->Cmd )) return token;
  }
  mprinterr("'%s': Command not found.\n", argIn.Command());
  return 0;
}

/// Strings that correspond to Cmd::Ctype
const char* Command::CommandTitle[] = { 0, "Topology", "Trajectory", "Coords",
  "Action", "Analysis", "General", "System", "Deprecated" };

/** List all commands of the given type, or all commands if type
  * is Cmd::NONE.
  */
void Command::ListCommands(Cmd::Ctype dtype) {
  std::string Line;
  Cmd::Ctype lastType = Cmd::NONE;
  for (Cmd::TokenPtr token = Commands; token->Type != Cmd::DEPRECATED; ++token)
  {
    Cmd::Ctype currentType = token->Type;
    if (currentType == Cmd::HIDDEN) continue;
    if (dtype != Cmd::NONE && dtype != currentType) continue;
    // Command group type title
    if (currentType != lastType) {
      if (!Line.empty()) {
        mprintf("%s\n", Line.c_str());
        Line.clear();
      }
      mprintf("%s Commands:\n", CommandTitle[currentType]);
      lastType = currentType;
    }
    if (Line.empty()) Line.assign("        ");
    std::string Command(token->Cmd);
    Command.append(" ");
    if ( Line.size() + Command.size() > 80 ) {
      mprintf("%s\n", Line.c_str());
      Line.assign("        ");
    }
    Line.append(Command);
  }
  if (!Line.empty())
    mprintf("%s\n", Line.c_str());
}

/** Search the Commands list for given command.
  * \return the token if found, 0 if not.
  */
Cmd::TokenPtr Command::SearchToken(ArgList& argIn) {
  // SPECIAL CASE: For backwards compat. remove analyze prefix
  if (argIn.CommandIs("analyze")) {
    argIn.RemoveFirstArg();
    argIn.MarkArg(0); // Mark new first arg as command
    return (SearchTokenType(Cmd::ANALYSIS, argIn));
  }
  // Search for command.
  for (Cmd::TokenPtr token = Commands; token->Type != Cmd::NONE; ++token)
    if (argIn.CommandIs( token->Cmd )) {
      if (token->Type == Cmd::DEPRECATED) {
        WarnDeprecated( token );
        return 0; 
      } else
        return token;
    }
  //mprinterr("'%s': Command not found.\n", argIn.Command());
  return 0;
}

/** Search for the given command and execute it. */
Cmd::RetType Command::Dispatch(CpptrajState& State,
                                   std::string const& commandIn)
{
  ArgList cmdArg( commandIn );
  cmdArg.MarkArg(0); // Always mark the first arg as the command 
  Cmd::TokenPtr cmdToken = SearchToken( cmdArg );
  Cmd::RetType ret_val = Cmd::OK;
  if (cmdToken == 0) {
    // Try to evaluate the expression.
    RPNcalc calc;
    calc.SetDebug( State.Debug() );
    if (calc.ProcessExpression( commandIn ))
      ret_val = Cmd::ERR;
    else {
      if (calc.Evaluate(*State.DSL()))
        ret_val = Cmd::ERR;
    }
    if (ret_val == Cmd::ERR)
      mprinterr("'%s': Invalid command or expression.\n", commandIn.c_str());
  } else
    ret_val = cmdToken->Fxn( State, cmdArg, cmdToken->Alloc );
  return ret_val;
}

/** Read commands from given input file.
  * \return 0 if successfully read, 1 on error.
  */
Cmd::RetType Command::ProcessInput(CpptrajState& State, 
                                       std::string const& inputFilename)
{
  BufferedLine infile;
  if (infile.OpenFileRead( inputFilename )) {
    if (!inputFilename.empty())
      mprinterr("Error: Could not open input file '%s'\n", inputFilename.c_str());
    return Cmd::ERR;
  }
  mprintf("INPUT: Reading input from '%s'\n", infile.Filename().full());
  // Read in each line of input.
  int nInputErrors = 0;
  Cmd::RetType cmode = Cmd::OK;
  CmdInput input;
  const char* ptr = infile.Line();
  while (ptr != 0) {
    bool moreInput = input.AddInput( ptr );
    while (moreInput) {
      ptr = infile.Line();
      moreInput = input.AddInput( ptr );
    }
    // Only attempt to execute if the command is not blank.
    if (!input.Empty()) {
      // Print the input line that will be sent to dispatch
      mprintf("  [%s]\n", input.str());
      // Call Dispatch to convert input to ArgList and process.
      cmode = Command::Dispatch(State, input.Str());
      if (cmode == Cmd::ERR) {
        nInputErrors++;
        if (State.ExitOnError()) break;
      } else if (cmode == Cmd::QUIT)
        break;
    }
    // Reset Input line
    input.Clear();
    ptr = infile.Line();
  }
  infile.CloseFile();
  if (nInputErrors > 0) {
    mprinterr("\t%i errors encountered reading input.\n", nInputErrors);
    return Cmd::ERR;
  }
  return cmode;
}

// ====================== CPPTRAJ COMMANDS HELP ================================
static void Help_System() { mprintf("  Call command from system.\n"); }

static void Help_NoProgress() {
  mprintf("  Do not print progress while reading in trajectories.\n");
}

static void Help_NoExitOnError() {
  mprintf("  Do not exit when errors are encountered. This is the default\n"
          "  in interactive mode.\n");
}

static void Help_Run() {
  mprintf("  Process all trajectories currently in input trajectory list.\n"
          "  All actions in action list will be run on each frame.\n"
          "  If not processing ensemble input, all analyses in analysis\n"
          "  list will be run after trajectory processing.\n");
}

static void Help_Quit() { mprintf("  Exit CPPTRAJ\n"); }

static void Help_List() {
  mprintf("\t[<type>] (<type> =%s)\n"
          "  List currently loaded objects of the specified type. If no type is given\n"
          "  then list all loaded objects.\n", CpptrajState::PrintListKeys().c_str());
}

static void Help_Debug() {
  mprintf("\t[<type>] <#> (<type> =%s)\n", CpptrajState::PrintListKeys().c_str());
  mprintf("  Set debug level for new objects of the specified type. If no type is given\n"
          "  then set debug level for all new objects. Does not affect current objects.\n");
}

static void Help_Clear() {
  mprintf("\t[ {all | <type>} ] (<type> =%s)\n", CpptrajState::PrintListKeys().c_str());
  mprintf("  Clear currently loaded objects of the specified type. If 'all' is specified\n"
          "  then clear all loaded objects.\n");
}

static void Help_RemoveData() {
  mprintf("\t[<arg>]\n"
          "  Remove data sets(s) corresponding to <arg> from data set list.\n");
}

static void Help_Create_DataFile() {
  mprintf("\t<filename> <dataset0> [<dataset1> ...]\n"
          "  Add a file with specified data sets to the data file list. Does not\n"
          "  immediately write the data.\n");
  DataFile::WriteHelp();
  DataFile::WriteOptions();
}

static void Help_DataFile() {
  mprintf("\t<data filename> <datafile cmd>\n"
          "  Pass <datafile cmd> to specified data file currently in data file list.\n");
  DataFile::WriteHelp();
  DataFile::WriteOptions();
}

static void Help_ReadData() {
  mprintf("\t<filename> [name <dsname>] [as <fmt>] [<format options>]\n"
          "  Read data from <filename> into data sets.\n");
  DataFile::ReadOptions();
}

static void Help_ReadInput() {
  mprintf("\t<filename>\n"
          "  Read commands from input file <filename>\n");
}

static void Help_Write_DataFile() {
  mprintf("\t[<filename> <dataset0> [<dataset1> ...]]\n");
  DataFile::WriteHelp();
  mprintf("  With no arguments, write all files currently in the data file list.\n"
          "  Otherwise, write specified data sets to <filename> immediately.\n");
  DataFile::WriteOptions();
}

static void Help_Precision() {
  mprintf("\t{<filename> | <dataset arg>} [<width>] [<precision>]\n"
          "  Set precision for all datasets in datafile <filename> or dataset(s)\n"
          "  specified by <dataset arg> to <width>.<precision>. If width/precision\n"
          "  is not specified then default to 12.4\n");
}

static void Help_Select() {
  mprintf("\t[<parmindex>] <mask>\n"
          "  Show atom numbers selected by <mask> for parm <parmindex>\n"
          "  (default first parm)\n");
}

static void Help_SelectDS() {
  mprintf("\t<dataset selection>\n"
          "  Show results of data set selection. Data set selection format is:\n"
          "\t<name>[<aspect]:<idx range>\n"
          "  Where '<name>' is the data set name, '[<aspect>]' is the data set aspect,\n"
          "  and <idx range> is a numerical range specifying data set indices (i.e. 2-5,7 etc).\n"
          "  The aspect and index portions may be optional. An asterisk '*' may be used as\n"
          "  a wildcard. E.g. 'selectds R2', 'selectds RoG[Max]', 'selectds PR[res]:2-12'\n");
}

static void Help_Trajin() {
  mprintf("\t<filename> {[<start>] [<stop> | last] [offset]} | lastframe\n"
          "\t           [%s]\n", DataSetList::TopArgs);
  mprintf("\t           [ <Format Options> ]\n"
          "\t           [ remdtraj [remdtrajtemp <T> | remdtrajidx <#>]\n"
          "\t             [trajnames <rep1>,<rep2>,...,<repN> ] ]\n"
          "  Load trajectory specified by <filename> to the input trajectory list.\n");
  TrajectoryFile::ReadOptions();
}

static void Help_Ensemble() {
  mprintf("\t<file0> {[<start>] [<stop> | last] [offset]} | lastframe\n"
          "\t        [%s]\n", DataSetList::TopArgs);
  mprintf("\t        [trajnames <file1>,<file2>,...,<fileN>\n"
          "\t        [remlog <remlogfile> [nstlim <nstlim> ntwx <ntwx>]]\n"
          "  Load an ensemble of trajectories starting with <file0> that will be\n"
          "  processed together as an ensemble.\n");
}

static void Help_Trajout() {
  mprintf("\t<filename> [<fileformat>] [append] [nobox]\n"
          "\t           [%s] [onlyframes <range>] [title <title>]\n", DataSetList::TopArgs);
  mprintf("\t           %s\n", ActionFrameCounter::HelpText);
  mprintf("\t           [ <Format Options> ]\n"
          "  Write frames after all actions have been processed to output trajectory\n"
          "  specified by <filename>.\n");
  TrajectoryFile::WriteOptions();
}

static void Help_Reference() {
  mprintf("\t<name> [<frame#>] [<mask>] [TAG] [lastframe] [crdset]\n"
          "\t       [%s]\n", DataSetList::TopArgs);
  mprintf("  Load trajectory file <name> as a reference frame.\n"
          "  If 'crdset' is specified use COORDS data set specified by <name> as reference.\n");
}

static void Help_RunAnalysis() {
  mprintf("\t[<analysis> [<analysis args>]]\n"
          "  If specified alone, run all analyses in the analysis list.\n"
          "  Otherwise run the specified analysis immediately.\n");
}

// ---------- Information on Deprecated commands -------------------------------
static void Deprecate_MinDist() {
  mprinterr("  Use the 'nativecontacts' action instead.\n");
}

static void Deprecate_Hbond() {
  mprinterr("  Hydrogen bond acceptors and donors are defined within the 'hbond' action.\n");
}

static void Deprecate_TopSearch() {
  mprinterr("  Bonds and/or molecules are automatically searched for if needed.\n");
}

static void Deprecate_ParmBondInfo() {
  mprinterr("  Use bonds, bondinfo, or printbonds instead.\n");
}

static void Deprecate_ParmResInfo() {
  mprinterr("  Use resinfo instead.\n");
}

static void Deprecate_ParmMolInfo() {
  mprinterr("  Use molinfo instead.\n");
}

static void Deprecate_AvgCoord() {
  mprinterr("  Use 'vector center' (optionally with keyword 'magnitude') instead.\n");
}

// ---------- GENERAL COMMANDS -------------------------------------------------
static void Help_ActiveRef() {
  mprintf("\t%s\n", DataSetList::RefArgs);
  mprintf("  Set the reference structure to be used for coordinate-based mask parsing.\n"
          "  <#> starts from 0 (first reference structure).\n");
}

/// Set active reference for distance-based masks etc.
Cmd::RetType ActiveRef(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  return (Cmd::RetType)State.DSL()->SetActiveReference( argIn );
}

/// Clear data in specified lists
Cmd::RetType ClearList(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  return (Cmd::RetType)State.ClearList( argIn );
}

Cmd::RetType RemoveData(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  return (Cmd::RetType)State.RemoveDataSet( argIn );
}

/// Set debug value for specified list(s)
Cmd::RetType SetListDebug(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  return (Cmd::RetType)State.SetListDebug( argIn );
}

/// List all members of specified list(s)
Cmd::RetType ListAll(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  return (Cmd::RetType)State.ListAll( argIn );
}

static void Help_SilenceActions() { mprintf("Silence Actions Init/Setup output.\n"); }
/// Silence Actions Init/Setup output.
Cmd::RetType SilenceActions(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{ State.SetActionSilence( true ); return Cmd::OK; }

// ========== COORDS ===========================================================
/// Perform action on given COORDS dataset
Cmd::RetType CrdActionCmd(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  std::string setname = argIn.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: %s: Specify COORDS dataset name.\n", argIn.Command());
    return Cmd::ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL()->FindCoordsSet( setname );
  if (CRD == 0) {
    mprinterr("Error: %s: No COORDS set with name %s found.\n", argIn.Command(), setname.c_str());
    return Cmd::ERR;
  }
  mprintf("\tUsing set '%s'\n", CRD->legend());
  // Start, stop, offset
  TrajFrameCounter frameCount;
  ArgList crdarg( argIn.GetStringKey("crdframes"), "," );
  if (frameCount.CheckFrameArgs( CRD->Size(), crdarg )) return Cmd::ERR;
  frameCount.PrintInfoLine(CRD->legend());
  ArgList actionargs = argIn.RemainingArgs();
  actionargs.MarkArg(0);
  Cmd::TokenPtr tkn = Command::SearchTokenType( Cmd::ACTION, actionargs);
  if ( tkn == 0 ) return Cmd::ERR;
  Action* act = (Action*)tkn->Alloc();
  if (act == 0) return Cmd::ERR;
  return CrdAction(State, actionargs, CRD, act, frameCount);
}

// -----------------------------------------------------------------------------
/// Add DataSets specified by arguments to given DataFile.
// NOTE: Used by Create_DataFile and Write_DataFile
// TODO: Put in DataFile?
static int AddSetsToDataFile(DataFile& df, ArgList const& dsetArgs, DataSetList& DSL)
{
  int err = 0;
  for (ArgList::const_iterator dsa = dsetArgs.begin(); dsa != dsetArgs.end(); ++dsa) {
    DataSetList Sets = DSL.GetMultipleSets( *dsa );
    if (Sets.empty())
      mprintf("Warning: %s does not correspond to any data sets.\n", dsa->c_str());
    for (DataSetList::const_iterator set = Sets.begin(); set != Sets.end(); ++set) {
      mprintf(" %s", (*set)->legend());
      if ( df.AddDataSet(*set) ) {
        mprinterr("Error: Could not add data set %s to file.\n", (*set)->legend());
        ++err;
      }
    }
  }
  mprintf("\n");
  return err;
}

/// Add a new DataFile to DFL with specified DataSets, to be written later.
Cmd::RetType Create_DataFile(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  // Next string is datafile that command pertains to.
  std::string name1 = argIn.GetStringNext();
  if (name1.empty()) {
    mprinterr("Error: No filename given.\n");
    return Cmd::ERR;
  }
  DataFile* df = State.DFL()->AddDataFile(name1, argIn);
  if (df == 0) return Cmd::ERR;
  return (Cmd::RetType)( AddSetsToDataFile(*df, argIn.RemainingArgs(), *(State.DSL())) );
}

/// Write DataFile with specified DataSets immediately, or force write of all DataFiles in State
Cmd::RetType Write_DataFile(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  // Next string is datafile that command pertains to.
  std::string name1 = argIn.GetStringNext();
  if (name1.empty()) {
    State.DFL()->ResetWriteStatus();
    State.MasterDataFileWrite();
    return Cmd::OK;
  }
  DataFile* df = new DataFile();
  if (df == 0) return Cmd::ERR;
  if (df->SetupDatafile( name1, argIn, State.Debug() )) {
    delete df;
    return Cmd::ERR;
  }
  mprintf("\tWriting sets to %s, format '%s'\n", df->DataFilename().full(), df->FormatString());
  int err = AddSetsToDataFile(*df, argIn.RemainingArgs(), *(State.DSL()));
  if (err == 0) df->WriteDataOut();
  delete df;
  return (Cmd::RetType)err;
}

/// Process DataFile-specific command
Cmd::RetType DataFileCmd(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  return (Cmd::RetType)( State.DFL()->ProcessDataFileArgs( argIn ) );
}

// -----------------------------------------------------------------------------
static void Help_DataSetCmd() {
  mprintf("\t{ legend <legend> <set> | makexy <Xset> <Yset> [name <name>] |\n"
          "\t  cat <set0> <set1> ... [name <name>] [nooffset]\n"
          "\t  [mode <mode>] [type <type>] <set arg1> [<set arg 2> ...] }\n"
          "\t<mode>: ");
  for (int i = 0; i != (int)MetaData::M_MATRIX; i++) // TODO: Allow matrix?
    mprintf(" %s", MetaData::ModeString((MetaData::scalarMode)i));
  mprintf("\n\t<type>: ");
  for (int i = 0; i != (int)MetaData::DIST; i++)
    mprintf(" %s", MetaData::TypeString((MetaData::scalarType)i));
  mprintf("\n\tOptions for 'type noe':\n"
          "\t  %s\n"
          "  legend: Set the legend for a single data set\n"
          "  makexy: Create new data set with X values from one set and Y values from another.\n"
          "  cat   : Concatenate 2 or more data sets.\n"
          "  Otherwise, change the mode/type for one or more data sets.\n",
          AssociatedData_NOE::HelpText);
}

/// Process DataSet-specific command
Cmd::RetType DataSetCmd(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  if (argIn.Contains("legend")) { // Set legend for one data set
    std::string legend = argIn.GetStringKey("legend");
    DataSet* ds = State.DSL()->GetDataSet( argIn.GetStringNext() );
    if (ds == 0) return Cmd::ERR;
    mprintf("\tChanging legend '%s' to '%s'\n", ds->legend(), legend.c_str());
    ds->SetLegend( legend );
  // ---------------------------------------------
  } else if (argIn.hasKey("makexy")) { // Combine values from two sets into 1
    std::string name = argIn.GetStringKey("name");
    DataSet* ds1 = State.DSL()->GetDataSet( argIn.GetStringNext() );
    DataSet* ds2 = State.DSL()->GetDataSet( argIn.GetStringNext() );
    if (ds1 == 0 || ds2 == 0) return Cmd::ERR;
    if (ds1->Ndim() != 1 || ds2->Ndim() != 1) {
      mprinterr("Error: makexy only works for 1D data sets.\n");
      return Cmd::ERR;
    }
    DataSet* ds3 = State.DSL()->AddSet( DataSet::XYMESH, name, "XY" );
    if (ds3 == 0) return Cmd::ERR;
    mprintf("\tUsing values from '%s' as X, values from '%s' as Y, output set '%s'\n",
            ds1->legend(), ds2->legend(), ds3->legend());
    DataSet_1D const& ds_x = static_cast<DataSet_1D const&>( *ds1 );
    DataSet_1D const& ds_y = static_cast<DataSet_1D const&>( *ds2 );
    DataSet_1D&       out  = static_cast<DataSet_1D&>( *ds3 );
    size_t nframes = std::min( ds_x.Size(), ds_y.Size() );
    if (ds_x.Size() != ds_y.Size())
      mprintf("Warning: Data sets do not have equal sizes, only using %zu frames.\n", nframes);
    double XY[2];
    for (size_t i = 0; i != nframes; i++) {
      XY[0] = ds_x.Dval(i);
      XY[1] = ds_y.Dval(i);
      out.Add( i, XY );
    }
  // ---------------------------------------------
  } else if (argIn.hasKey("cat")) { // Concatenate two or more data sets
    std::string name = argIn.GetStringKey("name");
    bool use_offset = !argIn.hasKey("nooffset");
    DataSet* ds3 = State.DSL()->AddSet( DataSet::XYMESH, name, "CAT" );
    if (ds3 == 0) return Cmd::ERR;
    DataSet_1D& out = static_cast<DataSet_1D&>( *ds3 );
    mprintf("\tConcatenating sets into '%s'\n", out.legend());
    if (use_offset)
      mprintf("\tX values will be offset.\n");
    else
      mprintf("\tX values will not be offset.\n");
    std::string dsarg = argIn.GetStringNext();
    double offset = 0.0;
    while (!dsarg.empty()) {
      DataSetList dsl = State.DSL()->GetMultipleSets( dsarg );
      double XY[2];
      for (DataSetList::const_iterator ds = dsl.begin(); ds != dsl.end(); ++ds)
      {
        if ( (*ds)->Type() != DataSet::INTEGER &&
             (*ds)->Type() != DataSet::DOUBLE &&
             (*ds)->Type() != DataSet::FLOAT &&
             (*ds)->Type() != DataSet::XYMESH )
        {
          mprintf("Warning: '%s': Concatenation only supported for 1D scalar data sets.\n",
                  (*ds)->legend());
        } else {
          DataSet_1D const& set = static_cast<DataSet_1D const&>( *(*ds) );
          mprintf("\t\t'%s'\n", set.legend());
          for (size_t i = 0; i != set.Size(); i++) {
            XY[0] = set.Xcrd( i ) + offset;
            XY[1] = set.Dval( i );
            out.Add( i, XY ); // NOTE: value of i does not matter for mesh
          }
          if (use_offset) offset = XY[0];
        }
      }
      dsarg = argIn.GetStringNext();
    }
  // ---------------------------------------------
  } else { // Default: change mode/type for one or more sets.
    std::string modeKey = argIn.GetStringKey("mode");
    std::string typeKey = argIn.GetStringKey("type");
    if (modeKey.empty() && typeKey.empty()) {
      mprinterr("Error: No valid keywords specified.\n");
      return Cmd::ERR;
    }
    // First determine mode if specified.
    MetaData::scalarMode dmode = MetaData::UNKNOWN_MODE;
    if (!modeKey.empty()) {
      dmode = MetaData::ModeFromKeyword( modeKey );
      if (dmode == MetaData::UNKNOWN_MODE) {
        mprinterr("Error: Invalid mode keyword '%s'\n", modeKey.c_str());
        return Cmd::ERR;
      }
    }
    // Next determine type if specified.
    MetaData::scalarType dtype = MetaData::UNDEFINED;
    if (!typeKey.empty()) {
      dtype = MetaData::TypeFromKeyword( typeKey, dmode );
      if (dtype == MetaData::UNDEFINED) {
        mprinterr("Error: Invalid type keyword '%s'\n", typeKey.c_str());
        return Cmd::ERR;
      }
    }
    // Additional options for type 'noe'
    AssociatedData_NOE noeData;
    if (dtype == MetaData::NOE) {
      if (noeData.NOE_Args(argIn))
        return Cmd::ERR;
    }
    if (dmode != MetaData::UNKNOWN_MODE)
      mprintf("\tDataSet mode = %s\n", MetaData::ModeString(dmode));
    if (dtype != MetaData::UNDEFINED)
      mprintf("\tDataSet type = %s\n", MetaData::TypeString(dtype));
    // Loop over all DataSet arguments 
    std::string ds_arg = argIn.GetStringNext();
    while (!ds_arg.empty()) {
      DataSetList dsl = State.DSL()->GetMultipleSets( ds_arg );
      for (DataSetList::const_iterator ds = dsl.begin(); ds != dsl.end(); ++ds)
      {
        if ( (*ds)->Ndim() != 1 ) // TODO remove restriction
          mprintf("Warning:\t\t'%s': Can only set mode/type for 1D data sets.\n",
                  (*ds)->legend());
        else {
          if ( dtype == MetaData::NOE ) (*ds)->AssociateData( &noeData );
          mprintf("\t\t'%s'\n", (*ds)->legend());
          MetaData md = (*ds)->Meta();
          md.SetScalarMode( dmode );
          md.SetScalarType( dtype );
          (*ds)->SetMeta( md );
        }
      }
      ds_arg = argIn.GetStringNext();
    }
  }
  return Cmd::OK;
}

// -----------------------------------------------------------------------------
static void Help_DataFilter() {
  mprintf("\t<dataset arg> min <min> max <max> [out <file> [name <setname>]]\n"
          "  Create a data set (optionally named <setname>) containing 1 for\n"
          "  data within given <min> and <max> criteria for each specified\n"
          "  data set. There must be at least one <min> and <max> argument,\n"
          "  and can be as many as there are specified data sets.\n");
}

/// Use the filter command on DataSets outside trajectory processing.
Cmd::RetType DataFilter(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  Action_FilterByData filterAction;
  ActionInit state(*State.DSL(), *State.DFL());
  if (filterAction.Init(argIn, state, State.Debug()) != Action::OK)
    return Cmd::ERR;
  size_t nframes = filterAction.DetermineFrames();
  if (nframes < 1) {
    mprinterr("Error: No data to filter. All sets must contain some data.\n");
    return Cmd::ERR;
  }
  ProgressBar progress( nframes );
  ActionFrame frm;
  for (size_t frame = 0; frame != nframes; frame++) {
    progress.Update( frame );
    filterAction.DoAction(frame, frm); // Filter does not need frame.
  }
  // Trigger master datafile write just in case
  State.MasterDataFileWrite();
  return Cmd::OK;
}
// -----------------------------------------------------------------------------

/// Read data from file into master DataSetList
Cmd::RetType ReadData(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  DataFile dataIn;
  dataIn.SetDebug( State.DFL()->Debug() );
  std::string filenameIn = argIn.GetStringNext();
  File::NameArray fnames = File::ExpandToFilenames( filenameIn );
  if (fnames.empty()) {
    mprinterr("Error: '%s' matches no files.\n", filenameIn.c_str());
    return Cmd::ERR;
  }
  int err = 0;
  for (File::NameArray::const_iterator fn = fnames.begin(); fn != fnames.end(); ++fn) {
    if (dataIn.ReadDataIn( *fn, argIn, *State.DSL() )!=0) {
      mprinterr("Error: Could not read data file '%s'.\n", fn->full());
      err++;
    }
  }
  if (err > 0) return Cmd::ERR;
  return Cmd::OK;
}

/// Exit
Cmd::RetType Quit(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  return Cmd::QUIT;
}

/// Run a system command
Cmd::RetType SystemCmd(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  int err = system( argIn.ArgLine() );
  if (err != 0) mprintf("Warning: '%s' returned %i\n", argIn.Command(), err);
  return Cmd::OK;
}

// -----------------------------------------------------------------------------
static void Help_Help() {
  mprintf("\t[{ <cmd> | <category>}]\n\tCategories:");
  for (int i = 1; i != (int)Cmd::DEPRECATED; i++)
    mprintf(" %s", Command::CommandCategoryKeyword(i));
  mprintf("\n");
  mprintf("  With no arguments list all known commands, otherwise display help for specified\n"
          "  command. If a category is specified list only commands in that category.\n");
}

/// Find help for command/topic
Cmd::RetType Help(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  ArgList arg = argIn;
  arg.RemoveFirstArg();
  if (arg.empty())
    // Cmd::NONE in this context means list all commands
    Command::ListCommands(Cmd::NONE);
  else {
    for (int i = 1; i != (int)Cmd::DEPRECATED; i++) {
      if (arg.CommandIs(Command::CommandCategoryKeyword(i))) {
        Command::ListCommands((Cmd::Ctype)i);
        return Cmd::OK;
      }
    }
    // Find help for specified command.
    Cmd::TokenPtr dispatchToken = Command::SearchToken( arg );
    if (dispatchToken == 0 || dispatchToken->Help == 0)
      mprinterr("No help found for %s\n", arg.Command());
    else
      dispatchToken->Help();
  }
  return Cmd::OK;
}

/// Run the current State
Cmd::RetType RunState(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  return (Cmd::RetType)State.Run();
}

/// Read input from a file.
Cmd::RetType ReadInput(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  // Next arg should be a filename. Not allowed to be blank in command.
  std::string inputFilename = argIn.GetStringNext();
  if (inputFilename.empty()) {
    mprinterr("Error: No input filename given.\n");
    return Cmd::ERR;
  }
  return Command::ProcessInput(State, inputFilename);
}

/// Tell CpptrajState to ignore errors if possible
Cmd::RetType NoExitOnError(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  State.SetNoExitOnError();
  mprintf("\tAttempting to ignore errors if possible.\n");
  return Cmd::OK;
}

/// Tell CpptrajState not to use a progress bar during Run.
Cmd::RetType NoProgress(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  State.SetNoProgress();
  mprintf("\tProgress bar will not be used during Run.\n");
  return Cmd::OK;
}

///  Set precision for specific set or all sets in specified DataFile
Cmd::RetType Precision(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  // Next string is DataSet(s)/DataFile that command pertains to.
  std::string name1 = argIn.GetStringNext();
  if (name1.empty()) {
    mprinterr("Error: No filename/setname given.\n");
    return Cmd::ERR;
  }
  // This will break if dataset name starts with a digit...
  int width = argIn.getNextInteger(12);
  if (width < 1) {
    mprintf("Error: Cannot set width < 1 (%i).\n", width);
    return Cmd::ERR;
  }
  int precision = argIn.getNextInteger(4);
  if (precision < 0) precision = 0;
  DataFile* df = State.DFL()->GetDataFile(name1);
  if (df != 0) {
    mprintf("\tSetting precision for all sets in %s to %i.%i\n", df->DataFilename().base(),
            width, precision);
    df->SetDataFilePrecision(width, precision);
  } else {
    State.DSL()->SetPrecisionOfDataSets( name1, width, precision );
  }
  return Cmd::OK;
}

/// Run specified analysis or all analyses in State.
Cmd::RetType RunAnalysis(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  // If only 1 arg (the command) run all analyses in list
  if (argIn.Nargs() == 1) {
    int eval = State.RunAnalyses();
    State.MasterDataFileWrite();
    if (eval == 0)
      return Cmd::OK;
    else
      return Cmd::ERR;
  }
  // Run specified analysis
  // FIXME: Use RemoveFirstArg
  ArgList analyzeargs = argIn.RemainingArgs();
  analyzeargs.MarkArg(0);
  Cmd::TokenPtr tkn = Command::SearchTokenType( Cmd::ANALYSIS, analyzeargs );
  if ( tkn == 0 ) return Cmd::ERR;
  Analysis* ana = (Analysis*)tkn->Alloc();
  if (ana == 0) return Cmd::ERR;
  Timer total_time;
  total_time.Start();
  Cmd::RetType err = Cmd::ERR;
  if ( ana->Setup( analyzeargs, State.DSL(), State.DFL(), State.Debug() ) == Analysis::OK )
  {
    analyzeargs.CheckForMoreArgs();
    if (ana->Analyze() != Analysis::ERR) {
      err = Cmd::OK;
      State.MasterDataFileWrite();
    }
  }
  delete ana;
  total_time.Stop();
  mprintf("TIME: Total analysis execution time: %.4f seconds.\n", total_time.Total());
  return err;
}

/// Show results of mask expression
Cmd::RetType SelectAtoms(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  AtomMask tempMask( argIn.GetMaskNext() );
  Topology* parm = State.DSL()->GetTopByIndex( argIn );
  if (parm == 0) return Cmd::ERR;
  if (parm->SetupIntegerMask( tempMask )) return Cmd::ERR;
  mprintf("Selected %i atoms.\n", tempMask.Nselected());
  if (!argIn.hasKey("total"))
    tempMask.PrintMaskAtoms("Selected");
  return Cmd::OK;
}

/// Show results of DataSet expression
Cmd::RetType SelectDataSets(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  std::string dsarg = argIn.GetStringNext();
  DataSetList dsets = State.DSL()->GetMultipleSets( dsarg );
  if (!dsets.empty()) {
    mprintf("SelectDS: Arg '%s':", dsarg.c_str());
    dsets.List();
  }
  return Cmd::OK;
}

// -----------------------------------------------------------------------------
static void Help_PrintData() {
  mprintf("\t<data set>\n"
          "  Print data from data set to screen.\n");
}

Cmd::RetType PrintData(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  DataFile ToStdout;
  ToStdout.SetupStdout(argIn, State.Debug());
  DataSetList selected; 
  std::string ds_arg = argIn.GetStringNext();
  while (!ds_arg.empty()) {
    selected += State.DSL()->GetMultipleSets( ds_arg );
    ds_arg = argIn.GetStringNext();
  }
  for (DataSetList::const_iterator ds = selected.begin(); ds != selected.end(); ++ds)
    ToStdout.AddDataSet( *ds );
  ToStdout.WriteDataOut();
  return Cmd::OK;
}

// -----------------------------------------------------------------------------
static void Help_Calc() {
  mprintf("\t<expression>\n"
          "  Evaluate the given mathematical expression.\n");
}

/// Parse a mathematical expression.
Cmd::RetType Calc(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  RPNcalc calc;
  calc.SetDebug( State.Debug() );
  // Do NOT include command in expression.
  if (calc.ProcessExpression( argIn.ArgString().substr(argIn[0].size()) ))
    return Cmd::ERR;
  if (calc.Evaluate(*State.DSL())) return Cmd::ERR;
  return Cmd::OK;
}
  
// ---------- TRAJECTORY COMMANDS ----------------------------------------------
/// Add output trajectory to State
Cmd::RetType Trajout(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  return (Cmd::RetType)State.AddOutputTrajectory( argIn );
}

/// Add input trajectory to State
Cmd::RetType Trajin(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  return (Cmd::RetType)State.AddTrajin( argIn, false );
}

/// Add ensemble of input trajectories to State.
Cmd::RetType Ensemble(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  return (Cmd::RetType)State.AddTrajin( argIn, true );
}

/// Add reference trajectory to State
Cmd::RetType Reference(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  return (Cmd::RetType)State.AddReference( argIn.GetStringNext(), argIn );
}

// ---------- TOPOLOGY COMMANDS ------------------------------------------------
static void Help_LoadParm() {
  mprintf("\t<filename> [<tag>] [nobondsearch | bondsearch [<offset>]]\n"
          "  Add <filename> to the topology list.\n");
  ParmFile::ReadOptions();
}
/// Load topology from file into to State
Cmd::RetType LoadParm(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  return (Cmd::RetType)State.AddTopology(argIn.GetStringNext(), argIn);
}

static void Help_ParmInfo() {
  mprintf("\t[%s] [<mask>]\n", DataSetList::TopIdxArgs);
  mprintf("  Print information on specfied topology (first by default).\n");
}
/// Print info for specified parm or atoms in specified parm.
Cmd::RetType ParmInfo(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  Topology* parm = State.DSL()->GetTopByIndex( argIn );
  if (parm == 0) return Cmd::ERR;
  parm->Summary();
  return Cmd::OK;
}

static void Help_AtomInfo() {
  mprintf("\t[%s] [<mask>]\n", DataSetList::TopIdxArgs);
  mprintf("  Print information on atoms in <mask> for specified topology (first by default).\n");
}
/// Print info for atoms in mask.
Cmd::RetType AtomInfo(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  Topology* parm = State.DSL()->GetTopByIndex( argIn );
  if (parm == 0) return Cmd::ERR;
  parm->PrintAtomInfo( argIn.GetMaskNext() );
  return Cmd::OK;
}

static void Help_BondInfo() {
  mprintf("\t[%s] [<mask>]\n", DataSetList::TopIdxArgs);
  mprintf("  Print bond info for atoms in <mask> for specified topology (first by default).\n");
}
/// Print bond info for atoms in mask.
Cmd::RetType BondInfo(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  Topology* parm = State.DSL()->GetTopByIndex( argIn );
  if (parm == 0) return Cmd::ERR;
  parm->PrintBondInfo( argIn.GetMaskNext() );
  return Cmd::OK;
}

static void Help_AngleInfo() {
  mprintf("\t[%s] [<mask>]\n", DataSetList::TopIdxArgs);
  mprintf("  Print angle info for atoms in <mask> for specified topology (first by default).\n");
}
/// Print angle info for atoms in mask.
Cmd::RetType AngleInfo(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  Topology* parm = State.DSL()->GetTopByIndex( argIn );
  if (parm == 0) return Cmd::ERR;
  parm->PrintAngleInfo( argIn.GetMaskNext() );
  return Cmd::OK;
}

static void Help_DihedralInfo() {
  mprintf("\t[%s] [<mask>] [and]\n", DataSetList::TopIdxArgs);
  mprintf("  Print dihedral info for atoms in <mask> for specified topology (first by default).\n"
          "  If 'and' is specified dihedral must include all atoms specfied by <mask>.\n");
}
/// Print dihedral info for atoms in mask.
Cmd::RetType DihedralInfo(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  Topology* parm = State.DSL()->GetTopByIndex( argIn );
  if (parm == 0) return Cmd::ERR;
  parm->PrintDihedralInfo( argIn.GetMaskNext(), !argIn.hasKey("and") );
  return Cmd::OK;
}

static void Help_ResInfo() {
  mprintf("\t[%s] [<mask>] [short]\n", DataSetList::TopIdxArgs);
  mprintf("  Print info for residues in <mask> for specified topology (first by default).\n"
          "  If 'short' is specified print residue info in shorter form.\n");
}
/// Print residue info for atoms in mask.
Cmd::RetType ResInfo(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  Topology* parm = State.DSL()->GetTopByIndex( argIn );
  if (parm == 0) return Cmd::ERR;
  bool printShort = argIn.hasKey("short");
  if (printShort)
    parm->PrintShortResInfo( argIn.GetMaskNext(), argIn.getKeyInt("maxwidth",50) );
  else
    parm->PrintResidueInfo( argIn.GetMaskNext() );
  return Cmd::OK;
}

static void Help_MolInfo() {
  mprintf("\t[%s] [<mask>]\n", DataSetList::TopIdxArgs);
  mprintf("  Print info for molecules in <mask> for specfied topology (first by default).\n");
}
/// Print molecule info for atoms in mask.
Cmd::RetType MolInfo(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  Topology* parm = State.DSL()->GetTopByIndex( argIn );
  if (parm == 0) return Cmd::ERR;
  parm->PrintMoleculeInfo( argIn.GetMaskNext() );
  return Cmd::OK;
}

static void Help_ChargeInfo() {
  mprintf("\t[%s] <mask>\n", DataSetList::TopIdxArgs);
  mprintf("  Print total charge of atoms in <mask> for specified topology (first by default).\n");
}
/// Print the total charge of atoms in mask
Cmd::RetType ChargeInfo(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  Topology* parm = State.DSL()->GetTopByIndex( argIn );
  if (parm == 0) return Cmd::ERR;
  if (parm->PrintChargeMassInfo( argIn.GetMaskNext(), 0 )) return Cmd::ERR;
  return Cmd::OK;
}

static void Help_MassInfo() {
  mprintf("\t[%s] <mask>\n", DataSetList::TopIdxArgs);
  mprintf("  Print total mass of atoms in <mask> for specified topology (first by default).\n");
}
/// Print the total mass of atoms in mask
Cmd::RetType MassInfo(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  Topology* parm = State.DSL()->GetTopByIndex( argIn );
  if (parm == 0) return Cmd::ERR;
  if (parm->PrintChargeMassInfo( argIn.GetMaskNext(), 1 )) return Cmd::ERR;
  return Cmd::OK;
}

static void Help_ParmBox() {
  mprintf("\t[%s] [nobox] [truncoct]\n", DataSetList::TopIdxArgs);
  mprintf("\t[x <xval>] [y <yval>] [z <zval>] [alpha <a>] [beta <b>] [gamma <g>]\n"
          "  Set the box info for specified topology (currently only relevant for Amber\n"
          "  Topology/ASCII coords). If 'nobox' is specified, remove box info. If\n"
          "  'truncoct' specified, set truncated octahedron with lengths = <xval>.\n");
}
/// Modify specified parm box info
Cmd::RetType ParmBox(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
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
  Topology* parm = State.DSL()->GetTopByIndex( argIn );
  if (parm == 0) return Cmd::ERR;
  if (nobox)
    mprintf("\tRemoving box information from parm %i:%s\n", parm->Pindex(), parm->c_str());
  else
    // Fill in missing parm box information from specified parm
    pbox.SetMissingInfo( parm->ParmBox() );
  if (argIn.hasKey("truncoct")) pbox.SetTruncOct();
  parm->SetParmBox( pbox );
  parm->ParmBox().PrintInfo();
  return Cmd::OK;
}

static void Help_ParmStrip() {
  mprintf("\t<mask> [%s]\n", DataSetList::TopIdxArgs);
  mprintf("  Strip atoms in mask from specified topology (first by default).\n");
}
/// Strip atoms from specified parm
Cmd::RetType ParmStrip(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  Topology* parm = State.DSL()->GetTopByIndex( argIn );
  if (parm == 0) return Cmd::ERR;
  // Check if this topology has already been used to set up an input
  // trajectory, as this will break the traj read.
  bool topology_in_use = false;
  const char* fname = 0;
  for (TrajinList::trajin_it tIn = State.InputTrajList().trajin_begin();
                             tIn != State.InputTrajList().trajin_end(); ++tIn)
    if ( (*tIn)->Traj().Parm() == parm ) {
      topology_in_use = true;
      fname = (*tIn)->Traj().Filename().full();
      break;
    }
  if (!topology_in_use) {
    for (TrajinList::ensemble_it eIn = State.InputTrajList().ensemble_begin();
                                 eIn != State.InputTrajList().ensemble_end(); ++eIn)
      if ( (*eIn)->Traj().Parm() == parm ) {
        topology_in_use = true;
        fname = (*eIn)->Traj().Filename().full();
        break;
      }
  }
  if (topology_in_use) {
    mprinterr("Error: Topology '%s' has already been used to set up trajectory '%s'.\n"
              "Error:   To strip this topology use the 'strip' action.\n",
              parm->c_str(), fname);
    return Cmd::ERR;
  }
  AtomMask tempMask( argIn.GetMaskNext() );
  // Since want to keep atoms outside mask, invert selection
  tempMask.InvertMaskExpression();
  if (parm->SetupIntegerMask( tempMask )) return Cmd::ERR;
  mprintf("\tStripping atoms in mask [%s] (%i) from %s\n",tempMask.MaskString(),
           parm->Natom() - tempMask.Nselected(), parm->c_str());
  Topology* tempParm = parm->modifyStateByMask(tempMask);
  if (tempParm==0) {
    mprinterr("Error: %s: Could not strip parm.\n", argIn.Command());
    return Cmd::ERR;
  } else {
    // Replace parm with stripped version
    *parm = *tempParm;
    parm->Brief("Stripped parm:");
    delete tempParm;
  }
  return Cmd::OK;
}

static void Help_ParmWrite() {
  mprintf("\tout <filename> [{%s | crdset <setname>}] [<fmt>] [nochamber]\n",
          DataSetList::TopIdxArgs);
  mprintf("  Write specified topology or topology from COORDS set to <filename>.\n");
  ParmFile::WriteOptions();
}
/// Write parm to Amber Topology file.
Cmd::RetType ParmWrite(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  std::string outfilename = argIn.GetStringKey("out");
  if (outfilename.empty()) {
    mprinterr("Error: No output filename specified (use 'out <filename>').\n");
    return Cmd::ERR;
  }
  int err = 0;
  ParmFile pfile;
  // Check if a COORDS data set was specified.
  std::string crdset = argIn.GetStringKey("crdset");
  if (crdset.empty()) {
    Topology* parm = State.DSL()->GetTopByIndex( argIn );
    if (parm == 0) return Cmd::ERR;
    err = pfile.WriteTopology( *parm, outfilename, argIn, ParmFile::UNKNOWN_PARM, State.Debug() );
  } else {
    DataSet_Coords* ds = (DataSet_Coords*)State.DSL()->FindCoordsSet(crdset);
    if (ds == 0) return Cmd::ERR;
    mprintf("\tUsing topology from data set '%s'\n", ds->legend());
    err = pfile.WriteTopology(ds->Top(), outfilename, argIn, ParmFile::UNKNOWN_PARM, State.Debug());
  }
  if (err != 0)
    return Cmd::ERR;
  return Cmd::OK;
}

static void Help_ParmSolvent() {
  mprintf("\t[%s] { <mask> | none }\n", DataSetList::TopIdxArgs);
  mprintf("  Set solvent for the specified topology (default first) based on <mask>.\n"
          "  If 'none' specified, remove all solvent information.\n");
}
/// Modify parm solvent information
Cmd::RetType ParmSolvent(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  std::string maskexpr;
  if (!argIn.hasKey("none")) {
    maskexpr = argIn.GetMaskNext();
    if ( maskexpr.empty() ) {
      mprinterr("Error: solvent: No mask specified.\n");
      return Cmd::ERR;
    }
  }
  // Get parm index
  Topology* parm = State.DSL()->GetTopByIndex( argIn );
  if (parm == 0) return Cmd::ERR;
  parm->SetSolvent( maskexpr );
  return Cmd::OK;
}

static void Help_ScaleDihedralK() {
  mprintf("\t[%s] <scale factor> [<mask> [useall]]\n", DataSetList::TopArgs);
}
/// Scale dihedral force constants in specfied parm by factor.
Cmd::RetType ScaleDihedralK(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  Topology* parm = State.DSL()->GetTopology( argIn );
  if (parm == 0) {
    mprinterr("Error: No topologies loaded.\n");
    return Cmd::ERR;
  }
  double scale_factor = argIn.getNextDouble(1.0);
  std::string maskexpr = argIn.GetMaskNext();
  bool useAll = argIn.hasKey("useall");
  mprintf("\tScaling dihedral force constants in %s by %f\n", parm->c_str(), scale_factor);
  if (!maskexpr.empty()) {
    if (useAll)
      mprintf("\tAll atoms in mask '%s' must be present to select dihedral.\n",maskexpr.c_str());
    else
      mprintf("\tAny atom in mask '%s' will select a dihedral.\n",maskexpr.c_str());
  }
  parm->ScaleDihedralK( scale_factor, maskexpr, useAll );
  return Cmd::OK;
}

// ---------- DISPATCHABLE COMMANDS --------------------------------------------
/// Add an action to the State ActionList
Cmd::RetType AddAction(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  return ( (Cmd::RetType)State.AddAction( Alloc, argIn ) );
}

/// Add an action to the State AnalysisList
Cmd::RetType AddAnalysis(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  if (State.AddAnalysis( Alloc, argIn )) {
#   ifndef MPI
    if (State.InputTrajList().Mode() == TrajinList::ENSEMBLE)
      mprinterr("Info: Data sets for ensemble members beyond the first (member 0) have not\n"
                "Info:   yet been created for current Actions. If any Analyses report warnings\n"
                "Info:   or errors related to missing data sets, try entering a 'run' command\n"
                "Info:   prior to any analysis commands.\n");
#   endif
    return Cmd::ERR;
  }
  return Cmd::OK;
}

// ================ LIST OF ALL COMMANDS =======================================
/** Ideally keep this array first sorted by type (1st field), then 
  * alphabetically by command string (2nd field).
  */
const Cmd::Token Command::Commands[] = {
  // GENERAL COMMANDS
  { Cmd::GENERAL, "activeref",     0, Help_ActiveRef,       ActiveRef       },
  { Cmd::GENERAL, "calc",          0, Help_Calc,            Calc            },
  { Cmd::GENERAL, "clear",         0, Help_Clear,           ClearList       },
  { Cmd::GENERAL, "create",        0, Help_Create_DataFile, Create_DataFile },
  { Cmd::GENERAL, "datafile",      0, Help_DataFile,        DataFileCmd     },
  { Cmd::GENERAL, "datafilter",    0, Help_DataFilter,      DataFilter      },
  { Cmd::GENERAL, "dataset",       0, Help_DataSetCmd,      DataSetCmd      },
  { Cmd::GENERAL, "debug",         0, Help_Debug,           SetListDebug    },
  { Cmd::GENERAL, "exit" ,         0, Help_Quit,            Quit            },
  { Cmd::GENERAL, "go",            0, Help_Run,             RunState        },
  { Cmd::GENERAL, "help",          0, Help_Help,            Help            },
  { Cmd::GENERAL, "list",          0, Help_List,            ListAll         },
  { Cmd::GENERAL, "noexitonerror", 0, Help_NoExitOnError,   NoExitOnError   },
  { Cmd::GENERAL, "noprogress",    0, Help_NoProgress,      NoProgress      },
  { Cmd::GENERAL, "precision",     0, Help_Precision,       Precision       },
  { Cmd::GENERAL, "printdata",     0, Help_PrintData,       PrintData       },
  { Cmd::GENERAL, "prnlev",        0, Help_Debug,           SetListDebug    },
  { Cmd::GENERAL, "quit" ,         0, Help_Quit,            Quit            },
  { Cmd::GENERAL, "readdata",      0, Help_ReadData,        ReadData        },
  { Cmd::GENERAL, "readinput",     0, Help_ReadInput,       ReadInput       },
  { Cmd::GENERAL, "removedata",    0, Help_RemoveData,      RemoveData      },
  { Cmd::GENERAL, "rst"   ,        0, Help_GenerateAmberRst,GenerateAmberRst},
  { Cmd::GENERAL, "run"   ,        0, Help_Run,             RunState        },
  { Cmd::GENERAL, "runanalysis",   0, Help_RunAnalysis,     RunAnalysis     },
  { Cmd::GENERAL, "select",        0, Help_Select,          SelectAtoms     },
  { Cmd::GENERAL, "selectds",      0, Help_SelectDS,        SelectDataSets  },
  { Cmd::HIDDEN,  "sequencealign", 0, Help_SequenceAlign,   SequenceAlign   },
  { Cmd::GENERAL, "silenceactions",0, Help_SilenceActions,  SilenceActions  },
  { Cmd::GENERAL, "write",         0, Help_Write_DataFile,  Write_DataFile  },
  { Cmd::GENERAL, "writedata",     0, Help_Write_DataFile,  Write_DataFile  },
  // SYSTEM COMMANDS
  { Cmd::SYSTEM, "gnuplot",        0, Help_System,          SystemCmd       },
  { Cmd::SYSTEM, "head",           0, Help_System,          SystemCmd       },
  { Cmd::SYSTEM, "less",           0, Help_System,          SystemCmd       },
  { Cmd::SYSTEM, "ls",             0, Help_System,          SystemCmd       },
  { Cmd::SYSTEM, "pwd",            0, Help_System,          SystemCmd       },
  { Cmd::SYSTEM, "xmgrace",        0, Help_System,          SystemCmd       },
  // COORDS COMMANDS
  { Cmd::COORDS,  "combinecrd",    0, Help_CombineCoords,   CombineCoords   },
  { Cmd::COORDS,  "crdaction",     0, Help_CrdAction,       CrdActionCmd    },
  { Cmd::COORDS,  "crdout",        0, Help_CrdOut,          CrdOut          },
  { Cmd::COORDS,  "loadcrd",       0, Help_LoadCrd,         LoadCrd         },
  { Cmd::COORDS,  "loadtraj",      0, Help_LoadTraj,        LoadTraj        },
  // TRAJECTORY COMMANDS
  { Cmd::TRAJ,    "ensemble",      0, Help_Ensemble,        Ensemble        },
  { Cmd::TRAJ,    "reference",     0, Help_Reference,       Reference       },
  { Cmd::TRAJ,    "trajin",        0, Help_Trajin,          Trajin          },
  { Cmd::TRAJ,    "trajout",       0, Help_Trajout,         Trajout         },
  // TOPOLOGY COMMANDS
  { Cmd::PARM,    "angleinfo",     0, Help_AngleInfo,       AngleInfo       },
  { Cmd::PARM,    "angles",        0, Help_AngleInfo,       AngleInfo       },
  { Cmd::PARM,    "atominfo",      0, Help_AtomInfo,        AtomInfo        },
  { Cmd::PARM,    "atoms",         0, Help_AtomInfo,        AtomInfo        },
  { Cmd::PARM,    "bondinfo",      0, Help_BondInfo,        BondInfo        },
  { Cmd::PARM,    "bonds",         0, Help_BondInfo,        BondInfo        },
  { Cmd::PARM,    "charge",        0, Help_ChargeInfo,      ChargeInfo      },
  { Cmd::HIDDEN,  "comparetop",    0, Help_CompareTop,      CompareTop      },
  { Cmd::PARM,    "dihedralinfo",  0, Help_DihedralInfo,    DihedralInfo    },
  { Cmd::PARM,    "dihedrals",     0, Help_DihedralInfo,    DihedralInfo    },
  { Cmd::PARM,    "mass",          0, Help_MassInfo,        MassInfo        },
  { Cmd::PARM,    "molinfo",       0, Help_MolInfo,         MolInfo         },
  { Cmd::PARM,    "parm",          0, Help_LoadParm,        LoadParm        },
  { Cmd::PARM,    "parmbox",       0, Help_ParmBox,         ParmBox         },
  { Cmd::PARM,    "parminfo",      0, Help_ParmInfo,        ParmInfo        },
  { Cmd::PARM,    "parmstrip",     0, Help_ParmStrip,       ParmStrip       },
  { Cmd::PARM,    "parmwrite",     0, Help_ParmWrite,       ParmWrite       },
  { Cmd::PARM,    "printangles",   0, Help_AngleInfo,       AngleInfo       },
  { Cmd::PARM,    "printatoms",    0, Help_AtomInfo,        AtomInfo        },
  { Cmd::PARM,    "printbonds",    0, Help_BondInfo,        BondInfo        },
  { Cmd::PARM,    "printdihedrals",0, Help_DihedralInfo,    DihedralInfo    },
  { Cmd::PARM,    "resinfo",       0, Help_ResInfo,         ResInfo         },
  { Cmd::PARM,    "scaledihedralk",0, Help_ScaleDihedralK,  ScaleDihedralK  },
  { Cmd::PARM,    "solvent",       0, Help_ParmSolvent,     ParmSolvent     },
  // INC_ACTION: ACTION COMMANDS
  { Cmd::ACTION, "angle", Action_Angle::Alloc, Action_Angle::Help, AddAction },
  { Cmd::ACTION, "areapermol", Action_AreaPerMol::Alloc, Action_AreaPerMol::Help, AddAction },
  { Cmd::ACTION, "atomiccorr", Action_AtomicCorr::Alloc, Action_AtomicCorr::Help, AddAction },
  { Cmd::ACTION, "atomicfluct", Action_AtomicFluct::Alloc, Action_AtomicFluct::Help, AddAction },
  { Cmd::ACTION, "atommap", Action_AtomMap::Alloc, Action_AtomMap::Help, AddAction },
  { Cmd::ACTION, "autoimage", Action_AutoImage::Alloc, Action_AutoImage::Help, AddAction },
  { Cmd::ACTION, "average", Action_Average::Alloc, Action_Average::Help, AddAction },
  { Cmd::ACTION, "bounds", Action_Bounds::Alloc, Action_Bounds::Help, AddAction },
  { Cmd::ACTION, "box", Action_Box::Alloc, Action_Box::Help, AddAction },
  { Cmd::ACTION, "center", Action_Center::Alloc, Action_Center::Help, AddAction },
  { Cmd::ACTION, "channel", Action_Channel::Alloc, Action_Channel::Help, AddAction },
  { Cmd::ACTION, "check", Action_CheckStructure::Alloc, Action_CheckStructure::Help, AddAction },
  { Cmd::ACTION, "checkchirality", Action_CheckChirality::Alloc, Action_CheckChirality::Help, AddAction },
  { Cmd::ACTION, "checkoverlap", Action_CheckStructure::Alloc, Action_CheckStructure::Help, AddAction },
  { Cmd::ACTION, "checkstructure", Action_CheckStructure::Alloc, Action_CheckStructure::Help, AddAction },
  { Cmd::ACTION, "closest", Action_Closest::Alloc, Action_Closest::Help, AddAction },
  { Cmd::ACTION, "closestwaters", Action_Closest::Alloc, Action_Closest::Help, AddAction },
  { Cmd::ACTION, "clusterdihedral", Action_ClusterDihedral::Alloc, Action_ClusterDihedral::Help, AddAction },
  { Cmd::ACTION, "contacts", Action_Contacts::Alloc, Action_Contacts::Help, AddAction },
  { Cmd::ACTION, "createcrd", Action_CreateCrd::Alloc, Action_CreateCrd::Help, AddAction },
  { Cmd::ACTION, "createreservoir", Action_CreateReservoir::Alloc, Action_CreateReservoir::Help, AddAction },
  { Cmd::ACTION, "density", Action_Density::Alloc, Action_Density::Help, AddAction },
  { Cmd::ACTION, "diffusion", Action_Diffusion::Alloc, Action_Diffusion::Help, AddAction },
  { Cmd::ACTION, "dihedral", Action_Dihedral::Alloc, Action_Dihedral::Help, AddAction },
  { Cmd::ACTION, "dihedralscan", Action_DihedralScan::Alloc, Action_DihedralScan::Help, AddAction },
  { Cmd::ACTION, "dipole", Action_Dipole::Alloc, Action_Dipole::Help, AddAction },
  { Cmd::ACTION, "distance", Action_Distance::Alloc, Action_Distance::Help, AddAction },
//  { Cmd::ACTION, "dnaiontracker", Action_DNAionTracker::Alloc, Action_DNAionTracker::Help, AddAction },
  { Cmd::ACTION, "drms", Action_DistRmsd::Alloc, Action_DistRmsd::Help, AddAction },
  { Cmd::ACTION, "drmsd", Action_DistRmsd::Alloc, Action_DistRmsd::Help, AddAction },
  { Cmd::ACTION, "dssp", Action_DSSP::Alloc, Action_DSSP::Help, AddAction },
  { Cmd::ACTION, "energy", Action_Energy::Alloc, Action_Energy::Help, AddAction },
  { Cmd::ACTION, "filter", Action_FilterByData::Alloc, Action_FilterByData::Help, AddAction },
  { Cmd::ACTION, "fixatomorder", Action_FixAtomOrder::Alloc, Action_FixAtomOrder::Help, AddAction },
  { Cmd::ACTION, "gist", Action_Gist::Alloc, Action_Gist::Help, AddAction },
//  { Cmd::ACTION, "gfe", Action_GridFreeEnergy::Alloc, Action_GridFreeEnergy::Help, AddAction },
  { Cmd::ACTION, "grid", Action_Grid::Alloc, Action_Grid::Help, AddAction },
  { Cmd::ACTION, "hbond", Action_Hbond::Alloc, Action_Hbond::Help, AddAction },
  { Cmd::ACTION, "image", Action_Image::Alloc, Action_Image::Help, AddAction },
  { Cmd::ACTION, "jcoupling", Action_Jcoupling::Alloc, Action_Jcoupling::Help, AddAction },
  { Cmd::ACTION, "lessplit", Action_LESsplit::Alloc, Action_LESsplit::Help, AddAction },
  { Cmd::ACTION, "lie", Action_LIE::Alloc, Action_LIE::Help, AddAction },
  { Cmd::ACTION, "lipidorder", Action_OrderParameter::Alloc, Action_OrderParameter::Help, AddAction },
  { Cmd::ACTION, "makestructure", Action_MakeStructure::Alloc, Action_MakeStructure::Help, AddAction },
  { Cmd::ACTION, "mask", Action_Mask::Alloc, Action_Mask::Help, AddAction },
  { Cmd::ACTION, "matrix", Action_Matrix::Alloc, Action_Matrix::Help, AddAction },
  { Cmd::ACTION, "minimage", Action_MinImage::Alloc, Action_MinImage::Help, AddAction },
  { Cmd::ACTION, "molsurf", Action_Molsurf::Alloc, Action_Molsurf::Help, AddAction },
  { Cmd::ACTION, "multidihedral", Action_MultiDihedral::Alloc, Action_MultiDihedral::Help, AddAction },
  { Cmd::ACTION, "multivector", Action_MultiVector::Alloc, Action_MultiVector::Help, AddAction },
  { Cmd::ACTION, "nastruct", Action_NAstruct::Alloc, Action_NAstruct::Help, AddAction },
  { Cmd::ACTION, "nativecontacts", Action_NativeContacts::Alloc, Action_NativeContacts::Help, AddAction },
  { Cmd::ACTION, "nmrrst", Action_NMRrst::Alloc, Action_NMRrst::Help, AddAction },
  { Cmd::ACTION, "outtraj", Action_Outtraj::Alloc, Action_Outtraj::Help, AddAction },
  { Cmd::ACTION, "pairdist", Action_PairDist::Alloc, Action_PairDist::Help, AddAction },
  { Cmd::ACTION, "pairwise", Action_Pairwise::Alloc, Action_Pairwise::Help, AddAction },
  { Cmd::ACTION, "principal", Action_Principal::Alloc, Action_Principal::Help, AddAction },
  { Cmd::ACTION, "projection", Action_Projection::Alloc, Action_Projection::Help, AddAction },
  { Cmd::ACTION, "pucker", Action_Pucker::Alloc, Action_Pucker::Help, AddAction },
  { Cmd::ACTION, "radgyr", Action_Radgyr::Alloc, Action_Radgyr::Help, AddAction },
  { Cmd::ACTION, "radial", Action_Radial::Alloc, Action_Radial::Help, AddAction },
  { Cmd::ACTION, "randomizeions", Action_RandomizeIons::Alloc, Action_RandomizeIons::Help, AddAction },
  { Cmd::ACTION, "rdf", Action_Radial::Alloc, Action_Radial::Help, AddAction },
  { Cmd::ACTION, "replicatecell", Action_ReplicateCell::Alloc, Action_ReplicateCell::Help, AddAction },
  { Cmd::ACTION, "rms", Action_Rmsd::Alloc, Action_Rmsd::Help, AddAction },
  { Cmd::ACTION, "rmsd", Action_Rmsd::Alloc, Action_Rmsd::Help, AddAction },
  { Cmd::ACTION, "rmsf", Action_AtomicFluct::Alloc, Action_AtomicFluct::Help, AddAction },
  { Cmd::ACTION, "rog", Action_Radgyr::Alloc, Action_Radgyr::Help, AddAction },
  { Cmd::ACTION, "rotate", Action_Rotate::Alloc, Action_Rotate::Help, AddAction },
  { Cmd::ACTION, "runavg", Action_RunningAvg::Alloc, Action_RunningAvg::Help, AddAction },
  { Cmd::ACTION, "runningaverage", Action_RunningAvg::Alloc, Action_RunningAvg::Help, AddAction },
  { Cmd::ACTION, "scale", Action_Scale::Alloc, Action_Scale::Help, AddAction },
  { Cmd::ACTION, "secstruct", Action_DSSP::Alloc, Action_DSSP::Help, AddAction },
  { Cmd::ACTION, "setvelocity", Action_SetVelocity::Alloc, Action_SetVelocity::Help, AddAction },
  { Cmd::ACTION, "spam", Action_Spam::Alloc, Action_Spam::Help, AddAction },
  { Cmd::ACTION, "stfcdiffusion", Action_STFC_Diffusion::Alloc, Action_STFC_Diffusion::Help, AddAction },
  { Cmd::ACTION, "strip", Action_Strip::Alloc, Action_Strip::Help, AddAction },
  { Cmd::ACTION, "surf", Action_Surf::Alloc, Action_Surf::Help, AddAction },
  { Cmd::ACTION, "symmrmsd", Action_SymmetricRmsd::Alloc, Action_SymmetricRmsd::Help, AddAction },
  { Cmd::ACTION, "temperature", Action_Temperature::Alloc, Action_Temperature::Help, AddAction },
  { Cmd::ACTION, "trans", Action_Translate::Alloc, Action_Translate::Help, AddAction },
  { Cmd::ACTION, "translate", Action_Translate::Alloc, Action_Translate::Help, AddAction },
  { Cmd::ACTION, "unstrip", Action_Unstrip::Alloc, Action_Unstrip::Help, AddAction },
  { Cmd::ACTION, "unwrap", Action_Unwrap::Alloc, Action_Unwrap::Help, AddAction },
  { Cmd::ACTION, "vector", Action_Vector::Alloc, Action_Vector::Help, AddAction },
  { Cmd::ACTION, "velocityautocorr", Action_VelocityAutoCorr::Alloc, Action_VelocityAutoCorr::Help, AddAction },
  { Cmd::ACTION, "volmap", Action_Volmap::Alloc, Action_Volmap::Help, AddAction},
  { Cmd::ACTION, "volume", Action_Volume::Alloc, Action_Volume::Help, AddAction},
  { Cmd::ACTION, "watershell", Action_Watershell::Alloc, Action_Watershell::Help, AddAction },
  // INC_ANALYSIS: ANALYSIS COMMANDS
  { Cmd::ANALYSIS, "2drms", Analysis_Rms2d::Alloc, Analysis_Rms2d::Help, AddAnalysis },
  { Cmd::ANALYSIS, "amdbias", Analysis_AmdBias::Alloc, Analysis_AmdBias::Help, AddAnalysis },
  { Cmd::ANALYSIS, "autocorr", Analysis_AutoCorr::Alloc, Analysis_AutoCorr::Help, AddAnalysis },
  { Cmd::ANALYSIS, "avg", Analysis_Average::Alloc, Analysis_Average::Help, AddAnalysis },
  { Cmd::ANALYSIS, "calcstate", Analysis_State::Alloc, Analysis_State::Help, AddAnalysis },
  { Cmd::ANALYSIS, "cluster", Analysis_Clustering::Alloc, Analysis_Clustering::Help, AddAnalysis },
  { Cmd::ANALYSIS, "corr", Analysis_Corr::Alloc, Analysis_Corr::Help, AddAnalysis },
  { Cmd::ANALYSIS, "correlationcoe", Analysis_Corr::Alloc, Analysis_Corr::Help, AddAnalysis },
  { Cmd::ANALYSIS, "crank", Analysis_CrankShaft::Alloc, Analysis_CrankShaft::Help, AddAnalysis },
  { Cmd::ANALYSIS, "crankshaft", Analysis_CrankShaft::Alloc, Analysis_CrankShaft::Help, AddAnalysis },
  { Cmd::ANALYSIS, "crdfluct", Analysis_CrdFluct::Alloc, Analysis_CrdFluct::Help, AddAnalysis },
  { Cmd::ANALYSIS, "crosscorr", Analysis_CrossCorr::Alloc, Analysis_CrossCorr::Help, AddAnalysis },
  { Cmd::ANALYSIS, "curvefit", Analysis_CurveFit::Alloc, Analysis_CurveFit::Help, AddAnalysis },
  { Cmd::ANALYSIS, "diagmatrix", Analysis_Matrix::Alloc, Analysis_Matrix::Help, AddAnalysis },
  { Cmd::ANALYSIS, "divergence", Analysis_Divergence::Alloc, Analysis_Divergence::Help, AddAnalysis },
  { Cmd::ANALYSIS, "fft", Analysis_FFT::Alloc, Analysis_FFT::Help, AddAnalysis },
  { Cmd::ANALYSIS, "hist", Analysis_Hist::Alloc, Analysis_Hist::Help, AddAnalysis },
  { Cmd::ANALYSIS, "histogram", Analysis_Hist::Alloc, Analysis_Hist::Help, AddAnalysis },
  { Cmd::ANALYSIS, "integrate", Analysis_Integrate::Alloc, Analysis_Integrate::Help, AddAnalysis },
  { Cmd::ANALYSIS, "ired", Analysis_IRED::Alloc, Analysis_IRED::Help, AddAnalysis },
  { Cmd::ANALYSIS, "kde", Analysis_KDE::Alloc, Analysis_KDE::Help, AddAnalysis },
  { Cmd::ANALYSIS, "lifetime", Analysis_Lifetime::Alloc, Analysis_Lifetime::Help, AddAnalysis },
  { Cmd::ANALYSIS, "lowestcurve", Analysis_LowestCurve::Alloc, Analysis_LowestCurve::Help, AddAnalysis },
  { Cmd::ANALYSIS, "matrix", Analysis_Matrix::Alloc, Analysis_Matrix::Help, AddAnalysis },
  { Cmd::ANALYSIS, "meltcurve", Analysis_MeltCurve::Alloc, Analysis_MeltCurve::Help, AddAnalysis },
  { Cmd::ANALYSIS, "modes", Analysis_Modes::Alloc, Analysis_Modes::Help, AddAnalysis },
  { Cmd::ANALYSIS, "multicurve", Analysis_Multicurve::Alloc, Analysis_Multicurve::Help, AddAnalysis },
  { Cmd::ANALYSIS, "multihist", Analysis_MultiHist::Alloc, Analysis_MultiHist::Help, AddAnalysis },
  { Cmd::ANALYSIS, "overlap", Analysis_Overlap::Alloc, Analysis_Overlap::Help, AddAnalysis },
  { Cmd::ANALYSIS, "phipsi", Analysis_PhiPsi::Alloc, Analysis_PhiPsi::Help, AddAnalysis },
  { Cmd::ANALYSIS, "regress", Analysis_Regression::Alloc, Analysis_Regression::Help, AddAnalysis },
  { Cmd::ANALYSIS, "remlog", Analysis_RemLog::Alloc, Analysis_RemLog::Help, AddAnalysis },
  { Cmd::ANALYSIS, "rms2d", Analysis_Rms2d::Alloc, Analysis_Rms2d::Help, AddAnalysis },
  { Cmd::ANALYSIS, "rmsavgcorr", Analysis_RmsAvgCorr::Alloc, Analysis_RmsAvgCorr::Help, AddAnalysis },
  { Cmd::ANALYSIS, "rotdif", Analysis_Rotdif::Alloc, Analysis_Rotdif::Help, AddAnalysis },
  { Cmd::ANALYSIS, "runningavg", Analysis_RunningAvg::Alloc, Analysis_RunningAvg::Help, AddAnalysis },
  { Cmd::ANALYSIS, "spline", Analysis_Spline::Alloc, Analysis_Spline::Help, AddAnalysis },
  { Cmd::ANALYSIS, "stat", Analysis_Statistics::Alloc, Analysis_Statistics::Help, AddAnalysis },
  { Cmd::ANALYSIS, "statistics", Analysis_Statistics::Alloc, Analysis_Statistics::Help, AddAnalysis },
  { Cmd::ANALYSIS, "ti", Analysis_TI::Alloc, Analysis_TI::Help, AddAnalysis },
  { Cmd::ANALYSIS, "timecorr", Analysis_Timecorr::Alloc, Analysis_Timecorr::Help, AddAnalysis },
  { Cmd::ANALYSIS, "vectormath", Analysis_VectorMath::Alloc, Analysis_VectorMath::Help, AddAnalysis },
  { Cmd::ANALYSIS, "wavelet", Analysis_Wavelet::Alloc, Analysis_Wavelet::Help, AddAnalysis },
  // DEPRECATED COMMANDS
  { Cmd::DEPRECATED, "acceptor",     0, Deprecate_Hbond,        0 },
  { Cmd::DEPRECATED, "avgcoord",     0, Deprecate_AvgCoord,     0 },
  { Cmd::DEPRECATED, "bondsearch",   0, Deprecate_TopSearch,    0 },
  { Cmd::DEPRECATED, "donor",        0, Deprecate_Hbond,        0 },
  { Cmd::DEPRECATED, "maxdist",      0, Deprecate_MinDist,      0 },
  { Cmd::DEPRECATED, "mindist",      0, Deprecate_MinDist,      0 },
  { Cmd::DEPRECATED, "molsearch",    0, Deprecate_TopSearch,    0 },
  { Cmd::DEPRECATED, "nobondsearch", 0, Deprecate_TopSearch,    0 },
  { Cmd::DEPRECATED, "nomolsearch",  0, Deprecate_TopSearch,    0 },
  { Cmd::DEPRECATED, "parmbondinfo", 0, Deprecate_ParmBondInfo, 0 },
  { Cmd::DEPRECATED, "parmmolinfo",  0, Deprecate_ParmMolInfo,  0 },
  { Cmd::DEPRECATED, "parmresinfo",  0, Deprecate_ParmResInfo,  0 },
  { Cmd::NONE      , 0,              0, 0,                      0 }
};
