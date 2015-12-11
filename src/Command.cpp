#include <cstdarg>
#include <algorithm> // std::sort()
#include "Command.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h" // ProcessInput()
#include "CmdInput.h"     // ProcessInput()
#include "Exec.h"
// ----- GENERAL ---------------------------------------------------------------
#include "Exec_Commands.h"
#include "Exec_DataFile.h"
#include "Exec_DataFilter.h"
#include "Exec_Help.h"
#include "Exec_ReadInput.h"
// ----- SYSTEM ----------------------------------------------------------------
#include "Exec_System.h"
// ----- COORDS ----------------------------------------------------------------
#include "Exec_CombineCoords.h"
#include "Exec_CrdAction.h"
#include "Exec_CrdOut.h"
#include "Exec_LoadCrd.h"
#include "Exec_LoadTraj.h"
// ----- ACTION ----------------------------------------------------------------
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
// ----- ANALYSIS --------------------------------------------------------------
#include "Analysis_AmdBias.h"

CmdList Command::commands_ = CmdList();

const Cmd Command::EMPTY_ = Cmd();

Command::Carray Command::names_ = Command::Carray();

/** Initialize all commands. */
void Command::Init() {
  // GENERAL
  Command::AddCmd( new Exec_ActiveRef(),     Cmd::EXE, 1, "activeref" );
  Command::AddCmd( new Exec_Clear(),         Cmd::EXE, 1, "clear" );
  Command::AddCmd( new Exec_CreateDataFile(),Cmd::EXE, 1, "create" );
  Command::AddCmd( new Exec_DataFileCmd(),   Cmd::EXE, 1, "datafile" );
  Command::AddCmd( new Exec_DataFilter(),    Cmd::EXE, 1, "datafilter" );
  Command::AddCmd( new Exec_Help(),          Cmd::EXE, 1, "help" );
  Command::AddCmd( new Exec_ListAll(),       Cmd::EXE, 1, "list" );
  Command::AddCmd( new Exec_NoExitOnError(), Cmd::EXE, 1, "noexitonerror" );
  Command::AddCmd( new Exec_NoProgress(),    Cmd::EXE, 1, "noprogress" );
  Command::AddCmd( new Exec_Quit(),          Cmd::EXE, 2, "exit", "quit" );
  Command::AddCmd( new Exec_ReadInput(),     Cmd::EXE, 1, "readinput" );
  Command::AddCmd( new Exec_RemoveData(),    Cmd::EXE, 1, "removedata" );
  Command::AddCmd( new Exec_Run(),           Cmd::EXE, 2, "go", "run" );
  Command::AddCmd( new Exec_SilenceActions(),Cmd::EXE, 1, "silenceactions" );
  Command::AddCmd( new Exec_SetListDebug(),  Cmd::EXE, 2, "debug", "prnlev" );
  Command::AddCmd( new Exec_WriteDataFile(), Cmd::EXE, 2, "write", "writedata" );
  // SYSTEM
  Command::AddCmd( new Exec_System(), Cmd::EXE, 6, "gnuplot", "head", "less", "ls", "pwd", "xmgrace" );
  // COORDS
  Command::AddCmd( new Exec_CombineCoords(),Cmd::EXE, 1, "combinecrd" ); 
  Command::AddCmd( new Exec_CrdAction(),    Cmd::EXE, 1, "crdaction" );
  Command::AddCmd( new Exec_CrdOut(),       Cmd::EXE, 1, "crdout" );
  Command::AddCmd( new Exec_LoadCrd(),      Cmd::EXE, 1, "loadcrd" );
  Command::AddCmd( new Exec_LoadTraj(),     Cmd::EXE, 1, "loadtraj" );
  // ACTION
  Command::AddCmd( new Action_Angle(),         Cmd::ACT, 1, "angle" );
  Command::AddCmd( new Action_AreaPerMol(),    Cmd::ACT, 1, "areapermol" );
  Command::AddCmd( new Action_AtomicCorr(),    Cmd::ACT, 1, "atomiccorr" );
  Command::AddCmd( new Action_AtomicFluct(),   Cmd::ACT, 2, "atomicfluct", "rmsf" );
  Command::AddCmd( new Action_AtomMap(),       Cmd::ACT, 1, "atommap" );
  Command::AddCmd( new Action_AutoImage(),     Cmd::ACT, 1, "autoimage" );
  Command::AddCmd( new Action_Average(),       Cmd::ACT, 1, "average" );
  Command::AddCmd( new Action_Bounds(),        Cmd::ACT, 1, "bounds" );
  Command::AddCmd( new Action_Box(),           Cmd::ACT, 1, "box" );
  Command::AddCmd( new Action_Center(),        Cmd::ACT, 1, "center" );
  Command::AddCmd( new Action_Channel(),       Cmd::ACT, 1, "channel" ); // HIDDEN
  Command::AddCmd( new Action_CheckStructure(),Cmd::ACT, 3,"check","checkoverlap","checkstructure");
  Command::AddCmd( new Action_CheckChirality(),Cmd::ACT, 1, "checkchirality" );
  Command::AddCmd( new Action_Closest(),       Cmd::ACT, 2, "closest", "closestwaters" );
  Command::AddCmd( new Action_ClusterDihedral(),Cmd::ACT,1, "clusterdihedral" );
  Command::AddCmd( new Action_Contacts(),      Cmd::ACT, 1, "contacts" );
  Command::AddCmd( new Action_CreateCrd(),     Cmd::ACT, 1, "createcrd" );
  Command::AddCmd( new Action_CreateReservoir(),Cmd::ACT,1, "createreservoir" );
  Command::AddCmd( new Action_Density(),       Cmd::ACT, 1, "density" );
  Command::AddCmd( new Action_Diffusion(),     Cmd::ACT, 1, "diffusion" );
  Command::AddCmd( new Action_Dihedral(),      Cmd::ACT, 1, "dihedral" );
  Command::AddCmd( new Action_DihedralScan(),  Cmd::ACT, 1, "dihedralscan" );
  Command::AddCmd( new Action_Dipole(),        Cmd::ACT, 1, "dipole" );
  Command::AddCmd( new Action_Distance(),      Cmd::ACT, 1, "distance" );
  Command::AddCmd( new Action_DNAionTracker(), Cmd::ACT, 1, "dnaiontracker" ); // HIDDEN
  Command::AddCmd( new Action_DistRmsd(),      Cmd::ACT, 2, "drms", "drmsd" );
  Command::AddCmd( new Action_DSSP(),          Cmd::ACT, 2, "dssp", "secstruct" );
  Command::AddCmd( new Action_Energy(),        Cmd::ACT, 1, "energy" );
  Command::AddCmd( new Action_FilterByData(),  Cmd::ACT, 1, "filter" );
  Command::AddCmd( new Action_FixAtomOrder(),  Cmd::ACT, 1, "fixatomorder" );
  Command::AddCmd( new Action_Gist(),          Cmd::ACT, 1, "gist" );
  Command::AddCmd( new Action_GridFreeEnergy(),Cmd::ACT, 1, "gfe" ); // HIDDEN
  Command::AddCmd( new Action_Grid(),          Cmd::ACT, 1, "grid" );
  Command::AddCmd( new Action_Hbond(),         Cmd::ACT, 1, "hbond" );
  Command::AddCmd( new Action_Image(),         Cmd::ACT, 1, "image" );
  Command::AddCmd( new Action_Jcoupling(),     Cmd::ACT, 1, "jcoupling" );
  Command::AddCmd( new Action_LESsplit(),      Cmd::ACT, 1, "lessplit" );
  Command::AddCmd( new Action_LIE(),           Cmd::ACT, 1, "lie" );
  Command::AddCmd( new Action_OrderParameter(),Cmd::ACT, 1, "lipidorder" );
  Command::AddCmd( new Action_MakeStructure(), Cmd::ACT, 1, "makestructure" );
  Command::AddCmd( new Action_Mask(),          Cmd::ACT, 1, "mask" );
  Command::AddCmd( new Action_Matrix(),        Cmd::ACT, 1, "matrix" );
  Command::AddCmd( new Action_MinImage(),      Cmd::ACT, 1, "minimage" );
  Command::AddCmd( new Action_Molsurf(),       Cmd::ACT, 1, "molsurf" );
  Command::AddCmd( new Action_MultiDihedral(), Cmd::ACT, 1, "multidihedral" );
  Command::AddCmd( new Action_MultiVector(),   Cmd::ACT, 1, "multivector" );
  Command::AddCmd( new Action_NAstruct(),      Cmd::ACT, 1, "nastruct" );
  Command::AddCmd( new Action_NativeContacts(),Cmd::ACT, 1, "nativecontacts" );
  Command::AddCmd( new Action_NMRrst(),        Cmd::ACT, 1, "nmrrst" );
  Command::AddCmd( new Action_Outtraj(),       Cmd::ACT, 1, "outtraj" );
  Command::AddCmd( new Action_PairDist(),      Cmd::ACT, 1, "pairdist" );
  Command::AddCmd( new Action_Pairwise(),      Cmd::ACT, 1, "pairwise" );
  Command::AddCmd( new Action_Principal(),     Cmd::ACT, 1, "principal" );
  Command::AddCmd( new Action_Projection(),    Cmd::ACT, 1, "projection" );
  Command::AddCmd( new Action_Pucker(),        Cmd::ACT, 1, "pucker" );
  Command::AddCmd( new Action_Radgyr(),        Cmd::ACT, 2, "radgyr", "rog" );
  Command::AddCmd( new Action_Radial(),        Cmd::ACT, 2, "radial", "rdf" );
  Command::AddCmd( new Action_RandomizeIons(), Cmd::ACT, 1, "randomizeions" );
  Command::AddCmd( new Action_ReplicateCell(), Cmd::ACT, 1, "replicatecell" );
  Command::AddCmd( new Action_Rmsd(),          Cmd::ACT, 2, "rms", "rmsd" );
  Command::AddCmd( new Action_Rotate(),        Cmd::ACT, 1, "rotate" );
  Command::AddCmd( new Action_RunningAvg(),    Cmd::ACT, 2, "runavg", "runningaverage" );
  Command::AddCmd( new Action_Scale(),         Cmd::ACT, 1, "scale" );
  Command::AddCmd( new Action_SetVelocity(),   Cmd::ACT, 1, "setvelocity" );
  Command::AddCmd( new Action_Spam(),          Cmd::ACT, 1, "spam" );
  Command::AddCmd( new Action_STFC_Diffusion(),Cmd::ACT, 1, "stfcdiffusion" );
  Command::AddCmd( new Action_Strip(),         Cmd::ACT, 1, "strip" );
  Command::AddCmd( new Action_Surf(),          Cmd::ACT, 1, "surf" );
  Command::AddCmd( new Action_SymmetricRmsd(), Cmd::ACT, 1, "symmrmsd" );
  Command::AddCmd( new Action_Temperature(),   Cmd::ACT, 1, "temperature" );
  Command::AddCmd( new Action_Translate(),     Cmd::ACT, 2, "trans", "translate" );
  Command::AddCmd( new Action_Unstrip(),       Cmd::ACT, 1, "unstrip" );
  Command::AddCmd( new Action_Unwrap(),        Cmd::ACT, 1, "unwrap" );
  Command::AddCmd( new Action_Vector(),        Cmd::ACT, 1, "vector" );
  Command::AddCmd( new Action_VelocityAutoCorr(),Cmd::ACT,1,"velocityautocorr" );
  Command::AddCmd( new Action_Volmap(),        Cmd::ACT, 1, "volmap" );
  Command::AddCmd( new Action_Volume(),        Cmd::ACT, 1, "volume" );
  Command::AddCmd( new Action_Watershell(),    Cmd::ACT, 1, "watershell" );
  // ANALYSIS
  Command::AddCmd( new Analysis_AmdBias(), Cmd::ANA, 1, "amdbias" );

  // Add null ptr to indicate end of command key addresses for ReadLine
  names_.push_back( 0 );
}

/** Free all commands. Should only be called just before program exit. */
void Command::Free() { commands_.Clear(); }

/** \param oIn Pointer to DispatchObject to add as command.
  * \param dIn Command destination
  * \param nKeys Number of command keywords associated with this command.
  * The remaining arguments are the nKeys command keywords.
  */
void Command::AddCmd(DispatchObject* oIn, Cmd::DestType dIn, int nKeys, ...) {
  Cmd::Sarray keys;
  va_list args;
  va_start(args, nKeys);
  for (int nk = 0; nk < nKeys; nk++) {
    char* key = va_arg(args, char*);
    keys.push_back( std::string(key) );
  }
  va_end(args);
  commands_.Add( Cmd(oIn, keys, dIn) );
  // Store memory addresses of command keys for ReadLine
  for (Cmd::key_iterator key = commands_.Back().keysBegin();
                         key != commands_.Back().keysEnd(); ++key)
    names_.push_back( key->c_str() );
}

/// Warn about deprecated commands.
void Command::WarnDeprecated(const char* keyIn, Cmd const& cmdIn)
{
  mprinterr("Error: '%s' is deprecated.\n", keyIn);
  cmdIn.Help();
}

/** Search Commands list for command with given keyword and object type. */
Cmd const& Command::SearchTokenType(DispatchObject::Otype catIn, const char* keyIn)
{
  for (CmdList::const_iterator cmd = commands_.begin(); cmd != commands_.end(); ++cmd)
  {
    bool match = cmd->KeyMatches(keyIn);
    if (cmd->Obj().Type() == DispatchObject::DEPRECATED && match) {
      WarnDeprecated( keyIn, *cmd );
      return EMPTY_;
    }
    if (catIn != cmd->Obj().Type()) continue;
    if (match) return *cmd;
  }
  mprinterr("'%s': Command not found.\n", keyIn);
  return EMPTY_;
}

/** Search the Commands list for given command.
  * \return the token if found, 0 if not.
  */
Cmd const& Command::SearchToken(ArgList& argIn) {
  // SPECIAL CASE: For backwards compat. remove analyze prefix
  if (argIn.CommandIs("analyze")) {
    argIn.RemoveFirstArg();
    argIn.MarkArg(0); // Mark new first arg as command
    return (SearchTokenType(DispatchObject::ANALYSIS, argIn.Command()));
  }
  // Search for command.
  for (CmdList::const_iterator cmd = commands_.begin(); cmd != commands_.end(); ++cmd)
  {
    if ( cmd->KeyMatches( argIn.Command() ) ) {
      if ( cmd->Obj().Type() == DispatchObject::DEPRECATED) {
        WarnDeprecated( argIn.Command(), *cmd );
        return EMPTY_;
      } else
        return *cmd;
    }
  }
  //mprinterr("'%s': Command not found.\n", argIn.Command());
  return EMPTY_;
}

/** First list the command category, then the commands for that category
  * in alphabetical order. Should not be called with NONE, HIDDEN, or
  * DEPRECATED.
  */
void Command::ListCommandsForType(DispatchObject::Otype typeIn) {
  std::vector< std::string > command_keys;
  mprintf("%s Commands:\n", DispatchObject::ObjKeyword(typeIn));
  for (CmdList::const_iterator cmd = commands_.begin(); cmd != commands_.end(); ++cmd)
  {
    if (cmd->Obj().Type() == typeIn)
      for (Cmd::key_iterator key = cmd->keysBegin(); key != cmd->keysEnd(); ++key)
        command_keys.push_back( *key );
  }
  std::sort( command_keys.begin(), command_keys.end() );
  std::string Line("        ");
  for (std::vector< std::string >::const_iterator key = command_keys.begin();
                                                  key != command_keys.end(); ++key)
  {
    if ( Line.size() + key->size() + 1 > 80 ) {
      mprintf("%s\n", Line.c_str());
      Line.assign("        ");
    }
    Line.append( *key + " " );
  }
  if (!Line.empty()) // TODO is it ever empty?
    mprintf("%s\n", Line.c_str());
}
    
/** List all commands of the given type, or all commands if type
  * is DispatchObject::NONE.
  */
void Command::ListCommands(DispatchObject::Otype typeIn) {
  if (typeIn == DispatchObject::NONE) {
    for (int idx = 1; idx != DispatchObject::HIDDEN; idx++)
      ListCommandsForType( (DispatchObject::Otype)idx );
  } else
    ListCommandsForType( typeIn );
}

/** Search for the given command and execute it. EXE commands are executed
  * immediately and then freed. ACT and ANA commands are sent to the
  * CpptrajState for later execution.
  */
CpptrajState::RetType Command::Dispatch(CpptrajState& State, std::string const& commandIn)
{
  ArgList cmdArg( commandIn );
  cmdArg.MarkArg(0); // Always mark the first arg as the command 
  Cmd const& cmd = SearchToken( cmdArg );
/*  Cmd::RetType ret_val = Cmd::OK;
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
*/
  CpptrajState::RetType ret_val = CpptrajState::OK;
  if (cmd.Empty()) {
    mprinterr("Command '%s' not found.\n", cmdArg.Command());
    ret_val = CpptrajState::ERR;
  } else {
    DispatchObject* obj = cmd.Alloc();
    switch (cmd.Destination()) {
      case Cmd::EXE:
        ret_val = ((Exec*)obj)->Execute( State, cmdArg );
        delete obj;
        break;
      case Cmd::ACT: State.AddToActionQueue( (Action*)obj, cmdArg ); break;
      case Cmd::ANA: State.AddToAnalysisQueue( (Analysis*)obj, cmdArg ); break; 
    }
  }
  return ret_val;
}

CpptrajState::RetType Command::ProcessInput(CpptrajState& State, std::string const& inputFilename)
{
  BufferedLine infile;
  if (infile.OpenFileRead( inputFilename )) {
    if (!inputFilename.empty())
      mprinterr("Error: Could not open input file '%s'\n", inputFilename.c_str());
    return CpptrajState::ERR;
  }
  mprintf("INPUT: Reading input from '%s'\n", infile.Filename().full());
  // Read in each line of input.
  int nInputErrors = 0;
  CpptrajState::RetType cmode = CpptrajState::OK;
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
      if (cmode == CpptrajState::ERR) {
        nInputErrors++;
        if (State.ExitOnError()) break;
      } else if (cmode == CpptrajState::QUIT)
        break;
    }
    // Reset Input line
    input.Clear();
    ptr = infile.Line();
  }
  infile.CloseFile();
  if (nInputErrors > 0) {
    mprinterr("\t%i errors encountered reading input.\n", nInputErrors);
    return CpptrajState::ERR;
  }
  return cmode;
}
