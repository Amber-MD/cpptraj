#include <cstdarg>
#include <algorithm> // std::sort()
#include "Command.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h" // ProcessInput()
#include "CmdInput.h"     // ProcessInput()
#include "Exec.h"
// ----- GENERAL ---------------------------------------------------------------
#include "Exec_Help.h"
#include "Exec_ReadInput.h"
// ----- COORDS ----------------------------------------------------------------
#include "Exec_CrdAction.h"
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
  Command::AddCmd( new Exec_Help(),      Cmd::EXE, 1, "help" );
  Command::AddCmd( new Exec_ReadInput(), Cmd::EXE, 1, "readinput" );

  Command::AddCmd( new Exec_CrdAction(), Cmd::EXE, 1, "crdaction" );

  Command::AddCmd( new Action_Angle(),       Cmd::ACT, 1, "angle" );
  Command::AddCmd( new Action_AreaPerMol(),  Cmd::ACT, 1, "areapermol" );
  Command::AddCmd( new Action_AtomicCorr(),  Cmd::ACT, 1, "atomiccorr" );
  Command::AddCmd( new Action_AtomicFluct(), Cmd::ACT, 2, "atomicfluct", "rmsf" );
  Command::AddCmd( new Action_AtomMap(),     Cmd::ACT, 1, "atommap" );
  Command::AddCmd( new Action_AutoImage(),   Cmd::ACT, 1, "autoimage" );
  Command::AddCmd( new Action_Average(),     Cmd::ACT, 1, "average" );
  Command::AddCmd( new Action_Bounds(),      Cmd::ACT, 1, "bounds" );
  Command::AddCmd( new Action_Box(),         Cmd::ACT, 1, "box" );
  Command::AddCmd( new Action_Center(),      Cmd::ACT, 1, "center" );
  Command::AddCmd( new Action_Channel(),     Cmd::ACT, 1, "channel" );
  Command::AddCmd( new Action_CreateCrd(),   Cmd::ACT, 1, "createcrd" );
  Command::AddCmd( new Action_Unstrip(),     Cmd::ACT, 1, "unstrip" );

  Command::AddCmd( new Analysis_AmdBias(), Cmd::ANA, 1, "amdbias" );

  // Add null ptr to indicate end of command key addresses for ReadLine
  names_.push_back( 0 );
}

void Command::Free() {
  commands_.Clear();
}

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
