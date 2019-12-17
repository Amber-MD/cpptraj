#include <algorithm> // std::min
#include "Control.h"
#include "CpptrajStdio.h"
#include "VariableArray.h"
#include "ArgList.h"
#include "StringRoutines.h"

void Control_Set::Help() const {
  mprintf("\t{ <variable> <OP> <value> |\n"
          "\t  <variable> <OP> {atoms|residues|molecules} inmask <mask>\n"
          "\t    [%s]\n"
          "\t  <variable> <OP> trajinframes }\n",
          DataSetList::TopIdxArgs);
  mprintf("  Set (<OP> = '=') or append (<OP> = '+=') a script variable.\n"
          "  - Set script variable <variable> to value <value>.\n"
          "  - Set script variable to the number of atoms/residues/molecules in\n"
          "     the given atom mask.\n"
          "  - Set script variable to the current number of frames that will\n"
          "     be read from all previous 'trajin' statements.\n");
}

/** Set up variable with value. In this case allow any amount of whitespace,
  * so re-tokenize the original argument line (minus the command).
  */
CpptrajState::RetType
  Control_Set::SetupControl(CpptrajState& State, ArgList& argIn, VariableArray& CurrentVars)
{
  ArgList remaining = argIn.RemainingArgs();
  size_t pos0 = remaining.ArgLineStr().find_first_of("=");
  if (pos0 == std::string::npos) {
    mprinterr("Error: Expected <var>=<value>\n");
    return CpptrajState::ERR;
  }
  size_t pos1 = pos0;
  bool append = false;
  if (pos0 > 0 && remaining.ArgLineStr()[pos0-1] == '+') {
    pos0--;
    append = true;
  }
  std::string variable = NoWhitespace( remaining.ArgLineStr().substr(0, pos0) );
  if (variable.empty()) {
    mprinterr("Error: No variable name.\n");
    return CpptrajState::ERR;
  }
  ArgList equals( NoLeadingWhitespace(remaining.ArgLineStr().substr(pos1+1)) );
  std::string value;
  if (equals.Contains("inmask")) {
    AtomMask mask( equals.GetStringKey("inmask") );
    Topology* top = State.DSL().GetTopByIndex( equals );
    if (top == 0) return CpptrajState::ERR;
    if (top->SetupIntegerMask( mask )) return CpptrajState::ERR;
    if (equals.hasKey("atoms"))
      value = integerToString( mask.Nselected() );
    else if (equals.hasKey("residues")) {
      int curRes = -1;
      int nres = 0;
      for (AtomMask::const_iterator at = mask.begin(); at != mask.end(); ++at) {
        int res = (*top)[*at].ResNum();
        if (res != curRes) {
          nres++;
          curRes = res;
        }
      }
      value = integerToString( nres );
    } else if (equals.hasKey("molecules")) {
      int curMol = -1;
      int nmol = 0;
      for (AtomMask::const_iterator at = mask.begin(); at != mask.end(); ++at) {
        int mol = (*top)[*at].MolNum();
        if (mol != curMol) {
          nmol++;
          curMol = mol;
        }
      }
      value = integerToString( nmol );
    } else {
      mprinterr("Error: Expected 'atoms', 'residues', or 'molecules'.\n");
      return CpptrajState::ERR;
    }
  } else if (equals.hasKey("trajinframes")) {
    value = integerToString(State.InputTrajList().MaxFrames());
  } else
    value = equals.ArgLineStr();
  if (append)
    CurrentVars.AppendVariable( "$" + variable, value );
  else
    CurrentVars.UpdateVariable( "$" + variable, value );
  mprintf("\tVariable '%s' set to '%s'\n", variable.c_str(), value.c_str());
  for (int iarg = 0; iarg < argIn.Nargs(); iarg++)
    argIn.MarkArg( iarg );
  return CpptrajState::OK;
}

// -----------------------------------------------------------------------------
void Control_Show::Help() const {
  mprintf("  Show all current script variables and their values.\n");
}

CpptrajState::RetType
  Control_Show::SetupControl(CpptrajState& State, ArgList& argIn, VariableArray& CurrentVars)
{
  for (VariableArray::const_iterator it = CurrentVars.begin(); it != CurrentVars.end(); ++it)
    mprintf("\t%s = '%s'\n", it->first.c_str(), it->second.c_str());
  return CpptrajState::OK;
}
