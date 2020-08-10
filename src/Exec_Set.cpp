#include "Exec_Set.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include "DataSet_StringVar.h"

void Exec_Set::Help() const {
  mprintf("\t{ <variable> <OP> <value> |\n"
          "\t  <variable> <OP> {atoms|residues|molecules|atomnums|\n"
          "\t                   resnums|oresnums|molnums} inmask <mask>\n"
          "\t    [%s]\n"
          "\t  <variable> <OP> trajinframes }\n",
          DataSetList::TopIdxArgs);
  mprintf("  Set (<OP> = '=') or append (<OP> = '+=') a script variable.\n"
          "  - Set script variable <variable> to value <value>.\n"
          "  - Set script variable to the number of atoms/residues/molecules\n"
          "    or selected atom #s/residue #s/original residue #s/molecule #s\n"
          "    in the given mask.\n"
          "  - Set script variable to the current number of frames that will\n"
          "     be read from all previous 'trajin' statements.\n");
}

/** Set up variable with value. In this case allow any amount of whitespace,
  * so re-tokenize the original argument line (minus the command).
  */
Exec::RetType Exec_Set::Execute(CpptrajState& State, ArgList& argIn)
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
    mprintf("\tUsing topology: %s\n", top->c_str());
    if (top->SetupIntegerMask( mask )) return CpptrajState::ERR;
    if (equals.hasKey("atoms"))
      value = integerToString( mask.Nselected() );
    else if (equals.hasKey("residues")) {
      std::vector<int> resnums = top->ResnumsSelectedBy( mask );
      value = integerToString( resnums.size() );
    } else if (equals.hasKey("molecules")) {
      std::vector<int> molnums = top->MolnumsSelectedBy( mask );
      value = integerToString( molnums.size() );
    } else if (equals.hasKey("atomnums")) {
      value = ArrayToRangeExpression( mask.Selected(), 1 );
    } else if (equals.hasKey("resnums")) {
      std::vector<int> resnums = top->ResnumsSelectedBy( mask );
      value = ArrayToRangeExpression( resnums, 1 );
    } else if (equals.hasKey("molnums")) {
      std::vector<int> molnums = top->MolnumsSelectedBy( mask );
      value = ArrayToRangeExpression( molnums, 1 );
    } else if (equals.hasKey("oresnums")) {
      std::vector<int> resnums = top->ResnumsSelectedBy( mask );
      std::vector<int> oresnums;
      oresnums.reserve( resnums.size() );
      for (std::vector<int>::const_iterator it = resnums.begin(); it != resnums.end(); ++it)
        oresnums.push_back( top->Res(*it).OriginalResNum() );
      value = ArrayToRangeExpression( oresnums, 0 );
    } else {
      mprinterr("Error: Expected one of: 'atoms', 'residues', 'molecules',\n"
                "Error:   'atomnums', 'resnums', 'oresnums', or 'molnums'.\n");
      return CpptrajState::ERR;
    }
  } else if (equals.hasKey("trajinframes")) {
    value = integerToString(State.InputTrajList().MaxFrames());
  } else
    value = equals.ArgLineStr();
  if (append) {
    AppendVariable( State, variable, value );
    mprintf("\tVariable '%s' appended with '%s'\n", variable.c_str(), value.c_str());
  } else {
    UpdateVariable( State, variable, value );
    mprintf("\tVariable '%s' set to '%s'\n", variable.c_str(), value.c_str());
  }
  for (int iarg = 0; iarg < argIn.Nargs(); iarg++)
    argIn.MarkArg( iarg );
  return CpptrajState::OK;
}

/** Create new variable set with given value. */
int Exec_Set::AddVariable(CpptrajState& State, std::string const& varname,
                              std::string const& value)
{
  DataSet* ds = State.DSL().AddSet( DataSet::STRINGVAR, MetaData(varname) );
  if (ds == 0) return 1;
  DataSet_StringVar& var = static_cast<DataSet_StringVar&>( *ds );
  var.assign( value );

  return 0;
}

/** Add variable/value if it doesnt exist, otherwise append value. */
int Exec_Set::AppendVariable(CpptrajState& State, std::string const& varname,
                              std::string const& value)
{
  DataSet* ds = State.DSL().CheckForSet( varname );
  if (ds == 0) {
    // Create new
    return AddVariable(State, varname, value);
  }
  // Set exists; append
  if (ds->Type() != DataSet::STRINGVAR) {
    mprinterr("Error: 'set +=' can only be done for string variable sets.\n");
    return 1;
  }
  DataSet_StringVar& var = static_cast<DataSet_StringVar&>( *ds );
  var.append( value );
  return 0;
}

/** Add variable/value if it doesnt exist, otherwise set value. */
int Exec_Set::UpdateVariable(CpptrajState& State, std::string const& varname,
                             std::string const& value)
{
  DataSet* ds = State.DSL().CheckForSet( varname );
  if (ds == 0) {
    // Create new
    return AddVariable(State, varname, value);
  }
  // Set exists; replace current value
  if (ds->Type() != DataSet::STRINGVAR) {
    mprinterr("Error: 'set +=' can only be done for string variable sets.\n");
    return 1;
  }
  DataSet_StringVar& var = static_cast<DataSet_StringVar&>( *ds );
  var.assign( value );
  return 0;
}
