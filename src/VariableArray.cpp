#include "VariableArray.h"
#include "CpptrajStdio.h"
#include "DataSet_string.h"
#include "DataSet_1D.h"
#include "StringRoutines.h"

/** Add variable/value to array if it doesnt exist, otherwise set value. */
void VariableArray::UpdateVariable(std::string const& varname, std::string const& value)
{
  Varray::iterator it = CurrentVars_.begin();
  for (; it != CurrentVars_.end(); ++it)
    if (it->first == varname)
      break;
  if (it == CurrentVars_.end())
    CurrentVars_.push_back( Vpair(varname, value) );
  else
    it->second = value;
}

/** Add variable/value to array if it doesnt exist, otherwise append value. */
void VariableArray::AppendVariable(std::string const& varname, std::string const& value)
{
  Varray::iterator it = CurrentVars_.begin();
  for (; it != CurrentVars_.end(); ++it)
    if (it->first == varname)
      break;
  if (it == CurrentVars_.end())
    CurrentVars_.push_back( Vpair(varname, value) );
  else
    it->second += value;
}

/** Replace all variables in given ArgList with their values. */
ArgList VariableArray::ReplaceVariables(ArgList const& argIn, DataSetList const& DSL) {
  ArgList modCmd = argIn;
  for (int n = 0; n < modCmd.Nargs(); n++) {
    size_t pos = modCmd[n].find("$");
    if (pos != std::string::npos) {
      // Argument is/contains a variable. Find first non-alphanumeric char
      size_t len = 1;
      for (size_t pos1 = pos+1; pos1 < modCmd[n].size(); pos1++, len++)
        if (!isalnum(modCmd[n][pos1])) break;
      std::string var_in_arg = modCmd[n].substr(pos, len);
      // See if variable occurs in CurrentVars_
      Varray::const_iterator vp = CurrentVars_.begin();
      for (; vp != CurrentVars_.end(); ++vp)
        if (vp->first == var_in_arg) break;
      // If found replace with value from CurrentVars_
      if (vp != CurrentVars_.end()) {
        std::string arg = modCmd[n];
        arg.replace(pos, vp->first.size(), vp->second);
        modCmd.ChangeArg(n, arg);
      } else {
        // Not found in CurrentVars_; see if this is a DataSet.
        for (size_t pos1 = pos+len; pos1 < modCmd[n].size(); pos1++, len++)
          if (!isalnum(modCmd[n][pos1]) &&
              modCmd[n][pos1] != '[' &&
              modCmd[n][pos1] != ':' &&
              modCmd[n][pos1] != ']' &&
              modCmd[n][pos1] != '_' &&
              modCmd[n][pos1] != '-' &&
              modCmd[n][pos1] != '%')
            break;
        var_in_arg = modCmd[n].substr(pos+1, len-1);
        DataSet* ds = DSL.GetDataSet( var_in_arg );
        if (ds == 0) {
          mprinterr("Error: Unrecognized variable in command: %s\n", var_in_arg.c_str());
          return ArgList();
        } else {
          if (ds->Type() != DataSet::STRING && ds->Group() != DataSet::SCALAR_1D) {
            mprinterr("Error: Only 1D data sets supported.\n");
            return ArgList();
          }
          if (ds->Size() < 1) {
            mprinterr("Error: Set is empty.\n");
            return ArgList();
          }
          if (ds->Size() > 1)
            mprintf("Warning: Only using first value.\n");
          std::string value;
          if (ds->Type() == DataSet::STRING)
            value = (*((DataSet_string*)ds))[0];
          else
            value = doubleToString(((DataSet_1D*)ds)->Dval(0));
          mprintf("\tReplaced variable '$%s' with value '%s' from DataSet '%s'\n",
                  var_in_arg.c_str(), value.c_str(), ds->legend());
          std::string arg = modCmd[n];
          arg.replace(pos, var_in_arg.size()+1, value);
          modCmd.ChangeArg(n, arg);
        }
      }
    }
  }
  return modCmd;
}

void VariableArray::PrintVariables() const {
  for (Varray::const_iterator vp = CurrentVars_.begin(); vp != CurrentVars_.end(); ++vp)
    mprintf(" %s=%s", vp->first.c_str(), vp->second.c_str());
  mprintf("\n");
}
