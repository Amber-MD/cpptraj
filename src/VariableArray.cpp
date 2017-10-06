#include "VariableArray.h"
#include "CpptrajStdio.h"

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
ArgList VariableArray::ReplaceVariables(ArgList const& argIn) {
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
        mprinterr("Error: Unrecognized variable in command: %s\n", var_in_arg.c_str());
        return ArgList();
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
