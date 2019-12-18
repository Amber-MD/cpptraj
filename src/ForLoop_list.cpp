#include "ForLoop_list.h"
#include "CpptrajStdio.h"
#include "ArgList.h"
#include "FileName.h"
#include "VariableArray.h"

int ForLoop_list::SetupFor(CpptrajState& State, std::string const& expr, ArgList& argIn) {
  int Niterations = -1;
  // <var> in <string0>[,<string1>...]
  //MH.varType_ = ftype;
  argIn.PrintDebug();
  // Variable name. Should be expr.
  SetVarName( expr );
  // Comma-separated list of strings.
  std::string listArg = argIn.GetStringKey("in");
  if (listArg.empty()) {
    mprinterr("Error: 'for in': missing ' in <comma-separated list of strings>'.\n");
    return 1;
  }
  ArgList list(listArg, ",");
  if (list.Nargs() < 1) {
    mprinterr("Error: Could not parse '%s' for 'for in'\n", listArg.c_str());
    return 1;
  }
  for (int il = 0; il != list.Nargs(); il++) {
    // Check if file name expansion should occur
    if (list[il].find_first_of("*?") != std::string::npos) {
      File::NameArray files = File::ExpandToFilenames( list[il] );
      for (File::NameArray::const_iterator fn = files.begin(); fn != files.end(); ++fn)
        List_.push_back( fn->Full() );
    } else
      List_.push_back( list[il] );
  }
  Niterations = (int)List_.size();
  // Description
  std::string description( VarName() + " in " + listArg );
  SetDescription( description );
  SetNiterations( Niterations );
  return 0;
}

int ForLoop_list::BeginFor() {
  sdx_ = List_.begin();
  return 0;
}

bool ForLoop_list::EndFor(VariableArray& CurrentVars) {
  if (sdx_ == List_.end()) return true;
  // Get variable value
  CurrentVars.UpdateVariable( VarName(), *(sdx_) );
  // Increment
  ++(sdx_);
  return false;
}
