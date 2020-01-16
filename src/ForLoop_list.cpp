#include "ForLoop_list.h"
#include "CpptrajStdio.h"
#include "ArgList.h"
#include "FileName.h"
#include "DataSetList.h"
#include "CpptrajState.h"

void ForLoop_list::helpText() {
  mprintf("\t<var> in <list>\n"
          "  Loop over values in comma-separated list. Values can be file names\n"
          "  containing wildcard characters ('*' or '?').\n");
}

int ForLoop_list::SetupFor(CpptrajState& State, ArgList& argIn) {
  // <var> in <string0>[,<string1>...]
  //MH.varType_ = ftype;
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
  // Variable name.
  if (SetupLoopVar( State.DSL(), argIn.GetStringNext() )) return 1;
  // Go through list of strings
  for (int il = 0; il != list.Nargs(); il++) {
    Names_.push_back( list[il] );
  }
  // Description
  std::string description( "(" + VarName() + " in " + listArg + ")" );
  SetDescription( description );
  return 0;
}

int ForLoop_list::BeginFor(DataSetList const& CurrentVars) {
  // Go through list of strings
  List_.clear();
  for (Sarray::const_iterator it = Names_.begin(); it != Names_.end(); ++it)
  {
    // Try variable replacement
    std::string listEntry;
    CurrentVars.ReplaceVariables(listEntry, *it);
    // Check if file name expansion should occur
    if (listEntry.find_first_of("*?") != std::string::npos) {
      File::NameArray files = File::ExpandToFilenames( listEntry );
      // DEBUG
      //for (File::NameArray::const_iterator fn = files.begin(); fn != files.end(); ++fn)
      //  mprintf("DEBUG: '%s'\n", fn->full());
      // Allow wildcard expansion to fail with a warning.
      if (!files.empty() && files.front().Full().compare( listEntry ) == 0) {
        mprintf("Warning: '%s' selects no files.\n", it->c_str());
      } else {
        for (File::NameArray::const_iterator fn = files.begin(); fn != files.end(); ++fn)
          List_.push_back( fn->Full() );
      }
    } else
      List_.push_back( listEntry );
  }

  sdx_ = List_.begin();
  return (int)List_.size();
}

bool ForLoop_list::EndFor(DataSetList& DSL) {
  if (sdx_ == List_.end()) return true;
  // Get variable value
  DSL.UpdateStringVar( VarName(), *(sdx_) );
  // Increment
  ++(sdx_);
  return false;
}
