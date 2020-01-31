#include "ForLoop.h"
#include "DataSetList.h"
#include "CpptrajStdio.h"
#include "DataSet.h"

int ForLoop::SetupLoopVar(DataSetList& DSL, std::string const& varname)
{
  loopvar_ = DSL.CheckForSet( varname );
  if (loopvar_ != 0) {
    if (loopvar_->Type() != DataSet::STRINGVAR) {
      mprinterr("Error: Set '%s' exists but is not a string variable, cannot be used for loop.\n",
                varname.c_str());
      return 1;
    }
  } else {
    loopvar_ = DSL.AddSet( DataSet::STRINGVAR, varname );
    if (loopvar_ == 0) return 1;
  }
  //mprintf("DEBUG: Loop variable: '%s'\n", loopvar_->legend());
  return 0;
}

std::string const& ForLoop::VarName() const {
  return loopvar_->Meta().Name();
}
