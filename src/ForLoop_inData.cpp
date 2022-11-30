#include "ForLoop_inData.h"
#include "ArgList.h"
#include "CpptrajState.h"
#include "CpptrajStdio.h"
#include "DataSet.h"

/** CONSTRUCTOR */
ForLoop_inData::ForLoop_inData() : set_(0) {}

/** Help text */
void ForLoop_inData::helpText() {
  mprintf("\t<var> indata <data set name>\n"
          "  Loop over elements of given data set.\n");
}

/** Set up 'indata' for loop.
  * <var> indata <set selection>
  */
int ForLoop_inData::SetupFor(CpptrajState& State, ArgList& argIn) {
  std::string setArg = argIn.GetStringKey("indata");
  if (setArg.empty()) {
    mprinterr("Error: 'for indata': missing ' indata <data set name>'.\n");
    return 1;
  }
  set_ = State.DSL().GetDataSet( setArg );
  if (set_ == 0) {
    mprinterr("Error: No data set selected by '%s'\n", setArg.c_str());
    return 1;
  }
  // Variable name.
  if (SetupLoopVar( State.DSL(), argIn.GetStringNext() )) return 1;
  // Description
  std::string description( "(" + VarName() + " indata " + set_->Meta().PrintName() + ")" );
  SetDescription( description );
  return 0;
}
