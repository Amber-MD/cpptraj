#include "ForLoop_inData.h"
#include "ArgList.h"
#include "CpptrajState.h"
#include "CpptrajStdio.h"
#include "DataSet.h"
#include "DataSet_1D.h"
#include "DataSet_integer.h"
#include "DataSet_string.h"
#include "DataSet_unsignedInt.h"
#include "StringRoutines.h"

/** CONSTRUCTOR */
ForLoop_inData::ForLoop_inData() : set_(0), idx_(0) {}

/** Help text */
void ForLoop_inData::helpText() {
  mprintf("\t<var> indata <data set name>\n"
          "  Loop over elements of given data set.\n");
}

/** \return true if set type is valid for 'indata' loop. */
static inline bool valid_inData_set(DataSet const* ds) {
  if (ds->Group() == DataSet::SCALAR_1D ||
      ds->Type() == DataSet::STRING
     )
  {
    return true;
  }
  return false;
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
  // Check that we have the right kind of set.
  if (!valid_inData_set(set_)) {
    mprinterr("Error: Set '%s' is not a valid set type for 'for indata' loop.\n", set_->legend());
    return 1;
  }
  // Variable name.
  if (SetupLoopVar( State.DSL(), argIn.GetStringNext() )) return 1;
  // Description
  std::string description( "(" + VarName() + " indata " + set_->Meta().PrintName() + ")" );
  SetDescription( description );
  return 0;
}

/** Begin an 'indata' for loop. */
int ForLoop_inData::BeginFor(DataSetList const& CurrentVars) {
  idx_ = 0;
  return (int)set_->Size();
}

/** Check for end condition / update variable with current value. */
bool ForLoop_inData::EndFor(DataSetList& DSL) {
  if (idx_ < set_->Size()) {
    // Not yet done
    if ( set_->Type() == DataSet::STRING ) {
      // ----- String ------------------
      DataSet_string const& strSet = static_cast<DataSet_string const&>( *set_ );
      DSL.UpdateStringVar( VarName(), strSet[idx_] );
    } else {
      // ----- DataSet_1D --------------
      if ( set_->Type() == DataSet::INTEGER ) {
        // Integer numbers
        DataSet_integer const& iset = static_cast<DataSet_integer const&>( *set_ );
        DSL.UpdateStringVar( VarName(), integerToString( iset[idx_] ) );
      } else if ( set_->Type() == DataSet::UNSIGNED_INTEGER ) {
        // Unsigned integer numbers
        DataSet_unsignedInt const& uset = static_cast<DataSet_unsignedInt const&>( *set_ );
        DSL.UpdateStringVar( VarName(), integerToString( uset[idx_] ) );
      } else {
        // Some form of floating point numbers
        DataSet_1D const& set1d = static_cast<DataSet_1D const&>( *set_ );
        DSL.UpdateStringVar( VarName(), doubleToString( set1d.Dval(idx_) ) );
      }
    }
    // Increment
    idx_++;
    return false;
  }
  // Done
  return true;
}
