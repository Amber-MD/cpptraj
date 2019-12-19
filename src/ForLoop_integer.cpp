#include "ForLoop_integer.h"
#include "CpptrajStdio.h"
#include "ArgList.h"
#include "StringRoutines.h"
#include "VariableArray.h"

const char* ForLoop_integer::OpStr_[] = {"+=", "-=", "<", ">"};

ForLoop_integer::ForLoop_integer() :
  endOp_(NO_OP),
  incOp_(NO_OP),
  start_(0),
  end_(0),
  inc_(0),
  currentVal_(0)
{}

/** Setup integer for loop. */
int ForLoop_integer::SetupFor(CpptrajState& State, std::string const& expr, ArgList& argIn) {
  // [<var>=<start>;[<var><OP><end>;]<var><OP>[<value>]]
  //SetType( INTEGER );
  ArgList varArg( expr, ";" );
  varArg.PrintDebug(); // DEBUG
  if (varArg.Nargs() < 2 || varArg.Nargs() > 3) {
    mprinterr("Error: Malformed 'for' loop variable.\n"
              "Error: Expected '[<var>=<start>;[<var><OP><end>;]<var><OP>[<value>]]'\n"
              "Error: Got '%s'\n", expr.c_str());
    return 1;
  }
  // First argument: <var>=<start>
  ArgList startArg( varArg[0], "=" );
  if (startArg.Nargs() != 2) {
    mprinterr("Error: Malformed 'start' argument.\n"
              "Error: Expected <var>=<start>, got '%s'\n", varArg[0].c_str());
    return 1;
  }
  SetVarName( startArg[0] );
  mprintf("DEBUG: Start argument: '%s' = '%s'\n", VarName().c_str(), startArg[1].c_str());
  if ( startArg[1][0] == '$' ) {
    // Variable name
    startVarName_ = startArg[1];
  } else if (!validInteger(startArg[1])) {
    // Not Integer or variable
    mprinterr("Error: Start argument must be an integer or variable name.\n");
    return 1;
  } else
    start_ = convertToInteger(startArg[1]);
  // Second argument: <var><OP><end>
  size_t pos0 = VarName().size() - 1; // Minus 1 to account for prepended '$'
  size_t pos1 = pos0 + 1;
  endOp_ = NO_OP;
  int iargIdx = 1;
  if (varArg.Nargs() == 3) {
    iargIdx = 2;
    if ( varArg[1][pos0] == '<' )
      endOp_ = LESS_THAN;
    else if (varArg[1][pos0] == '>')
      endOp_ = GREATER_THAN;
    if (endOp_ == NO_OP) {
      mprinterr("Error: Unrecognized end op: '%s'\n",
                varArg[1].substr(pos0, pos1-pos0).c_str());
      return 1;
    }
    std::string endStr = varArg[1].substr(pos1);
    if (!validInteger(endStr)) {
      // TODO allow variables
      mprinterr("Error: End argument must be an integer.\n");
      return 1;
    } else
      end_ = convertToInteger(endStr);
  }
  // Third argument: <var><OP>[<value>]
  pos1 = pos0 + 2;
  incOp_ = NO_OP;
  bool needValue = false;
  if ( varArg[iargIdx][pos0] == '+' ) {
    if (varArg[iargIdx][pos0+1] == '+') {
      incOp_ = INCREMENT;
      inc_ = 1;
    } else if (varArg[iargIdx][pos0+1] == '=') {
      incOp_ = INCREMENT;
      needValue = true;
    }
  } else if ( varArg[iargIdx][pos0] == '-' ) {
    if (varArg[iargIdx][pos0+1] == '-' ) {
      incOp_ = DECREMENT;
      inc_ = 1;
    } else if (varArg[iargIdx][pos0+1] == '=') {
      incOp_ = DECREMENT;
      needValue = true;
    }
  }
  if (incOp_ == NO_OP) {
    mprinterr("Error: Unrecognized increment op: '%s'\n",
              varArg[iargIdx].substr(pos0, pos1-pos0).c_str());
    return 1;
  }
  if (needValue) {
    std::string incStr = varArg[iargIdx].substr(pos1);
    if (!validInteger(incStr)) {
      mprinterr("Error: increment value is not a valid integer.\n");
      return 1;
    }
    inc_ = convertToInteger(incStr);
    if (inc_ < 1) {
      mprinterr("Error: Extra '-' detected in increment.\n");
      return 1;
    }
  }
  // Description
  std::string sval;
  if (startVarName_.empty())
    sval = integerToString(start_);
  else
    sval = startVarName_;
  std::string description("(" + VarName() + "=" + sval + "; ");
  std::string eval;
  if (startVarName_.empty() && iargIdx == 2) {
    // Start value specified and End argument present
    eval = integerToString(end_);
    description.append(VarName() + std::string(OpStr_[endOp_]) + eval + "; ");

  }
  description.append( VarName() + std::string(OpStr_[incOp_]) +
                       integerToString(inc_) + ")" );
  SetDescription( description );
  // If decrementing just negate value
  if (incOp_ == DECREMENT)
    inc_ = -inc_;
  // DEBUG
  //mprintf("DEBUG: start=%i endOp=%i end=%i incOp=%i val=%i startArg=%s endArg=%s\n",
  //        start_, (int)endOp_, end_, (int)incOp_, inc_,
  //        startArg_.c_str(), endArg_.c_str());

  return 0;
}

/** Calculate the number of times the loop will execute. Also check
  * that start and end values are sane.
  */
int ForLoop_integer::calcNumIterations() const {
  // Check end > start for increment, start > end for decrement
  int maxval, minval;
  if (incOp_ == INCREMENT) {
    if (start_ >= end_) {
      mprinterr("Error: start must be less than end for increment.\n");
      return NITERATIONS_UNKNOWN;
    }
    minval = start_;
    maxval = end_;
  } else {
    if (end_ >= start_) {
      mprinterr("Error: end must be less than start for decrement.\n");
      return NITERATIONS_UNKNOWN;
    }
    minval = end_;
    maxval = start_;
  }
  // Figure out number of iterations
  int Niterations = (maxval - minval) / inc_;
  // Sanity check
  if (Niterations < 0) return NITERATIONS_UNKNOWN;
  if (((maxval-minval) % inc_) > 0) Niterations++;
  return Niterations;
}

/** Set start value and determine number of iterations. */
int ForLoop_integer::BeginFor(VariableArray const& CurrentVars) {
  if (!startVarName_.empty()) {
    std::string sval = CurrentVars.GetVariable( startVarName_ );
    if (sval.empty()) {
      mprinterr("Error: Start variable '%s' does not exist.\n", startVarName_.c_str());
      return NITERATIONS_UNKNOWN;
    }
    if (!validInteger(sval)) {
      mprinterr("Error: Variable '%s' does not contain a valid integer (%s)\n",
                startVarName_.c_str(), sval.c_str());
      return NITERATIONS_UNKNOWN;
    }
    start_ = convertToInteger(sval);
  }
  currentVal_ = start_;
  return calcNumIterations();
}

/** Check if integer for loop is done, increment if not. */
bool ForLoop_integer::EndFor(VariableArray& CurrentVars) {
  if (endOp_ == LESS_THAN) {
    if (currentVal_ >= end_) return true;
  } else if (endOp_ == GREATER_THAN) {
    if (currentVal_ <= end_) return true;
  }
  // Get variable value and update CurrentVars
  CurrentVars.UpdateVariable( VarName(), integerToString( currentVal_ ));
  // Increment
  currentVal_ += inc_;
  return false;
}
