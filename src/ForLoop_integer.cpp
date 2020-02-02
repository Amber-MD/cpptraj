#include "ForLoop_integer.h"
#include "CpptrajStdio.h"
#include "ArgList.h"
#include "StringRoutines.h"
#include "CpptrajState.h"
#include "DataSetList.h"

const char* ForLoop_integer::OpStr_[NO_OP] = {"+=", "-=", "<", ">", "<=", ">=" };

ForLoop_integer::ForLoop_integer() :
  endOp_(NO_OP),
  incOp_(NO_OP),
  start_(0),
  end_(0),
  inc_(0),
  currentVal_(0),
  endArgPresent_(false)
{}

void ForLoop_integer::helpText() {
  mprintf("\t<var>=<start>;[<var><end OP><end>;]<var><increment OP>[<value>]\n"
          "    Available 'end OP'       : '<' '>' '<=' '>='\n"
          "    Available 'increment OP' : '++', '--', '+=', '-='\n"
          "  Loop over an integer value.\n");
}

/** Setup integer for loop. */
int ForLoop_integer::SetupFor(CpptrajState& State, ArgList& argIn) {
  // [<var>=<start>;[<var><OP><end>;]<var><OP>[<value>]]
  //SetType( INTEGER );
  // Sanity check
  if (argIn.Nargs() != 1) {
    mprinterr("Internal Error: Too many arguments for integer for loop.\n");
    return 1;
  }
  ArgList varArg( argIn.GetStringNext(), ";" );
  //varArg.PrintDebug(); // DEBUG
  if (varArg.Nargs() < 2 || varArg.Nargs() > 3) {
    mprinterr("Error: Malformed 'for' loop variable.\n"
              "Error: Expected '[<var>=<start>;[<var><OP><end>;]<var><OP>[<value>]]'\n"
              "Error: Got '%s'\n", argIn[0].c_str());
    return 1;
  }
  // First argument: <var>=<start>
  ArgList startArg( varArg[0], "=" );
  if (startArg.Nargs() != 2) {
    mprinterr("Error: Malformed 'start' argument.\n"
              "Error: Expected <var>=<start>, got '%s'\n", varArg[0].c_str());
    return 1;
  }
  if (SetupLoopVar( State.DSL(),startArg[0] )) return 1;
  //mprintf("DEBUG: Start argument: '%s' = '%s'\n", VarName().c_str(), startArg[1].c_str());
  if ( startArg[1][0] == '$' ) {
    // Variable name 
    startVarName_ = RemoveLeadingChars(startArg[1], 1);
  } else if (!validInteger(startArg[1])) {
    // Not Integer or variable
    mprinterr("Error: Start argument must be an integer or variable name.\n");
    return 1;
  } else
    start_ = convertToInteger(startArg[1]);
  // Second argument: <var><OP><end>
  size_t pos0 = VarName().size();
  size_t pos1 = pos0 + 1;
  endOp_ = NO_OP;
  int iargIdx = 1;
  if (varArg.Nargs() == 3) {
    iargIdx = 2;
    if ( varArg[1][pos0] == '<' ) {
      if ( varArg[1][pos0+1] == '=' ) {
        endOp_ = LT_EQUALS;
        pos1 = pos0 + 2;
      } else
        endOp_ = LESS_THAN;
    } else if ( varArg[1][pos0] == '>' ) {
      if ( varArg[1][pos0+1] == '=' ) {
        endOp_ = GT_EQUALS;
        pos1 = pos0 + 2;
      } else
        endOp_ = GREATER_THAN;
    }
    if (endOp_ == NO_OP) {
      mprinterr("Error: Unrecognized end op: '%s'\n",
                varArg[1].substr(pos0, pos1-pos0).c_str());
      return 1;
    }
    std::string endStr = varArg[1].substr(pos1);
    //mprintf("DEBUG: endStr= '%s'\n", endStr.c_str());
    if ( endStr[0] == '$' ) {
      // Variable name
      endVarName_ = RemoveLeadingChars(endStr, 1);
    } else if (!validInteger(endStr)) {
      // Not Integer of variable
      mprinterr("Error: End argument must be an integer or variable name.\n");
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
  if (iargIdx == 2) {
    // End argument present
    endArgPresent_ = true;
    if (endVarName_.empty())
      eval = integerToString(end_);
    else
      eval = endVarName_;
    description.append(VarName() + std::string(OpStr_[endOp_]) + eval + "; ");
  } else
    endArgPresent_ = false;
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
  if (!endArgPresent_) return NITERATIONS_UNKNOWN;
  // Check end > start for increment, start > end for decrement
  int maxval, minval;
  if (incOp_ == INCREMENT) {
    // When incrementing, start must be less than end, so only valid end
    // ops are < and <=.
    if (endOp_ != LESS_THAN && endOp_ != LT_EQUALS) {
      mprinterr("Error: For increment, only valid end ops are < and <=.\n");
      return LOOP_ERROR;
    }
    int offset = 0;
    if (endOp_ == LT_EQUALS) offset = inc_;
    if ( start_ >= end_ + offset )
    {
      mprintf("Warning: start (%i) and end (%i) values result in no iterations.\n",
                start_, end_ + offset);
      return 0;
    }
    minval = start_;
    maxval = end_ + offset;
  } else {
    // When decrementing, end must be less than start, so only valid end
    // ops are > and >=.
    if (endOp_ != GREATER_THAN && endOp_ != GT_EQUALS) {
      mprinterr("Error: For decrement, only valid end ops are > and >=.\n");
      return LOOP_ERROR;
    }
    int offset = 0;
    if (endOp_ == GT_EQUALS) offset = inc_;
    if ( (endOp_ == GREATER_THAN && end_ >= start_) ||
         (endOp_ == GT_EQUALS    && end_ > start_) )
    {
      mprintf("Warning: end (%i) and start (%i) values result in no iterations.\n",
                end_, start_);
      return 0;
    }
    minval = end_ + offset;
    maxval = start_;
  }
  // Figure out number of iterations
  int Niterations = (maxval - minval) / inc_;
  // Sanity check
  if (Niterations < 0) return LOOP_ERROR;
  if (((maxval-minval) % inc_) > 0) Niterations++;
  return Niterations;
}

/** Set start value and determine number of iterations. */
int ForLoop_integer::BeginFor(DataSetList const& DSL) {
  if (!startVarName_.empty()) {
    std::string sval = DSL.GetVariable( startVarName_ );
    if (sval.empty()) {
      mprinterr("Error: Start variable '%s' does not exist.\n", startVarName_.c_str());
      return LOOP_ERROR;
    }
    if (!validInteger(sval)) {
      mprinterr("Error: Variable '%s' does not contain a valid integer (%s)\n",
                startVarName_.c_str(), sval.c_str());
      return LOOP_ERROR;
    }
    start_ = convertToInteger(sval);
  }
  currentVal_ = start_;
  // Determine end arg if necessary
  if (endArgPresent_ && !endVarName_.empty()) {
    std::string eval = DSL.GetVariable( endVarName_ );
    if (eval.empty()) {
      mprinterr("Error: End variable '%s' does not exist.\n", endVarName_.c_str());
      return LOOP_ERROR;
    }
    if (!validInteger(eval)) {
      mprinterr("Error: Variable '%s' does not contain a valid integer (%s)\n",
                endVarName_.c_str(), eval.c_str());
      return LOOP_ERROR;
    }
    end_ = convertToInteger(eval);
  }
  return calcNumIterations();
}

/** Check if integer for loop is done, increment if not. */
bool ForLoop_integer::EndFor(DataSetList& DSL) {
  bool retval = false;
  if (endOp_ == LESS_THAN) {
    if (currentVal_ >= end_) retval = true;
  } else if (endOp_ == GREATER_THAN) {
    if (currentVal_ <= end_) retval = true;
  } else if (endOp_ == LT_EQUALS) {
    if ( currentVal_ > end_ ) retval = true;
  } else if (endOp_ == GT_EQUALS) {
    if ( currentVal_ < end_ ) retval = true;
  }
  // Get variable value and update CurrentVars
  DSL.UpdateStringVar( VarName(), integerToString( currentVal_ ) );
  // Increment
  currentVal_ += inc_;
  return retval;
}
