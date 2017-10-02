#include <algorithm> // std::min
#include "Control.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"

void Control_For::Help() const {
  mprintf("\t[ {atoms|residues|molecules} <var> inmask <mask> %s ...\n"
          "\t  <var>=<start>;[<var><OP><end>;]<var><OP>[<value>] ... ]\n",
          DataSetList::TopIdxArgs);
}

/** Set up each mask. */
int Control_For::SetupControl(CpptrajState& State, ArgList& argIn) {
  mprintf("    Setting up 'for' loop.\n");
  Masks_.clear();
  Topology* currentTop = 0;
  static const char* TypeStr[] = { "ATOMS ", "RESIDUES ", "MOLECULES " };
  static const char* OpStr[] = {"+=", "-=", "<", ">"};
  description_.assign("for (");
  int MaxIterations = -1;
  int iarg = 0;
  while (iarg < argIn.Nargs())
  {
    // Advance to next unmarked argument.
    while (iarg < argIn.Nargs() && argIn.Marked(iarg)) iarg++;
    if (iarg == argIn.Nargs()) break;
    // Determine 'for' type
    ForType ftype = UNKNOWN;
    if      ( argIn[iarg] == "atoms"     ) ftype = ATOMS;
    else if ( argIn[iarg] == "residues"  ) ftype = RESIDUES;
    else if ( argIn[iarg] == "molecules" ) ftype = MOLECULES;
    else if ( argIn[iarg].find(";") != std::string::npos )
      ftype = INTEGER;
    if (ftype == UNKNOWN) {
      mprinterr("Error: One of {atoms|residues|molecules} not specfied.\n");
      return 1;
    }
    argIn.MarkArg(iarg);
    Masks_.push_back( MaskHolder() );
    MaskHolder& MH = Masks_.back();
    int Niterations = -1;
    // Set up for specific type
    if (description_ != "for (") description_.append(", ");
    // -------------------------------------------
    if (ftype == ATOMS || ftype == RESIDUES || ftype == MOLECULES)
    {
      // {atoms|residues|molecules} <var> inmask <mask> [TOP KEYWORDS]
      if (argIn[iarg+2] != "inmask") {
        mprinterr("Error: Expected 'inmask', got %s\n", argIn[iarg+2].c_str());
        return 1;
      }
      AtomMask currentMask;
      if (currentMask.SetMaskString( argIn.GetStringKey("inmask") )) return 1;
      MH.varType_ = ftype;
      Topology* top = State.DSL().GetTopByIndex( argIn );
      if (top != 0) currentTop = top;
      if (currentTop == 0) return 1;
      MH.varname_ = argIn.GetStringNext();
      if (MH.varname_.empty()) {
        mprinterr("Error: 'for inmask': missing variable name.\n");
        return 1;
      }
      MH.varname_ = "$" + MH.varname_;
      // Set up mask
      if (currentTop->SetupIntegerMask( currentMask )) return 1;
      currentMask.MaskInfo();
      if (currentMask.None()) return 1;
      // Set up indices
      if (MH.varType_ == ATOMS)
        MH.Idxs_ = currentMask.Selected();
      else if (MH.varType_ == RESIDUES) {
        int curRes = -1;
        for (AtomMask::const_iterator at = currentMask.begin(); at != currentMask.end(); ++at) {
          int res = (*currentTop)[*at].ResNum();
          if (res != curRes) {
            MH.Idxs_.push_back( res );
            curRes = res;
          }
        }
      } else if (MH.varType_ == MOLECULES) {
        int curMol = -1;
        for (AtomMask::const_iterator at = currentMask.begin(); at != currentMask.end(); ++at) {
          int mol = (*currentTop)[*at].MolNum();
          if (mol != curMol) {
            MH.Idxs_.push_back( mol );
            curMol = mol;
          }
        }
      }
      Niterations = (int)MH.Idxs_.size();
      description_.append(std::string(TypeStr[MH.varType_]) +
                        MH.varname_ + " inmask " + currentMask.MaskExpression());
    // -------------------------------------------
    } else if (ftype == INTEGER) {
      // [<var>=<start>;[<var><OP><end>;]<var><OP>[<value>]]
      MH.varType_ = ftype;
      ArgList varArg( argIn[iarg], ";" );
      if (varArg.Nargs() < 2 || varArg.Nargs() > 3) {
        mprinterr("Error: Malformed 'for' loop variable.\n"
                  "Error: Expected '[<var>=<start>;[<var><OP><end>;]<var><OP>[<value>]]'\n"
                  "Error: Got '%s'\n", argIn[iarg].c_str());
        return 1;
      }
      // First argument: <var>=<start>
      ArgList startArg( varArg[0], "=" );
      if (startArg.Nargs() != 2) {
        mprinterr("Error: Malformed 'start' argument.\n"
                  "Error: Expected <var>=<start>, got '%s'\n", varArg[0].c_str());
        return 1;
      }
      MH.varname_ = startArg[0];
      if (!validInteger(startArg[1])) {
        if (startArg[1][0] != '$') {
          mprinterr("Error: Expected start argument to be integer or variable.\n");
          return 1;
        }
        MH.startArg_ = startArg[1];
      } else
        MH.start_ = convertToInteger(startArg[1]);
      // Second argument: <var><OP><end>
      size_t pos0 = MH.varname_.size();
      size_t pos1 = pos0 + 1;
      MH.endOp_ = NO_OP;
      int iargIdx = 1;
      if (varArg.Nargs() == 3) {
        iargIdx = 2;
        if ( varArg[1][pos0] == '<' )
          MH.endOp_ = LESS_THAN;
        else if (varArg[1][pos0] == '>')
          MH.endOp_ = GREATER_THAN;
        if (MH.endOp_ == NO_OP) {
          mprinterr("Error: Unrecognized end op: '%s'\n",
                    varArg[1].substr(pos0, pos1-pos0).c_str());
          return 1;
        }
        std::string endStr = varArg[1].substr(pos1);
        if (!validInteger(endStr)) {
          if (endStr[0] != '$') {
            mprinterr("Error: Expected end argument to be integer or variable.\n");
            return 1;
          }
          MH.endArg_ = endStr;
        } else
          MH.end_ = convertToInteger(endStr);
      }
      // Third argument: <var><OP>[<value>]
      pos1 = pos0 + 2;
      MH.incOp_ = NO_OP;
      bool needValue = false;
      if ( varArg[iargIdx][pos0] == '+' ) {
        if (varArg[iargIdx][pos0+1] == '+') {
          MH.incOp_ = INCREMENT;
          MH.inc_ = 1;
        } else if (varArg[iargIdx][pos0+1] == '=') {
          MH.incOp_ = INCREMENT;
          needValue = true;
        }
      } else if ( varArg[iargIdx][pos0] == '-' ) {
        if (varArg[iargIdx][pos0+1] == '-' ) {
          MH.incOp_ = DECREMENT;
          MH.inc_ = 1;
        } else if (varArg[iargIdx][pos0+1] == '=') {
          MH.incOp_ = DECREMENT;
          needValue = true;
        }
      }
      if (MH.incOp_ == NO_OP) {
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
        MH.inc_ = convertToInteger(incStr);
        if (MH.inc_ < 1) {
          mprinterr("Error: Extra '-' detected in increment.\n");
          return 1;
        }
      }
      // Description
      MH.varname_ = "$" + MH.varname_;
      std::string sval, eval;
      if (MH.startArg_.empty())
        sval = integerToString(MH.start_);
      else
        sval = MH.startArg_;
      description_.append("(" + MH.varname_ + "=" + sval + "; ");
      if (iargIdx == 2) {
        // End argument present
        if (MH.endArg_.empty())
          eval = integerToString(MH.end_);
        else
          eval = MH.endArg_;
        description_.append(MH.varname_ + std::string(OpStr[MH.endOp_]) + eval + "; ");
        // Check end > start for increment, start > end for decrement
        int maxval, minval;
        if (MH.incOp_ == INCREMENT) {
          if (MH.start_ >= MH.end_) {
            mprinterr("Error: start must be less than end for increment.\n");
            return 1;
          }
          minval = MH.start_;
          maxval = MH.end_;
        } else {
          if (MH.end_ >= MH.start_) {
            mprinterr("Error: end must be less than start for decrement.\n");
            return 1;
          }
          minval = MH.end_;
          maxval = MH.start_;
        }
        // Figure out number of iterations
        Niterations = (maxval - minval) / MH.inc_;
        if (((maxval-minval) % MH.inc_) > 0) Niterations++;
      }
      description_.append( MH.varname_ + std::string(OpStr[MH.incOp_]) +
                           integerToString(MH.inc_) + ")" );
      // If decrementing just negate value
      if (MH.incOp_ == DECREMENT)
        MH.inc_ = -MH.inc_;
      // DEBUG
      //mprintf("DEBUG: start=%i endOp=%i end=%i incOp=%i val=%i startArg=%s endArg=%s\n",
      //        MH.start_, (int)MH.endOp_, MH.end_, (int)MH.incOp_, MH.inc_,
      //        MH.startArg_.c_str(), MH.endArg_.c_str());
    }
    // Check number of values
    if (MaxIterations == -1)
      MaxIterations = Niterations;
    else if (Niterations != -1 && Niterations != MaxIterations) {
      mprintf("Warning: # iterations %i != previous # iterations %i\n",
              Niterations, MaxIterations);
      MaxIterations = std::min(Niterations, MaxIterations);
    }
  }
  mprintf("\tLoop will execute for %i iterations.\n", MaxIterations);
  if (MaxIterations < 1) {
    mprinterr("Error: Loop has less than 1 iteration.\n");
    return 1; 
  }
  description_.append(") do");

  return 0;
}

/** For each mask add variable to CurrentVars and initialize iterator. */
void Control_For::Start(Varray& CurrentVars) {
  for (Marray::iterator MH = Masks_.begin(); MH != Masks_.end(); ++MH) {
    if (MH->varType_ == INTEGER)
      MH->currentVal_ = MH->start_; // TODO search currentvars
    else
      MH->idx_ = MH->Idxs_.begin();
    // Init CurrentVars
    CurrentVars.push_back( VarPair(MH->varname_, "") );
  }
}

static inline void UpdateCurrentVars(Control::Varray& CurrentVars, std::string const& varname,
                                     std::string const& value)
{
  //mprintf("DEBUG: UpdateCurrentVars: %s = %s\n", varname.c_str(), value.c_str());
  Control::Varray::iterator it = CurrentVars.begin();
  for (; it != CurrentVars.end(); ++it) {
    if (it->first == varname) {
      it->second = value;
      break;
    }
  }
}

/** For each mask check if done, then update CurrentVars, then increment. */
Control::DoneType Control_For::CheckDone(Varray& CurrentVars) {
  static const char* prefix[] = {"@", ":", "^"};
  for (Marray::iterator MH = Masks_.begin(); MH != Masks_.end(); ++MH) {
    // Exit as soon as one is done TODO check all?
    if (MH->varType_ == INTEGER) {
      if (MH->endOp_ == LESS_THAN) {
        if (MH->currentVal_ >= MH->end_) return DONE;
      } else if (MH->endOp_ == GREATER_THAN) {
        if (MH->currentVal_ <= MH->end_) return DONE;
      }
      // Get variable value and update CurrentVars
      UpdateCurrentVars(CurrentVars, MH->varname_, integerToString( MH->currentVal_ ));
      // Increment
      MH->currentVal_ += MH->inc_;
    } else {
      if (MH->idx_ == MH->Idxs_.end()) return DONE;
      // Get variable value
      std::string maskStr = prefix[MH->varType_] + integerToString(*(MH->idx_) + 1);
      //mprintf("DEBUG: Control_For: %s\n", maskStr.c_str());
      // Update CurrentVars
      UpdateCurrentVars(CurrentVars, MH->varname_, maskStr);
      // Increment
      ++(MH->idx_);
    }
  }
  return NOT_DONE;
}
