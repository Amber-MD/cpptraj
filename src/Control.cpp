#include <algorithm> // std::min
#include "Control.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"

void Control_For::Help() const {
  mprintf("\t{ {atoms|residues|molecules|molfirstres|mollastres}\n"
          "\t    <var> inmask <mask> %s ... |\n"
          "\t  <var>=<start>;[<var><end OP><end>;]<var><increment OP>[<value>] ... }\n",
          DataSetList::TopIdxArgs);
  mprintf("\tEND KEYWORD: 'done'\n");
  mprintf("  Create a 'for' loop around specified mask(s) and/or integer value(s).\n"
          "  Any number and combination of masks and integers can be specified.\n"
          "  Variables created in the for loop can be referenced by prefacing\n"
          "  the name with a '$' character.\n"
          "  Available 'end OP'       : '<' '>'\n"
          "  Available 'increment OP' : '++', '--', '+=', '-='\n"
          "  Note that variables are NOT incremented after the final loop iteration,\n"
          "  i.e. loop variables always retain their final value.\n"
          "  Example:\n"
          "\tfor atoms A0 inmask :1-27&!@H= i=1;i++\n"
          "\t  distance d$i :TCS $A0 out $i.dat\n"
          "\tdone\n");
}

/** Set up each mask. */
int Control_For::SetupControl(CpptrajState& State, ArgList& argIn, Varray& CurrentVars) {
  mprintf("    Setting up 'for' loop.\n");
  Vars_.clear();
  Topology* currentTop = 0;
  static const char* TypeStr[] = { "ATOMS ", "RESIDUES ", "MOLECULES ",
                                   "MOL_FIRST_RES ", "MOL_LAST_RES " };
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
    bool isMaskFor = true;
    if      ( argIn[iarg] == "atoms"     ) ftype = ATOMS;
    else if ( argIn[iarg] == "residues"  ) ftype = RESIDUES;
    else if ( argIn[iarg] == "molecules" ) ftype = MOLECULES;
    else if ( argIn[iarg] == "molfirstres" ) ftype = MOLFIRSTRES;
    else if ( argIn[iarg] == "mollastres"  ) ftype = MOLLASTRES;
    else if ( argIn[iarg].find(";") != std::string::npos ) {
      isMaskFor = false;
      ftype = INTEGER;
    }
    if (ftype == UNKNOWN) {
      mprinterr("Error: for loop type not specfied.\n");
      return 1;
    }
    argIn.MarkArg(iarg);
    Vars_.push_back( LoopVar() );
    LoopVar& MH = Vars_.back();
    int Niterations = -1;
    // Set up for specific type
    if (description_ != "for (") description_.append(", ");
    // -------------------------------------------
    if (isMaskFor)
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
      } else if (MH.varType_ == MOLECULES ||
                 MH.varType_ == MOLFIRSTRES ||
                 MH.varType_ == MOLLASTRES)
      {
        int curMol = -1;
        for (AtomMask::const_iterator at = currentMask.begin(); at != currentMask.end(); ++at) {
          int mol = (*currentTop)[*at].MolNum();
          if (mol != curMol) {
            if (MH.varType_ == MOLECULES)
              MH.Idxs_.push_back( mol );
            else {
              int res;
              if (MH.varType_ == MOLFIRSTRES)
                res = (*currentTop)[ currentTop->Mol( mol ).BeginAtom() ].ResNum();
              else // MOLLASTRES
                res = (*currentTop)[ currentTop->Mol( mol ).EndAtom()-1 ].ResNum();
              MH.Idxs_.push_back( res );
            }
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
        // TODO allow variables
        mprinterr("Error: Start argument must be an integer.\n");
        return 1;
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
          // TODO allow variables
          mprinterr("Error: End argument must be an integer.\n");
          return 1;
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
      std::string sval = integerToString(MH.start_);
      description_.append("(" + MH.varname_ + "=" + sval + "; ");
      std::string eval;
      if (iargIdx == 2) {
        // End argument present
        eval = integerToString(MH.end_);
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

/** For each mask initialize iterator. For each integer set to start value. */
void Control_For::Start() {
  for (Marray::iterator MH = Vars_.begin(); MH != Vars_.end(); ++MH) {
    if (MH->varType_ == INTEGER)
      MH->currentVal_ = MH->start_; // TODO search currentvars
    else
      MH->idx_ = MH->Idxs_.begin();
  }
}

/** For each mask check if done, then update CurrentVars, then increment. */
Control::DoneType Control_For::CheckDone(Varray& CurrentVars) {
  static const char* prefix[] = {"@", ":", "^", ":", ":"};
  for (Marray::iterator MH = Vars_.begin(); MH != Vars_.end(); ++MH) {
    // Exit as soon as one is done TODO check all?
    if (MH->varType_ == INTEGER) {
      if (MH->endOp_ == LESS_THAN) {
        if (MH->currentVal_ >= MH->end_) return DONE;
      } else if (MH->endOp_ == GREATER_THAN) {
        if (MH->currentVal_ <= MH->end_) return DONE;
      }
      // Get variable value and update CurrentVars
      CurrentVars.UpdateVariable( MH->varname_, integerToString( MH->currentVal_ ));
      // Increment
      MH->currentVal_ += MH->inc_;
    } else {
      if (MH->idx_ == MH->Idxs_.end()) return DONE;
      // Get variable value
      std::string maskStr = prefix[MH->varType_] + integerToString(*(MH->idx_) + 1);
      //mprintf("DEBUG: Control_For: %s\n", maskStr.c_str());
      // Update CurrentVars
      CurrentVars.UpdateVariable( MH->varname_, maskStr );
      // Increment
      ++(MH->idx_);
    }
  }
  return NOT_DONE;
}

// =============================================================================
void Control_Set::Help() const {
  mprintf("\t<variable>=<value>\n");
}

int Control_Set::SetupControl(CpptrajState& State, ArgList& argIn, Varray& CurrentVars)
{
  int iarg = 0;
  while (iarg < argIn.Nargs()) {
    while (iarg < argIn.Nargs() && argIn.Marked(iarg)) iarg++;
    if (iarg == argIn.Nargs()) break;
    if ( argIn[iarg].find("=") != std::string::npos ) {
      // <var>=<value>
      ArgList equals( argIn[iarg], "=" );
      if (equals.Nargs() != 2) {
        mprinterr("Error: Expected <var>=<value>\n");
        return 1;
      }
      CurrentVars.UpdateVariable( "$" + equals[0], equals[1] );
      mprintf("\tVariable '%s' set to '%s'\n", equals[0].c_str(), equals[1].c_str());
      argIn.MarkArg(iarg);
    } else {
      mprinterr("Error: Unrecognized syntax: %s\n", argIn[iarg].c_str());
      return 1;
    }
  }
  return 0;
}
