#include "ControlBlock_For.h"
#include "CpptrajStdio.h"
#include "ArgList.h"
#include "DataSetList.h" // TODO move into ForLoop help?
// For loop types
#include "ForLoop.h"
#include "ForLoop_integer.h"
#include "ForLoop_mask.h"
#include "ForLoop_list.h"

/// DESTRUCTOR
ControlBlock_For::~ControlBlock_For() {
  for (Marray::iterator it = Vars_.begin(); it != Vars_.end(); ++it)
    delete *it;
  Vars_.clear();
}

void ControlBlock_For::Help() const {
  mprintf("\t{ {atoms|residues|molecules|molfirstres|mollastres}\n"
          "\t    <var> inmask <mask> [%s] ... |\n"
          "\t    <var> in <list>\n"
          "\t  <var>=<start>;[<var><end OP><end>;]<var><increment OP>[<value>] ... }\n",
          DataSetList::TopIdxArgs);
  mprintf("\tEND KEYWORD: 'done'\n");
  mprintf("  Create a 'for' loop around specified mask(s), comma-separated list of\n"
          "  strings and/or integer value(s). Any number and combination of masks,\n"
          "  lists, and integers can be specified. Strings in lists can be file\n"
          "  names containing wildcard characters ('*' or '?').\n"
          "  Variables created in the for loop can be referenced by prefacing\n"
          "  the name with a '$' character.\n"
          "  Available 'end OP'       : '<' '>'\n"
          "  Available 'increment OP' : '++', '--', '+=', '-='\n"
          "  Note that variables are NOT incremented after the final loop iteration,\n"
          "  i.e. loop variables always retain their final value.\n"
          "  Examples:\n"
          "\tfor atoms A0 inmask :1-27&!@H= i=1;i++\n"
          "\t  distance d$i :TCS $A0 out $i.dat\n"
          "\tdone\n"
          "\tfor TRAJ in trajA*.nc,trajB*.nc\n"
          "\t  trajin $TRAJ 1 last 10\n"
          "\tdone\n"
          );
}

/** Set up each mask/integer loop. */
int ControlBlock_For::SetupBlock(CpptrajState& State, ArgList& argIn) {
  mprintf("    Setting up 'for' loop.\n");
  Vars_.clear();
  description_.assign("for (");
  int MaxIterations = -1;
  int iarg = 0;
  while (iarg < argIn.Nargs())
  {
    // Advance to next unmarked argument.
    while (iarg < argIn.Nargs() && argIn.Marked(iarg)) iarg++;
    if (iarg == argIn.Nargs()) break;
    // Determine 'for' type
    int argToMark = iarg;
    if      ( argIn[iarg] == "atoms"       || 
              argIn[iarg] == "residues"    || 
              argIn[iarg] == "molecules"   || 
              argIn[iarg] == "molfirstres" || 
              argIn[iarg] == "mollastres" )
    {
      Vars_.push_back( static_cast<ForLoop*>( new ForLoop_mask() ) );
    } else if ( argIn[iarg].find(";") != std::string::npos ) {
      Vars_.push_back( static_cast<ForLoop*>( new ForLoop_integer() ) );
    } else if (iarg+1 < argIn.Nargs() && argIn[iarg+1] == "in") {
      Vars_.push_back( static_cast<ForLoop*>( new ForLoop_list() ) );
    } else {
      // Exit if type could not be determined.
      mprinterr("Error: for loop type not specfied.\n");
      return 1;
    }
    argIn.MarkArg(argToMark);
    ForLoop& forloop = static_cast<ForLoop&>( *(Vars_.back()) );
    if ( forloop.SetupFor( State, argIn[iarg], argIn ) ) {
      mprinterr("Error: For loop setup failed.\n");
      return 1;
    }
    // Check number of values
    int Niterations = forloop.Niterations();
    if (MaxIterations == -1)
      MaxIterations = Niterations;
    else if (Niterations != -1 && Niterations != MaxIterations) {
      mprintf("Warning: # iterations %i != previous # iterations %i\n",
              Niterations, MaxIterations);
      MaxIterations = std::min(Niterations, MaxIterations);
    }
    // Append description
    description_.append( forloop.Description() );
    // TODO check variable name
  }
  mprintf("\tLoop will execute for %i iterations.\n", MaxIterations);
  if (MaxIterations < 1) {
    mprinterr("Error: Loop has less than 1 iteration.\n");
    return 1; 
  }
  description_.append(") do");

  return 0;
}

/** \return true if the given command will finish this block. */
bool ControlBlock_For::EndBlock(ArgList const& a) const {
  return (a.CommandIs("done"));
}

/** For each mask initialize iterator. For each integer set to start value. */
int ControlBlock_For::Start(VariableArray const& CurrentVars) {
  // TODO search currentvars
  for (Marray::iterator MH = Vars_.begin(); MH != Vars_.end(); ++MH)
  {
    if ((*MH)->BeginFor(CurrentVars)) return 1;
  }
  return 0;
}

/** For each mask check if done, then update CurrentVars, then increment. */
ControlBlock::DoneType ControlBlock_For::CheckDone(VariableArray& CurrentVars) {
  for (Marray::iterator MH = Vars_.begin(); MH != Vars_.end(); ++MH) {
    // Exit as soon as one is done TODO check all?
    if ( (*MH)->EndFor(CurrentVars) ) return DONE;
  }

  return NOT_DONE;
}
