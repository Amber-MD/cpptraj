#include "ControlBlock_For.h"
#include "CpptrajStdio.h"
#include "ArgList.h"
#include "CpptrajState.h" // TODO move into ForLoop help?
// For loop types
#include "ForLoop.h"
#include "ForLoop_integer.h"
#include "ForLoop_mask.h"
#include "ForLoop_list.h"
#include "ForLoop_dataSetBlocks.h"

/// DESTRUCTOR
ControlBlock_For::~ControlBlock_For() {
  for (Marray::iterator it = Vars_.begin(); it != Vars_.end(); ++it)
    delete *it;
  Vars_.clear();
}

void ControlBlock_For::Help() const {
  mprintf("\t{ {atoms|residues|molecules|molfirstres|mollastres}\n"
          "\t    <var> inmask <mask> [%s] ... |\n"
          "\t  <var> in <list> |\n"
          "\t  <var>=<start>;[<var><end OP><end>;]<var><increment OP>[<value>] ... }\n",
          DataSetList::TopIdxArgs);
  mprintf("\tEND KEYWORD: 'done'\n");
  mprintf("  Create a 'for' loop around specified mask(s), comma-separated list of\n"
          "  strings and/or integer value(s). Any number and combination of masks,\n"
          "  lists, and integers can be specified. Strings in lists can be file\n"
          "  names containing wildcard characters ('*' or '?').\n"
          "  Variables created in the for loop can be referenced by prefacing\n"
          "  the name with a '$' character.\n"
          "  Available 'end OP'       : '<' '>' '<=' '>='\n"
          "  Available 'increment OP' : '++', '--', '+=', '-='\n"
          "  Note that non-integer variables (e.g. for mask loops) are NOT incremented\n"
          "  after the final loop iteration, i.e. these loop variables always retain\n"
          "  their final value.\n"
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
  // Parse out each component of this for block. E.g.
  //               .               .                 .          .
  //   for atoms X inmask :30 name in var1,var2,var3 i=0;i++ DS datasetblocks blocksize 10
  //   0   1     2 3      4   5    6  7              8       9  10            11        12
  //       *                  *                      *       *
  // component 1: 1,2,3,4
  // component 2: 5,6,7
  // component 3: 8
  argIn.MarkArg(0);
  std::vector<int> forLoopIdxs;
  for (int iarg = 1; iarg < argIn.Nargs(); iarg++) {
    if (argIn[iarg] == "inmask") {
      int idx = iarg - 2;
      if (idx < 1) {
        mprinterr("Error: Malformed 'inmask' for loop.\n");
        return 1;
      }
      forLoopIdxs.push_back( idx );
      Vars_.push_back( static_cast<ForLoop*>( new ForLoop_mask() ) );
    } else if (argIn[iarg] == "in") {
      int idx = iarg - 1;
      if (idx < 1) {
        mprinterr("Error: Malformed 'in' for loop.\n");
        return 1;
      }
      forLoopIdxs.push_back( idx );
      Vars_.push_back( static_cast<ForLoop*>( new ForLoop_list() ) );
    } else if (argIn[iarg] == "datasetblocks") {
      int idx = iarg - 1;
      if (idx < 1) {
        mprinterr("Error: Malformed 'datasetblocks' for loop.\n");
        return 1;
      }
      forLoopIdxs.push_back( idx );
      Vars_.push_back( static_cast<ForLoop*>( new ForLoop_dataSetBlocks() ) );
    } else if ( argIn[iarg].find(";") != std::string::npos ) {
      forLoopIdxs.push_back( iarg );
      Vars_.push_back( static_cast<ForLoop*>( new ForLoop_integer() ) );
    }
    argIn.MarkArg(iarg);
  }
  forLoopIdxs.push_back( argIn.Nargs() );
  mprintf("DEBUG: For loop indices:");
  for (std::vector<int>::const_iterator it = forLoopIdxs.begin();
                                        it != forLoopIdxs.end(); ++it)
    mprintf(" %i", *it);
  mprintf("\n");
  if (forLoopIdxs.size() < 2) {
    mprinterr("Error: No recognized for loop arguments.\n");
    return 1;
  }
  // Set up each individual for loop.
  for (unsigned int idx = 1; idx < forLoopIdxs.size(); ++idx)
  {
    ArgList forLoopArgs;
    for (int jdx = forLoopIdxs[idx-1]; jdx < forLoopIdxs[idx]; ++jdx)
      forLoopArgs.AddArg( argIn[jdx] );
    forLoopArgs.PrintDebug();
    ForLoop& forloop = static_cast<ForLoop&>( *(Vars_[idx-1]) );
    if ( forloop.SetupFor( State, forLoopArgs ) ) {
      mprinterr("Error: For loop setup failed.\n");
      return 1;
    }
    if ( !forloop.IsSetup() ) {
      mprinterr("Internal Error: For loop variable was not properly set up.\n");
      return 1;
    }
    if (forLoopArgs.CheckForMoreArgs()) {
      if (State.ExitOnError()) {
        mprinterr("Error: Not all for loop arguments handled.\n");
        return 1;
      } else
        mprintf("Warning: Not all for loop arguments handled.\n");
    }
    // Append description
    description_.append( forloop.Description() );
    // TODO check variable name
  }

  description_.append(") do");

  return 0;
}

/** \return true if the given command will finish this block. */
bool ControlBlock_For::EndBlock(ArgList const& a) const {
  return (a.CommandIs("done"));
}

/** Set all loops to their initial values and determine # iterations
  * \return 0 if all loops were set up, 1 if an error occurred.
  */
int ControlBlock_For::Start(DataSetList const& DSL) {
  int MaxIterations = -1;
  for (Marray::iterator MH = Vars_.begin(); MH != Vars_.end(); ++MH)
  {
    int Niterations = ((*MH)->BeginFor(DSL));
    if (Niterations == ForLoop::LOOP_ERROR) return 1;
    // Check number of values
    if (Niterations != ForLoop::NITERATIONS_UNKNOWN)
    {
      if (MaxIterations == -1)
        MaxIterations = Niterations;
      else {
        if (Niterations != MaxIterations)
          mprintf("Warning: # iterations %i != previous # iterations %i\n",
                  Niterations, MaxIterations);
        MaxIterations = std::min(Niterations, MaxIterations);
      }
      mprintf("\tLoop over '%s' will execute for %i iterations.\n",
              (*MH)->VarName().c_str(), Niterations);
    } else
      mprintf("\tLoop over '%s' has unknown # iterations.\n", (*MH)->VarName().c_str());
  }
  if (Vars_.size() > 1)
    mprintf("\tLoop will execute for %i iterations.\n", MaxIterations);
  if (MaxIterations < 1) {
    mprintf("Warning: Loop has less than 1 iteration.\n");
    return 0; 
  }

  return 0;
}

/** For each mask check if done, then update variables, then increment. */
ControlBlock::DoneType ControlBlock_For::CheckDone(DataSetList& DSL) {
  // NOTE: DataSetList is not const ref so dataSetBlocks loop can create sets
  DoneType retval = NOT_DONE;
  for (Marray::iterator MH = Vars_.begin(); MH != Vars_.end(); ++MH) {
    // Exit as soon as one is done, but check all so that
    // all loops are properly incremented if necessary.
    if ( (*MH)->EndFor(DSL) ) retval = DONE;
  }

  return retval;
}
