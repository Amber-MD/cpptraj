#include "Exec_Flatten.h"
#include "CpptrajStdio.h"
#include "DataSet_2D.h"
#include "DataSet_1D.h"

// Exec_Flatten::Help()
void Exec_Flatten::Help() const
{

}

// Exec_Flatten::Execute()
Exec::RetType Exec_Flatten::Execute(CpptrajState& State, ArgList& argIn)
{
  std::vector<DataSet_2D*> inpSets;
  std::vector<DataSet_1D*> outSets;

  std::string dsname = argIn.GetStringKey("name");
  if (dsname.empty()) {
    mprinterr("Error: Must specify output set name with 'name'\n");
    return CpptrajState::ERR;
  }
  mprintf("\tOutput set name: %s\n", dsname.c_str());

  const char* modeStr[] = {"sum", "avg"};
  enum ModeType { SUM = 0, AVG };
  ModeType mode;
  std::string modeArg = argIn.GetStringKey("mode");
  if (!modeArg.empty()) {
    if (modeArg == "sum")
      mode = SUM;
    else if (modeArg == "avg")
      mode = AVG;
    else {
      mprinterr("Error: Unrecognized keyword for 'mode': %s\n", modeArg.c_str());
      return CpptrajState::ERR;
    }
  } else
    mode = SUM;
  mprintf("\tFlatten mode: %s\n", modeStr[mode]);

  // Get input sets
  std::string setarg = argIn.GetStringNext();
  while (!setarg.empty()) {
    DataSetList dsl = State.DSL().GetMultipleSets( setarg );
    for (DataSetList::const_iterator ds = dsl.begin(); ds != dsl.end(); ++ds)
    {
      if ( (*ds)->Group() != DataSet::MATRIX_2D ) {
        mprintf("Warning: Set '%s' is not a matrix, skipping.\n", (*ds)->legend());
      } else {
        inpSets.push_back( (DataSet_2D*)(*ds) );
      }
    }
    setarg = argIn.GetStringNext();
  }
  mprintf("\t%zu matrices to flatten.\n", inpSets.size());

  return CpptrajState::OK;  
}
