#include "Exec_ClusterMap.h"

void Exec_ClusterMap::Help() const {
  mprintf("\t<2d set>\n");
}


Exec::RetType Exec_ClusterMap::Execute(CpptrajState& State, ArgList& argIn)
{
  DataSet* dsIn = State.DSL().GetDataSet( argIn.GetStringNext() );
  if (dsIn == 0) return CpptrajState::ERR;
  mprintf("\tSet '%s'\n", dsIn->legend());
  if (dsIn->Group() != DataSet::MATRIX_2D) {
    mprinterr("Error: Set is not 2D.\n");
    return CpptrajState::ERR;
  }
  if (dsIn->Size() < 1) {
    mprinterr("Error: Set is empty.\n");
    return CpptrajState::ERR;
  }

  DataSet_2D const& matrix = static_cast<DataSet_2D const&>( *dsIn );
  // Go through the set, calculate the average. Also determine the max.
  double maxVal = matrix.GetElement(0);
  double avg = 0.0;
  for (unsigned int i = 1; i != matrix.Size(); i++) {
    double val = matrix.GetElement(i);
    maxVal = std::max(val, maxVal);
    avg += val;
  }
  avg /= (double)matrix.Size();
  mprintf("\t%zu elements, max= %f, avg= %f\n", matrix.Size(), maxVal, avg);

  return CpptrajState::OK;
}
