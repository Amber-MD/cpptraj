#include "Analysis_HausdorffDistance.h"
#include "CpptrajStdio.h"

// Analysis_HausdorffDistance::Help()
void Analysis_HausdorffDistance::Help() const {

}

/** Assume input matrix contains all distances beween sets A and B.
  * The Hausdorff distance will be the maximum of all minimums in each row.
  */
double Analysis_HausdorffDistance::h_Matrix(DataSet_2D* m1) {
  if (m1 == 0) {
    mprinterr("Internal Error: Analysis_HausdorffDistance::h_Matrix(): Matrix set is null.\n");
    return -1.0;
  }
  if (m1->Size() < 1) {
    mprinterr("Error: '%s' is empty.\n", m1->legend());
    return -1.0;
  }
  double hd = 0.0;
  // Row 1 - initial value.
  for (unsigned int row = 0; row != m1->Nrows(); row++)
  {
    double minRow = m1->GetElement(0, row);
    for (unsigned int col = 1; col != m1->Ncols(); col++)
      minRow = std::min( minRow, m1->GetElement(col, row) );
    hd = std::max( hd, minRow );
  }
  return hd;
}

// Analysis_HausdorffDistance::Setup()
Analysis::RetType Analysis_HausdorffDistance::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  std::string dsarg = analyzeArgs.GetStringNext();
  while (!dsarg.empty()) {
    inputSets_ += setup.DSL().GetMultipleSets( dsarg );
    dsarg = analyzeArgs.GetStringNext();
  }
  if (inputSets_.empty()) {
    mprinterr("Error: No data sets specified.\n");
    return Analysis::ERR;
  }
  mprintf("    HAUSDORFF:\n");
  for (DataSetList::const_iterator it = inputSets_.begin(); it != inputSets_.end(); ++it)
    mprintf("\t%s\n", (*it)->legend());

  return Analysis::OK;
}

// Analysis_HausdorffDistance::Analyze()
Analysis::RetType Analysis_HausdorffDistance::Analyze() {
  // Get the Hausdorff distance for each set.
  double hd = 0.0;
  for (DataSetList::const_iterator it = inputSets_.begin(); it != inputSets_.end(); ++it)
  {
    hd = -1.0;
    if ( (*it)->Group() == DataSet::MATRIX_2D )
      hd = h_Matrix( static_cast<DataSet_2D*>( *it ) );
    else
      mprintf("Warning: '%s' type not yet supported for Hausdorff\n", (*it)->legend());
    mprintf("%12.4f %s\n", hd, (*it)->legend());
  }
  return Analysis::OK;
}
