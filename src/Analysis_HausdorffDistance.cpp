#include "Analysis_HausdorffDistance.h"
#include "CpptrajStdio.h"
#include "DataSet_MatrixFlt.h"

Analysis_HausdorffDistance::Analysis_HausdorffDistance() :
  outType_(BASIC),
  out_(0)
{}

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
    //mprintf("DEBUG: Min row %6u is %12.4f\n", row, minRow);
    hd = std::max( hd, minRow );
  }
  return hd;
}

// Analysis_HausdorffDistance::Help()
void Analysis_HausdorffDistance::Help() const {
  mprintf("\t<set arg0> [<set arg1> ...] [outtype {basic|trimatrix nrows <#>}]\n"
          "\t[name <output set name>] [out <filename>]\n");
}


// Analysis_HausdorffDistance::Setup()
Analysis::RetType Analysis_HausdorffDistance::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  // Keywords
  int nrows = -1;
  std::string outtypearg = analyzeArgs.GetStringKey("outtype");
  if (!outtypearg.empty()) {
    if (outtypearg == "basic")
      outType_ = BASIC;
    else if (outtypearg == "trimatrix") {
      outType_ = UPPER_TRI_MATRIX;
      nrows = analyzeArgs.getKeyInt("nrows", -1);
      if (nrows < 1) {
        mprinterr("Error: 'nrows' must be specified and > 0 for 'trimatrix'\n");
        return Analysis::ERR;
      }
    } else {
      mprinterr("Error: Unrecognized keyword for 'outtype': %s\n", outtypearg.c_str());
      return Analysis::ERR;
    }
  } else
    outType_ = BASIC;
  std::string dsname = analyzeArgs.GetStringKey("name");
  DataFile* df = setup.DFL().AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs );
  // Get input data sets
  std::string dsarg = analyzeArgs.GetStringNext();
  while (!dsarg.empty()) {
    inputSets_ += setup.DSL().GetMultipleSets( dsarg );
    dsarg = analyzeArgs.GetStringNext();
  }
  if (inputSets_.empty()) {
    mprinterr("Error: No data sets specified.\n");
    return Analysis::ERR;
  }
  // Output data set
  out_ = 0;
  if (outType_ == BASIC) {
    out_ = setup.DSL().AddSet(DataSet::FLOAT, dsname, "HAUSDORFF");
    if (out_ == 0) return Analysis::ERR;
  } else if (outType_ == UPPER_TRI_MATRIX) {
    out_ = setup.DSL().AddSet(DataSet::MATRIX_FLT, dsname, "HAUSDORFF");
    if (out_ == 0 || ((DataSet_2D*)out_)->AllocateTriangle( nrows )) return Analysis::ERR;
    if (out_->Size() != inputSets_.size())
      mprintf("Warning: Number of input data sets (%zu) != number of expected sets in matrix (%zu)\n",
              inputSets_.size(), out_->Size());
  }
  if (df != 0) df->AddDataSet( out_ );

  mprintf("    HAUSDORFF:\n");
  for (DataSetList::const_iterator it = inputSets_.begin(); it != inputSets_.end(); ++it)
    mprintf("\t%s\n", (*it)->legend());

  return Analysis::OK;
}

// Analysis_HausdorffDistance::Analyze()
Analysis::RetType Analysis_HausdorffDistance::Analyze() {
  // Get the Hausdorff distance for each set.
  double hd = 0.0;
  int idx = 0;
  for (DataSetList::const_iterator it = inputSets_.begin(); it != inputSets_.end(); ++it)
  {
    hd = -1.0;
    if ( (*it)->Group() == DataSet::MATRIX_2D )
      hd = h_Matrix( static_cast<DataSet_2D*>( *it ) );
    else
      mprintf("Warning: '%s' type not yet supported for Hausdorff\n", (*it)->legend());
    mprintf("%12.4f %s\n", hd, (*it)->legend());
    float fhd = (float)hd;
    switch (outType_) {
      case BASIC : out_->Add(idx++, &fhd); break;
      case UPPER_TRI_MATRIX : ((DataSet_MatrixFlt*)out_)->AddElement( fhd ); break;
    }
  }
  return Analysis::OK;
}
