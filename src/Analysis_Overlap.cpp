#include <cmath> // fabs
#include "Constants.h" // SMALL
#include "Analysis_Overlap.h"
#include "CpptrajStdio.h"

Analysis_Overlap::Analysis_Overlap() : ds1_(0), ds2_(0) {}

void Analysis_Overlap::Help() {
  mprintf("\tds1 <ds1> ds2 <ds2>\n");
}

static inline bool check_type(DataSet* ds, int n_ds) {
  if (ds == 0) {
    mprinterr("Error: Data set ds%i not found.\n", n_ds);
    return true;
  }
  if (ds->Type() != DataSet::HIST && 
      ds->Type() != DataSet::FLOAT &&
      ds->Type() != DataSet::DOUBLE &&
      ds->Type() != DataSet::INT) {
    mprinterr("Error: %s: bad set type for overlap.\n", ds->Legend().c_str());
    return true;
  }
  return false;
}

Analysis::RetType Analysis_Overlap::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  // Keywords
  ds1_ = datasetlist->GetDataSet( analyzeArgs.GetStringKey("ds1") );
  if (check_type(ds1_,1)) return Analysis::ERR;
  ds2_ = datasetlist->GetDataSet( analyzeArgs.GetStringKey("ds2") );
  if (check_type(ds2_,2)) return Analysis::ERR;

  mprintf("    OVERLAP: Between %s and %s\n", ds1_->Legend().c_str(),
          ds2_->Legend().c_str());

  return Analysis::OK;
}

Analysis::RetType Analysis_Overlap::Analyze() {
  if (ds1_->Size() < 1 || ds2_->Size() < 1) {
    mprinterr("Error: One or both data sets empty (ds1=%i, ds2=%i)\n",
              ds1_->Size(), ds2_->Size());
    return Analysis::ERR;
  }
  if (ds1_->Size() != ds2_->Size()) {
    mprinterr("Error: Data set sizes do not match (ds1=%i, ds2=%i)\n",
              ds1_->Size(), ds2_->Size());
    return Analysis::ERR;
  }
  int Npoints = 0;
  double sum = 0.0;
  for (int i = 0; i < ds1_->Size(); i++) {
    double val1 = ds1_->Dval(i);
    double val2 = ds2_->Dval(i);
    if (fabs(val1) < SMALL && fabs(val2) < SMALL) {
      // No data in either set, do not process;
      continue;
    }
    double denominator = val1 + val2;
    if (fabs(denominator) < SMALL) {
      // Complete opposite, no overlap, but process
      ++Npoints;
      continue;
    }
    //mprintf("\t%4i %8.3f %8.3f %8.3f %8.3f\n",Npoints,val1,val2,denominator,(1.0 - ( fabs(val1 - val2) / denominator ))); // DEBUG
    sum += (1.0 - ( fabs(val1 - val2) / denominator ));
    ++Npoints;
  }
  if (Npoints < 1)
    sum = 0.0;
  else
    sum /= (double)Npoints;
  mprintf("\t%i of %i points had no data.\n", ds1_->Size() - Npoints, ds1_->Size());
  mprintf("\tPercent overlap between %s and %s is %f\n", ds1_->Legend().c_str(),
          ds2_->Legend().c_str(), sum);
  return Analysis::OK;
}
