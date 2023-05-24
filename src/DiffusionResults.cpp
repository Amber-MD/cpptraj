#include "DiffusionResults.h"
#include "CpptrajStdio.h"
#include "DataSet.h"
#include "DataSet_1D.h"
#include "DataSetList.h"
#include "DataFileList.h"
#include "ArgList.h"

using namespace Cpptraj;
/** CONSTRUCTOR */
DiffusionResults::DiffusionResults() :
  diffConst_(0),
  diffLabel_(0),
  diffSlope_(0),
  diffInter_(0),
  diffCorrl_(0),
  diffout_(0)
{}

/** Process arguments and create sets. */
int DiffusionResults::InitDiffusionResults(DataSetList& DSL, DataFileList& DFL,
                                           ArgList& actionArgs, std::string const& dsname_)
{
  diffout_ = DFL.AddDataFile(actionArgs.GetStringKey("diffout"));
  MetaData::tsType ts = MetaData::NOT_TS;
  diffConst_ = DSL.AddSet(DataSet::DOUBLE, MetaData(dsname_, "D", ts));
  diffLabel_ = DSL.AddSet(DataSet::STRING, MetaData(dsname_, "Label", ts));
  diffSlope_ = DSL.AddSet(DataSet::DOUBLE, MetaData(dsname_, "Slope", ts));
  diffInter_ = DSL.AddSet(DataSet::DOUBLE, MetaData(dsname_, "Intercept", ts));
  diffCorrl_ = DSL.AddSet(DataSet::DOUBLE, MetaData(dsname_, "Corr", ts));
  if (diffConst_ == 0 || diffLabel_ == 0 || diffSlope_ == 0 || diffInter_ == 0 ||
      diffCorrl_ == 0)
  {
    mprinterr("Error: Could not allocate 1 or more diffusion constant sets.\n");
    return 1;
  }
# ifdef MPI
  // No sync needed since these are not time series
  diffConst_->SetNeedsSync( false );
  diffLabel_->SetNeedsSync( false );
  diffSlope_->SetNeedsSync( false );
  diffInter_->SetNeedsSync( false );
  diffCorrl_->SetNeedsSync( false );
# endif
  if (diffout_ != 0) {
    diffout_->AddDataSet( diffConst_ );
    diffout_->AddDataSet( diffSlope_ );
    diffout_->AddDataSet( diffInter_ );
    diffout_->AddDataSet( diffCorrl_ );
    diffout_->AddDataSet( diffLabel_ );
  }
  Dimension Ddim( 1, 1, "Set" );
  diffConst_->SetDim(Dimension::X, Ddim);
  diffLabel_->SetDim(Dimension::X, Ddim);
  diffSlope_->SetDim(Dimension::X, Ddim);
  diffInter_->SetDim(Dimension::X, Ddim);
  diffCorrl_->SetDim(Dimension::X, Ddim);
  return 0;
}

/** Print info to stdout. */
void DiffusionResults::Info() const {
  if (diffout_ != 0)
    mprintf("\tDiffusion constant output to '%s'\n", diffout_->DataFilename().full());
  else
    mprintf("\tDiffusion constant output to STDOUT.\n");
}

/** Calculate diffusion constant from slope of MSD vs time data. */
void DiffusionResults::CalcDiffusionConst(unsigned int& set, DataSet* ds, int Ndim,
                                          std::string const& label) const
{ // TODO check 1D
  DataSet_1D const& data = static_cast<DataSet_1D const&>( *ds );
  double Factor = 10.0 / ((double)Ndim * 2.0);
  double slope, intercept, corr;
  double Dval = 0.0;
  double Fval = 0;
  if (data.LinearRegression( slope, intercept, corr, Fval, 0 ) == 0)
    Dval = slope * Factor;
  if (diffout_ == 0)
    mprintf("\t'%s' D= %g  Slope= %g  Int= %g  Corr= %g\n", data.legend(), Dval,
            slope, intercept, corr);
  diffConst_->Add(set  , &Dval);
  diffSlope_->Add(set  , &slope);
  diffInter_->Add(set  , &intercept);
  diffCorrl_->Add(set  , &corr);
  diffLabel_->Add(set++, label.c_str());
}
