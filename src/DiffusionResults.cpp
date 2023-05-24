#include "DiffusionResults.h"
#include "CpptrajStdio.h"
#include "DataSet.h"
#include "DataSet_1D.h"

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
