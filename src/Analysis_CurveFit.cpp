#include <cmath> // sqrt
#include "Analysis_CurveFit.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include "CurveFit.h"
#include "RPNcalc.h"
#include "DataSet_Mesh.h"

/// The RPN calculator is static so it can be used in Equation.
static RPNcalc Calc_;

/// This function will be passed to CurveFit
int Equation(CurveFit::Darray const& Xvals, CurveFit::Darray const& Params,
                      CurveFit::Darray& Yvals)
{
  for (unsigned int n = 0; n != Xvals.size(); n++)
    Calc_.Evaluate(Params, Xvals[n], Yvals[n]);
  return 0;
}
    
// CONSTRUCTOR
Analysis_CurveFit::Analysis_CurveFit() :
  dset_(0),
  finalY_(0),
  tolerance_(0.0),
  maxIt_(0)
{} 

// Analysis_CurveFit::Help()
void Analysis_CurveFit::Help() {
  mprintf("\t<dset> <equation> [out <outfile>]\n"
          "\t[tol <tolerance>] [maxit <max iterations>]\n"
          "  Fit data set <dset> to <equation>. The equation must have form:\n"
          "    <var> = <expression>\n"
          "  where <expression> can contain variable 'X' and parameters A<n>.\n");
}

Analysis::RetType Analysis_CurveFit::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  // First argument should be DataSet to fit to.
  std::string dsinName = analyzeArgs.GetStringNext();
  dset_ = datasetlist->GetDataSet( dsinName );
  if (dset_ == 0) {
    mprinterr("Error: Data set '%s' not found.\n", dsinName.c_str());
    return Analysis::ERR;
  }
  if (dset_->Ndim() != 1) {
    mprinterr("Error: Curve fitting can only be done with 1D data sets.\n");
    return Analysis::ERR;
  }
  // Second argument should be the equation to fit DataSet to.
  std::string equation = analyzeArgs.GetStringNext();
  if (equation.empty()) {
    mprinterr("Error: Must specify an equation.\n");
    return Analysis::ERR;
  }
  Calc_.SetDebug(debugIn);
  if (Calc_.ProcessExpression( equation )) return Analysis::ERR;
  // Equation must have an assignment.
  if ( Calc_.AssignStatus() != RPNcalc::YES_ASSIGN ) {
    mprinterr("Error: No assignment '=' in equation.\n");
    return Analysis::ERR;
  }
  std::string dsoutName = Calc_.FirstTokenName();
  if (dsoutName.empty()) {
    mprinterr("Error: Invalid assignment in equation.\n");
    return Analysis::ERR;
  } 
  // Get keywords
  DataFile* outfile = DFLin->AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs );
  tolerance_ = analyzeArgs.getKeyDouble("tol", 0.0001);
  if (tolerance_ < 0.0) {
    mprinterr("Error: Tolerance must be greater than or equal to 0.0\n");
    return Analysis::ERR;
  }
  maxIt_ = analyzeArgs.getKeyInt("maxit", 50);
  if (maxIt_ < 1) {
    mprinterr("Error: Max iterations must be greater than or equal to 1.\n");
    return Analysis::ERR;
  }
  // Now get all parameters
  int n_expected_params = Calc_.Nparams();
  if (n_expected_params < 0) return Analysis::ERR;
  
  Params_.resize( n_expected_params, 0.0 );
  for (int p = 0; p != n_expected_params; p++) {
    std::string parameterArg = analyzeArgs.GetStringNext();
    if (parameterArg.empty())
      break;
    ArgList parameter(parameterArg, " =");
    if (parameter.Nargs() != 2) {
      mprinterr("Error: Invalid parameter argument. Expected 'A<n>=<value>'\n");
      return Analysis::ERR;
    }
    std::string parameterName = parameter.GetStringNext();
    if (parameterName[0] != 'A') {
      mprinterr("Error: Invalid parameter name (expected A<n>): %s\n", parameterName.c_str());
      return Analysis::ERR;
    }
    int pnum = convertToInteger(parameterName.substr(1));
    Params_[pnum] = parameter.getNextDouble(0.0);
  }
  // Set up output data set.
  finalY_ = datasetlist->AddSet(DataSet::XYMESH, dsoutName, "FIT");
  if (finalY_ == 0) return Analysis::ERR;
  if (outfile != 0) outfile->AddSet( finalY_ );

  mprintf("    CURVEFIT: Fitting set '%s' to equation '%s'\n",
          dset_->Legend().c_str(), equation.c_str());
  mprintf("\tFinal Y values will be saved in set '%s'\n", finalY_->Legend().c_str());
  mprintf("\tTolerance= %g, maximum iterations= %i\n", tolerance_, maxIt_);
  mprintf("\tInitial parameters:\n");
  for (Darray::const_iterator ip = Params_.begin(); ip != Params_.end(); ++ip)
    mprintf("\t  A%u = %g\n", ip - Params_.begin(), *ip);

  return Analysis::OK;
}

// Analysis_CurveFit::Analyze()
Analysis::RetType Analysis_CurveFit::Analyze() {
  if (dset_->Size() < 1) {
    mprinterr("Error: Set %s is empty.\n", dset_->Legend().c_str());
    return Analysis::ERR;
  }

  // Set up initial Y and X values.
  DataSet_1D& Set = static_cast<DataSet_1D&>( *dset_ );
  CurveFit::Darray Xvals, Yvals;
  Xvals.reserve( Set.Size() );
  Yvals.reserve( Set.Size() );
  bool setHasZero = false;
  for (unsigned int i = 0; i != Set.Size(); i++) {
    Xvals.push_back( Set.Xcrd(i) );
    if (Set.Dval(i) == 0.0)
      setHasZero = true;
    Yvals.push_back( Set.Dval(i) );
  }

  // Perform curve fitting.
  CurveFit fit;
  int info = fit.LevenbergMarquardt( Equation, Xvals, Yvals, Params_, tolerance_, maxIt_ );
  mprintf("\t%s\n", fit.Message(info));
  if (info == 0) {
    mprinterr("Error: %s\n", fit.ErrorMessage());
    return Analysis::ERR;
  }
  for (Darray::const_iterator ip = Params_.begin(); ip != Params_.end(); ++ip)
    mprintf("\t\tFinal Param %u = %g\n", ip - Params_.begin(), *ip);

  // FIXME Should probably have better error handling here.
  // Construct output data.
  DataSet_Mesh& Yout = static_cast<DataSet_Mesh&>( *finalY_ );
  CurveFit::Darray::const_iterator ny = fit.FinalY().begin();
  Yout.Allocate1D( dset_->Size() );
  for (CurveFit::Darray::const_iterator x = Xvals.begin(); x != Xvals.end(); ++x, ++ny)
    Yout.AddXY( *x, *ny );

  // Statistics
  // TODO: Move to CurveFit?
  double corr_coeff = Yout.CorrCoeff( Set );
  mprintf("\tCorrelation coefficient: %g\n", corr_coeff);
  double ChiSq = 0.0;
  double Y2 = 0.0;
  for (unsigned int i = 0; i != Set.Size(); i++) {
    double diff = Yout.Dval(i) - Set.Dval(i);
    ChiSq += (diff * diff);
    Y2 += (Set.Dval(i) * Set.Dval(i));
  }
  double TheilU = sqrt(ChiSq / Y2);
  mprintf("\tChi squared: %g\n", ChiSq);
  mprintf("\tUncertainty coefficient: %g\n", TheilU);
  if (!setHasZero) {
    double rms_percent_error = 0.0;
    for (unsigned int i = 0; i != Set.Size(); i++) {
      double diff = Yout.Dval(i) - Set.Dval(i);
      rms_percent_error += (diff * diff) / (Set.Dval(i) * Set.Dval(i));
    }
    rms_percent_error = sqrt( rms_percent_error / (double)Set.Size() );
    mprintf("\tRMS percent error: %g\n", rms_percent_error);
  } else
    mprintf("Warning: Input set Y values contain zero, cannot calculate RMS percent error\n");
  
  return Analysis::OK;
}
