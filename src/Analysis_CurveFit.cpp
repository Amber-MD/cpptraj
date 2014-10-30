#include <cmath> // sqrt
#include "Analysis_CurveFit.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include "CurveFit.h"
#include "RPNcalc.h"
#include "DataSet_Mesh.h"

// -----------------------------------------------------------------------------
/// The RPN calculator is static so it can be used in Equation.
static RPNcalc Calc_;

/// This is for generic equations with RPNcalc.
int Equation(CurveFit::Darray const& Xvals, CurveFit::Darray const& Params,
                      CurveFit::Darray& Yvals)
{
  for (unsigned int n = 0; n != Xvals.size(); n++)
    Calc_.Evaluate(Params, Xvals[n], Yvals[n]);
  return 0;
}

/// Multi exponential of form Y = SUM[Bi * exp(X * Bi+1)] 
int EQ_MultiExp(CurveFit::Darray const& Xvals, CurveFit::Darray const& Params,
                CurveFit::Darray& Yvals)
{
  for (unsigned int n = 0; n != Xvals.size(); n++) {
    double X = Xvals[n];
    double Y = 0.0;
    for (unsigned int i = 0; i < Params.size(); i += 2)
    {
      double expBx = exp( Params[i+1] * X );
      Y += Params[i] * expBx;
      //dYdP[i  ] = expBx;
      //dYdP[i+1] = Params[i] * X * expBx;
      //printf("DEBUG: MultiExponential: dYdP[%i]= %g\n", i, dYdP[i]);
      //printf("DEBUG: MultiExponential: dYdP[%i]= %g\n", i+1, dYdP[i+1]);
    }
    Yvals[n] = Y;
  }
  return 1;
}

/// Multi exponential of form Y =  B0 + SUM[Bi * exp(X * Bi+1)]
int EQ_MultiExpK(CurveFit::Darray const& Xvals, CurveFit::Darray const& Params,
                 CurveFit::Darray& Yvals)
{
  for (unsigned int n = 0; n != Xvals.size(); n++) {
    double X = Xvals[n];
    double Y = Params[0];
    // dYdP[0] = 1.0; 
    for (unsigned int i = 1; i < Params.size(); i += 2)
    {
      double expBx = exp( Params[i+1] * X );
      Y += Params[i] * expBx;
      //dYdP[i  ] = expBx;
      //dYdP[i+1] = Params[i] * X * expBx;
      //printf("DEBUG: MultiExponential: dYdP[%i]= %g\n", i, dYdP[i]);
      //printf("DEBUG: MultiExponential: dYdP[%i]= %g\n", i+1, dYdP[i+1]);
    }
    Yvals[n] = Y;
  }
  return 0;
}

// -----------------------------------------------------------------------------
// CONSTRUCTOR
Analysis_CurveFit::Analysis_CurveFit() :
  dset_(0),
  finalY_(0),
  tolerance_(0.0),
  maxIt_(0),
  nexp_(-1),
  eqForm_(GENERAL)
{} 

// Analysis_CurveFit::Help()
void Analysis_CurveFit::Help() {
  mprintf("\t<dset> {<equation> | nexp <n> [form {mexp|mexpk}} [out <outfile>]\n"
          "\t[tol <tolerance>] [maxit <max iterations>]\n"
          "  Fit data set <dset> to <equation>. The equation must have form:\n"
          "    <var> = <expression>\n"
          "  where <expression> can contain variable 'X' and parameters A<n>.\n"
          "  Alternatively, multi-exponential equations can be used via 'nexp' and 'form'\n");
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
  std::string dsoutName, equation;
  int n_expected_params = 0;
  // Determine if special equation is being used.
  nexp_ = analyzeArgs.getKeyInt("nexp", -1);
  if ( nexp_ > 0 ) {
    // Multi-exponential specialized equation form.
    dsoutName = analyzeArgs.GetStringKey("name");
    if (dsoutName.empty()) {
      mprinterr("Error: 'name <OutputSetName>' must be used with 'nexp <n>'\n");
      return Analysis::ERR;
    }
    equation = dsoutName + " = ";
    // Determine form
    std::string formStr = analyzeArgs.GetStringKey("form");
    if (formStr == "mexp") eqForm_ = MEXP;
    else if (formStr == "mexpk") eqForm_ = MEXP_K;
    else eqForm_ = MEXP;
    // Set up equation
    int nparam = 0;
    if (eqForm_ != MEXP) {
      equation.append("A0 +");
      nparam = 1;
    }
    for (int ie = 0; ie != nexp_; ie++, nparam += 2) {
      if (ie > 0)
        equation.append(" + ");
      equation.append("(A" + integerToString(nparam) + " * exp(X * A" + 
                             integerToString(nparam+1) + "))");
    }
    n_expected_params = nparam;
  } else {
    // Any equation form, solve with RPNcalc.
    eqForm_ = GENERAL;
    // Second argument should be the equation to fit DataSet to.
    equation = analyzeArgs.GetStringNext();
    if (equation.empty()) {
      mprinterr("Error: Must specify an equation if 'nexp <n>' not specified.\n");
      return Analysis::ERR;
    }
    Calc_.SetDebug(debugIn);
    if (Calc_.ProcessExpression( equation )) return Analysis::ERR;
    // Equation must have an assignment.
    if ( Calc_.AssignStatus() != RPNcalc::YES_ASSIGN ) {
      mprinterr("Error: No assignment '=' in equation.\n");
      return Analysis::ERR;
    }
    dsoutName = Calc_.FirstTokenName();
    if (dsoutName.empty()) {
      mprinterr("Error: Invalid assignment in equation.\n");
      return Analysis::ERR;
    } 
    n_expected_params = Calc_.Nparams();
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
  if (nexp_ > 0)
    mprintf("\tMulti-exponential form with %i exponentials.\n", nexp_);
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

  // Set up function to use
  CurveFit::FitFunctionType fxn = 0;
  switch (eqForm_) {
    case GENERAL: fxn = Equation; break;
    case MEXP_K:
//      pstart = 1;
      fxn = EQ_MultiExpK;
      break;
//    case MEXP_K_PENALTY:
//      pstart = 1;
//      fxn = MultiExpK_WithPenalty;
//      break;
    case MEXP: fxn = EQ_MultiExp; break;
    default: return Analysis::ERR;
  }

  // TODO: Set initial parameters and bounds if necessary.


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
  int info = fit.LevenbergMarquardt( fxn, Xvals, Yvals, Params_, tolerance_, maxIt_ );
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

  if (eqForm_ != GENERAL) {
    // Report decay constants from exponential parameters.
    unsigned int pstart = 0;
    if (eqForm_ != MEXP)
      pstart = 1;
    mprintf("\tTime constants:");
    for (unsigned int p = pstart + 1; p < Params_.size(); p += 2) {
      if (Params_[p] != 0.0) {
        double tau = 1.0 / -Params_[p];
        mprintf(" %12.6g", tau);
      } else
        mprintf("Warning: exp parameter %u is zero.\n", ((p - pstart) / 2) + 1);
    }
    mprintf("\n");
  }
  
  return Analysis::OK;
}
