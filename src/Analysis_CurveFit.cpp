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

/** Multi exponential of form Y = B0 + SUM[Bi * exp(X * Bi+1)] subject to the
  * constraints that B0 + {Bi} = 1.0 and {Bi+1} < 0.0
  */
int EQ_MultiExpK_Penalty(CurveFit::Darray const& Xvals, CurveFit::Darray const& Params,
                         CurveFit::Darray& Yvals)
{
  // Calculate penalty for coefficients
  double penalty1 = Params[0]; 
  for (unsigned int i = 1; i < Params.size(); i += 2)
    penalty1 += Params[i]; 
  penalty1 = 1000.0 * (1.0 - penalty1);
  // Calculate penalty for exponents
  double pvalue = 1000.0 / (double)((Params.size() - 1) / 2);
  double penalty2 = 0.0;
  for (unsigned int i = 2; i < Params.size(); i += 2)
    if (Params[i] > 0.0) penalty2 += pvalue;
  // Calculate Y values 
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
    Yvals[n] = Y + penalty1 + penalty2;
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
  mprintf("\t<dset> {<equation> | nexp <m> [form {mexp|mexpk|mexpk_penalty}}\n"
          "\t[out <outfile>] [tol <tolerance>] [maxit <max iterations>]\n"
          "  Fit data set <dset> to <equation>. The equation must have form:\n"
          "    <var> = <expression>\n"
          "  where <var> is the output data set name and <expression> can contain\n"
          "  variable 'X' and parameters A<n>.\n"
          "  Alternatively, multi-exponential equations can be used via 'nexp' and 'form':\n"
          "    mexp:  SUM(m)[ An * exp(An+1 * X)]\n"
          "    mexpk: A0 + SUM(m)[An * exp(An+1 * X)]\n"
          "    mexpk_penalty: Same as mexpk except sum of prefactors constrained to 1.0 and\n"
          "                   exp. constants constrained to < 0.0.\n");
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
    eqForm_ = MEXP;
    std::string formStr = analyzeArgs.GetStringKey("form");
    if (!formStr.empty()) {
      if (formStr == "mexp") eqForm_ = MEXP;
      else if (formStr == "mexpk") eqForm_ = MEXP_K;
      else if (formStr == "mexpk_penalty") eqForm_ = MEXP_K_PENALTY;
      else {
        mprinterr("Error: Multi-exponential form '%s' not recognized.\n", formStr.c_str());
        return Analysis::ERR;
      }
    }
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
      mprinterr("Error: No assignment '=' in equation '%s'.\n", equation.c_str());
      return Analysis::ERR;
    }
    dsoutName = Calc_.FirstTokenName();
    if (dsoutName.empty()) {
      mprinterr("Error: Invalid assignment in equation '%s'.\n", equation.c_str());
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
  int n_specified_params = 0;
  for (int p = 0; p != n_expected_params; p++, n_specified_params++) {
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
  // Check if all params specified.
  if (n_specified_params != n_expected_params) {
    mprintf("Warning: # specified params (%i) less than # expected params (%i)\n",
            n_specified_params, n_expected_params);
    if (eqForm_ != GENERAL && n_specified_params == 0) {
      mprintf("Warning: For multi-exponential using default params.\n");
      int pnum = 0;
      if (eqForm_ != MEXP) {
        pnum = 1;
        Params_[0] = 1.0;
      }
      for (int p = pnum; p != n_expected_params; p+=2) {
        Params_[p] = 1.0;
        Params_[p+1] = -1.0;
      }
    }
  }
  // Set up output data set.
  finalY_ = datasetlist->AddSet(DataSet::XYMESH, dsoutName, "FIT");
  if (finalY_ == 0) return Analysis::ERR;
  if (outfile != 0) outfile->AddSet( finalY_ );

  mprintf("    CURVEFIT: Fitting set '%s' to equation '%s'\n",
          dset_->Legend().c_str(), equation.c_str());
  if (nexp_ > 0) {
    mprintf("\tMulti-exponential form with %i exponentials.\n", nexp_);
    if (eqForm_ == MEXP_K_PENALTY)
      mprintf("\tMulti-exponential equation constraints: sum of prefactors = 1.0,"
              " exponent parameters < 0.0\n");
  }
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
    case GENERAL:        fxn = Equation; break;
    case MEXP_K:         fxn = EQ_MultiExpK; break;
    case MEXP_K_PENALTY: fxn = EQ_MultiExpK_Penalty; break;
    case MEXP:           fxn = EQ_MultiExp; break;
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

  // Stats specific to multi-exp forms.
  if (eqForm_ != GENERAL) {
    unsigned int pstart = 0;
    if (eqForm_ != MEXP)
      pstart = 1;
    // Report Sum of prefactors
    double sumB = 0.0;
    if (eqForm_ != MEXP)
      sumB = Params_[0];
    for (unsigned int p = pstart; p < Params_.size(); p += 2)
      sumB += Params_[p];
    mprintf("\tSum of prefactors = %g\n", sumB);
    // Report decay constants from exponential parameters.
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
