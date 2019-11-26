#include "Analysis_EvalEquilibration.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"
#include "DataSet_Mesh.h"
#include "CurveFit.h"

Analysis_EvalEquilibration::Analysis_EvalEquilibration() :
  Analysis(HIDDEN),
  statsout_(0),
  tolerance_(0),
  maxIt_(0),
  debug_(0)
{}

// Analysis_EvalEquilibration::Help()
void Analysis_EvalEquilibration::Help() const {
  mprintf("\n");
}

// Analysis_EvalEquilibration::Setup()
Analysis::RetType Analysis_EvalEquilibration::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  debug_ = debugIn;

  dsname_ = analyzeArgs.GetStringKey("name");
  if (dsname_.empty())
    dsname_ = setup.DSL().GenerateDefaultName("EvalEquil");

  tolerance_ = analyzeArgs.getKeyDouble("tol", 0.00001);
  if (tolerance_ < 0.0) {
    mprinterr("Error: Tolerance must be greater than or equal to 0.0\n");
    return Analysis::ERR;
  }
  maxIt_ = analyzeArgs.getKeyInt("maxit", 500);
  if (maxIt_ < 1) {
    mprinterr("Error: Max iterations must be greater than or equal to 1.\n");
    return Analysis::ERR;
  }

  // Curves out
  DataFile* outfile = setup.DFL().AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs );
  // Results out
  DataFile* resultsOut = setup.DFL().AddDataFile( analyzeArgs.GetStringKey("resultsout"), analyzeArgs );

  // Allocate output stats file. Allow STDOUT.
  statsout_ = setup.DFL().AddCpptrajFile( analyzeArgs.GetStringKey("statsout"),
                                          "EvalEquil stats",
                                          DataFileList::TEXT, true );
  if (statsout_ == 0) return Analysis::ERR;

  // get input data sets
  if (inputSets_.AddSetsFromArgs( analyzeArgs.RemainingArgs(), setup.DSL() ))
    return Analysis::ERR;

  // Create output sets.
  int idx = 0;
  for (Array1D::const_iterator it = inputSets_.begin(); it != inputSets_.end(); ++it, ++idx)
  {
    DataSet* setOut = setup.DSL().AddSet( DataSet::XYMESH, MetaData( dsname_, idx ) );
    if (setOut == 0) return Analysis::ERR;
    outputSets_.push_back( setOut );
    if (outfile != 0) {
      outfile->AddDataSet( *it );
      outfile->AddDataSet( setOut );
    }
  }
  dataChiSq_ = setup.DSL().AddSet( DataSet::DOUBLE, MetaData( dsname_, "chisq" ) );
  dataName_  = setup.DSL().AddSet( DataSet::STRING, MetaData( dsname_, "name"  ) );
  DataSet::SizeArray nData(1, inputSets_.size());
  dataChiSq_->Allocate( nData );
  dataName_->Allocate( nData );
  if (resultsOut != 0) {
    resultsOut->AddDataSet( dataChiSq_ );
    resultsOut->AddDataSet( dataName_ );
  }

  mprintf("    EVALEQUILIBRATION: Evaluate equilibration of %zu sets.\n", inputSets_.size());
  mprintf("\tOutput set name: %s\n", dsname_.c_str());
  mprintf("\tTolerance for curve fit: %g\n", tolerance_);
  mprintf("\tMax iterations for curve fit: %i\n", maxIt_);
  if (outfile != 0)
    mprintf("\tFit curve output to '%s'\n", outfile->DataFilename().full());
  mprintf("\tStatistics output to '%s'\n", statsout_->Filename().full());
  if (resultsOut != 0)
    mprintf("\tResults output to '%s'\n", resultsOut->DataFilename().full());

  return Analysis::OK;
}

/** Exponential (high relax to low)
  * A2 + (A0*(exp(A1*X)))
  */
int EQ_relax(CurveFit::Darray const& Xvals, CurveFit::Darray const& Params,
             CurveFit::Darray& Yvals)
{
  double A0 = Params[0];
  double A1 = Params[1];
  double A2 = Params[2];
  for (unsigned int n = 0; n != Xvals.size(); ++n)
    Yvals[n] = A2 + ( A0 * exp( -A1 * Xvals[n] ) );
  return 0;
}

/** Inverse exponential (low relax to high).
  * A2 - (A0*(exp(A1*X)))
  */
int EQ_invRelax(CurveFit::Darray const& Xvals, CurveFit::Darray const& Params,
                CurveFit::Darray& Yvals)
{
  double A0 = Params[0];
  double A1 = Params[1];
  double A2 = Params[2];
  for (unsigned int n = 0; n != Xvals.size(); ++n)
    Yvals[n] = A2 - ( A0 * exp( -A1 * Xvals[n] ) );
  return 0;
}


// Analysis_EvalEquilibration::Analyze()
Analysis::RetType Analysis_EvalEquilibration::Analyze() {
  std::vector<DataSet*>::const_iterator ot = outputSets_.begin();
  for (Array1D::const_iterator it = inputSets_.begin(); it != inputSets_.end(); ++it, ++ot)
  {
    mprintf("\tEvaluating: %s\n", (*it)->legend());
    if (!statsout_->IsStream())
      statsout_->Printf("# %s\n", (*it)->legend());
    DataSet_1D const& DS = static_cast<DataSet_1D const&>( *(*it) );
    // First do a linear fit.
    statsout_->Printf("\t----- Linear Fit -----\n");
    if (DS.Size() < 2) {
      mprintf("Warning: Not enough data in '%s' to evaluate.\n", DS.legend());
      continue;
    }
    double slope, intercept, correl;
    int err = DS.LinearRegression( slope, intercept, correl, statsout_ );
    if (err != 0) {
      mprinterr("Error: Could not perform linear regression fit.\n");
      return Analysis::ERR;
    }

    // Process expression to fit
    //if (calc.ProcessExpression( dsname_ + "=A0*(exp(A1*(1/X)))" )) {
    //  mprinterr("Error: Unable to process equation expression.\n");
    //  return Analysis::ERR;
    //}

    statsout_->Printf("\t----- Nonlinear Fit -----\n");
    // Determine relaxation direction
    CurveFit::FitFunctionType fxn = 0;
    int relaxationDir = 0;
    if (slope < 0) {
      mprintf("\tUsing relaxation form: A0*exp(A1*x)\n");
      relaxationDir = -1;
      fxn = EQ_relax;
    } else if (slope > 0) {
      mprintf("\tUsing inverse relaxation form: A0*exp(A1*(1/x))\n");
      relaxationDir = 1;
      fxn = EQ_invRelax;
    } else {
      // Special case: if slope was exactly 0 (should be rare). Consider this
      // equilibrated.
      mprintf("\tSlope of linear fit is 0.\n");
      continue;
    }

    // Set up initial X and Y values.
    double offset = 0.0; // TODO remove
    CurveFit::Darray Xvals, Yvals;
    Xvals.reserve( DS.Size() );
    Yvals.reserve( DS.Size() );
    for (unsigned int i = 0; i != DS.Size(); i++) {
      double xval = DS.Xcrd(i);
      if (xval <= 0) {
        mprintf("Warning: Ignoring X value <= 0: %g\n", xval);
      } else {
        Xvals.push_back( xval - offset );
        Yvals.push_back( DS.Dval(i) );
      }
    }
/*
    // Set up initial X and Y values. Offset the X values so we start from
    // 1.
    double offset = DS.Xcrd(0) - 1.0;
    mprintf("\tFirst X value is %g, Using an offset of %g\n", DS.Xcrd(0), offset);
    CurveFit::Darray Xvals, Yvals;
    Xvals.reserve( DS.Size() );
    Yvals.reserve( DS.Size() );
    for (unsigned int i = 0; i != DS.Size(); i++) {
      double xval = DS.Xcrd(i);
      Xvals.push_back( xval - offset );
      Yvals.push_back( DS.Dval(i) );
    }
*/
    // Determine the average value of the last half of the data.
    unsigned int halfwayPt = (Yvals.size() / 2);
    double Yavg = 0;
    for (unsigned int hidx = halfwayPt; hidx < Yvals.size(); hidx++)
      Yavg += Yvals[hidx];
    Yavg /= (double)(Yvals.size() - halfwayPt);
    mprintf("\tLast half <Y> = %g\n", Yavg);

    // Set initial guesses for parameters.
    CurveFit::Darray Params(3);
    //Params[0] = 0.01; // TODO long avg minus first?
    //Params[1] = 0.1; // TODO abs slope?
    //Params[2] = intercept;   // TODO long avg?
    Params[0] = Yavg - DS.Dval(0);
    if (Params[0] < 0) Params[0] = -Params[0];
    Params[1] = 0.1;
    Params[2] = Yavg;

/*
    // Set initial guesses for parameters: A0 = intercept, A2 = offset
    // A1 could be slope, but if it is too small this can lead to convergence
    // issues. Just use -1.
    CurveFit::Darray Params(3);
    Params[0] = intercept;
    Params[1] = -1.0;
    // For exponential fit, the A1 param should be < 0
    //Params[1] = slope;
    //if (Params[1] > 0)
    //  Params[1] = -Params[1];
    Params[2] = 0;
*/
    for (CurveFit::Darray::const_iterator ip = Params.begin(); ip != Params.end(); ++ip) {
      statsout_->Printf("\tInitial Param A%li = %g\n", ip - Params.begin(), *ip);
    }

    // Perform curve fitting
    CurveFit fit;
    int info = fit.LevenbergMarquardt( fxn, Xvals, Yvals, Params, tolerance_, maxIt_ );
    mprintf("\t%s\n", fit.Message(info));
    if (info == 0) {
      mprinterr("Error: %s\n", fit.ErrorMessage());
      return Analysis::ERR;
    }
    for (CurveFit::Darray::const_iterator ip = Params.begin(); ip != Params.end(); ++ip) {
      statsout_->Printf("\tFinal Param A%li = %g\n", ip - Params.begin(), *ip);
    }

    // Params[0] = A0 = 
    // Params[1] = A1 = 
    // Params[2] = A2 = 
    // Determine the absolute difference of the long-time estimated value
    // from the average value of the last half of the data.
    double ValA = Yavg - Params[2];
    if (ValA < 0) ValA = -ValA;
    mprintf("\tValA = %g\n", ValA);

    // Create output curve
    DataSet_Mesh& OUT = static_cast<DataSet_Mesh&>( *(*ot) );
    for (unsigned int i = 0; i != Xvals.size(); i++)
      OUT.AddXY( Xvals[i] + offset, fit.FinalY()[i] );

    // Statistics
    double corr_coeff, ChiSq, TheilU, rms_percent_error;
    err = fit.Statistics( Yvals, corr_coeff, ChiSq, TheilU, rms_percent_error);
    if (err != 0) mprintf("Warning: %s\n", fit.Message(err));
    statsout_->Printf("\tCorrelation coefficient: %g\n"
                      "\tChi squared: %g\n"
                      "\tUncertainty coefficient: %g\n"
                      "\tRMS percent error: %g\n",
                      corr_coeff, ChiSq, TheilU, rms_percent_error);

    long int oidx = (it - inputSets_.begin());
    dataChiSq_->Add(oidx, &ChiSq);
    dataName_->Add(oidx, (*it)->legend());

    statsout_->Printf("\n");
  }
  return Analysis::OK;
}
