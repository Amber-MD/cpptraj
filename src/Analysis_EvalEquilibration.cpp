#include "Analysis_EvalEquilibration.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"
#include "DataSet_Mesh.h"
#include "CurveFit.h"

Analysis_EvalEquilibration::Analysis_EvalEquilibration() :
  Analysis(HIDDEN),
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

  DataFile* outfile = setup.DFL().AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs );
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

  mprintf("    EVALEQUILIBRATION: Evaluate equilibration of %zu sets.\n", inputSets_.size());
  mprintf("\tOutput set name: %s\n", dsname_.c_str());
  mprintf("\tTolerance for curve fit: %g\n", tolerance_);
  mprintf("\tMax iterations for curve fit: %i\n", maxIt_);
  if (outfile != 0)
    mprintf("\tFit curve output to '%s'\n", outfile->DataFilename().full());

  return Analysis::OK;
}

// A0*(exp(A1*(X)))
int EQ_relax(CurveFit::Darray const& Xvals, CurveFit::Darray const& Params,
             CurveFit::Darray& Yvals)
{
  double A0 = Params[0];
  double A1 = Params[1];
  for (unsigned int n = 0; n != Xvals.size(); ++n)
    Yvals[n] = A0 * ( exp( A1 * Xvals[n] ) );
  return 0;
}

// A0*(exp(A1*(1/X)))
int EQ_invRelax(CurveFit::Darray const& Xvals, CurveFit::Darray const& Params,
                CurveFit::Darray& Yvals)
{
  double A0 = Params[0];
  double A1 = Params[1];
  for (unsigned int n = 0; n != Xvals.size(); ++n)
    Yvals[n] = A0 * ( exp( A1 * (1/Xvals[n]) ) );
  return 0;
}


// Analysis_EvalEquilibration::Analyze()
Analysis::RetType Analysis_EvalEquilibration::Analyze() {
  CpptrajFile statsout;
  statsout.OpenWrite("");

  std::vector<DataSet*>::const_iterator ot = outputSets_.begin();
  for (Array1D::const_iterator it = inputSets_.begin(); it != inputSets_.end(); ++it, ++ot)
  {
    mprintf("\tEvaluating '%s'\n", (*it)->legend());
    DataSet_1D const& DS = static_cast<DataSet_1D const&>( *(*it) );
    // First do a linear fit.
    if (DS.Size() < 2) {
      mprintf("Warning: Not enough data in '%s' to evaluate.\n", DS.legend());
      continue;
    }
    double slope, intercept, correl;
    int err = DS.LinearRegression( slope, intercept, correl, &statsout );
    if (err != 0) {
      mprinterr("Error: Could not perform linear regression fit.\n");
      return Analysis::ERR;
    }

    // Process expression to fit
    //if (calc.ProcessExpression( dsname_ + "=A0*(exp(A1*(1/X)))" )) {
    //  mprinterr("Error: Unable to process equation expression.\n");
    //  return Analysis::ERR;
    //}

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

    // Set up initial X and Y values
    CurveFit::Darray Xvals, Yvals;
    Xvals.reserve( DS.Size() );
    Yvals.reserve( DS.Size() );
    for (unsigned int i = 0; i != DS.Size(); i++) {
      double xval = DS.Xcrd(i);
      if (xval <= 0) {
        mprintf("Warning: Ignoring X value <= 0: %g\n", xval);
      } else {
        Xvals.push_back( xval );
        Yvals.push_back( DS.Dval(i) );
      }
    }

    // Set initial guesses for parameters: A0 = intercept, A1 = slope
    CurveFit::Darray Params(2);
    Params[0] = intercept;
    Params[1] = slope;

    // Perform curve fitting
    CurveFit fit;
    int info = fit.LevenbergMarquardt( fxn, Xvals, Yvals, Params, tolerance_, maxIt_ );
    mprintf("\t%s\n", fit.Message(info));
    if (info == 0) {
      mprinterr("Error: %s\n", fit.ErrorMessage());
      return Analysis::ERR;
    }
    for (CurveFit::Darray::const_iterator ip = Params.begin(); ip != Params.end(); ++ip) {
      statsout.Printf("\tFinal Param A%li = %g\n", ip - Params.begin(), *ip);
    }

    // Create output curve
    DataSet_Mesh& OUT = static_cast<DataSet_Mesh&>( *(*ot) );
    for (unsigned int i = 0; i != Xvals.size(); i++)
      OUT.AddXY( Xvals[i], fit.FinalY()[i] );

    
  }
  return Analysis::OK;
}
