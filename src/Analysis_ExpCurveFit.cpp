#include <cmath> // exp
#include <cfloat> // DBL_MAX
#include "Analysis_ExpCurveFit.h"
#include "CpptrajStdio.h"
#include "CurveFit.h"
#include "DataSet_Mesh.h"

/** Multi exponential of form Y =  B0 + SUM[Bi * exp(X * Bi+1)] */
int MultiExponentialK(CurveFit::Darray const& Xvals, CurveFit::Darray const& Params,
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
/*
double TestFxn(double X, CurveFit::Darray const& Params)
{
  double Y = Params[0];
  double penalty = Params[0];
  for (unsigned int i = 1; i < Params.size(); i += 2)
  {
    penalty += Params[i];
    double expBx = exp( Params[i+1] * X );
    Y += Params[i] * expBx;
  }
  penalty = 1000.0 * (1.0 - penalty);
  return Y + penalty;
}
*/
/** Multi exponential of form Y = SUM[Bi * exp(X * Bi+1)] */
int MultiExponential(CurveFit::Darray const& Xvals, CurveFit::Darray const& Params,
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

void Analysis_ExpCurveFit::Help() {
  mprintf("\t<dset> [name <output setname>] [out <outfile>] [nexp <n>]\n"
          "\t[tol <tolerance>] [maxit <max iterations>] [useconstant]\n");
}

// Analysis_ExpCurveFit::Setup()
Analysis::RetType Analysis_ExpCurveFit::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  // Get keywords
  DataFile* outfile = DFLin->AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs );
  std::string dsoutName = analyzeArgs.GetStringKey("name");
  useConstant_ = analyzeArgs.hasKey("useconstant");
  usePrefactorBounds_ = analyzeArgs.hasKey("prefactorbounds");
  nexp_ = analyzeArgs.getKeyInt("nexp", 2);
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
  // DataSet name should be first argument
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
  // Set up output data set.
  finalY_ = datasetlist->AddSet(DataSet::XYMESH, dsoutName, "EXPFIT");
  if (finalY_ == 0) return Analysis::ERR;
  if (outfile != 0) outfile->AddSet( finalY_ );

  mprintf("    EXPCURVEFIT: Fitting multi-exponential curve using %i exponentials.\n", nexp_);
  mprintf("\tUsing data in set '%s'.\n", dset_->Legend().c_str());
  mprintf("\tFinal Y values will be saved in set '%s'\n", finalY_->Legend().c_str());
  mprintf("\tFitting to form:");
  if (useConstant_)
    mprintf(" Y = B0 + SUM[Bi * exp(X * Bi+1)\n");
  else
    mprintf(" Y = SUM[Bi * exp(X * Bi+1)\n");
  mprintf("\tTolerance= %g, maximum iterations= %i\n", tolerance_, maxIt_);
  if (usePrefactorBounds_)
    mprintf("\tExponential prefactors will be bound between 0.0 and 1.0.\n");
  if (outfile != 0)
    mprintf("\tFinal Y values will be written in file '%s'\n", outfile->DataFilename().full());

  return Analysis::OK;
}

// Analysis_ExpCurveFit::Analyze()
Analysis::RetType Analysis_ExpCurveFit::Analyze() {
  if (dset_->Size() < 1) {
    mprinterr("Error: Set %s is empty.\n", dset_->Legend().c_str());
    return Analysis::ERR;
  }

  // Set up function to use
  CurveFit::Darray Params;
  CurveFit::FitFunctionType fxn = 0;
  unsigned int pstart = 0;
  if (useConstant_) {
    pstart = 1;
    fxn = MultiExponentialK; 
    Params.resize( 1 + (nexp_ * 2) );
//    if (usePrefactorBounds_) {
//      Params[0] = 1.0 / (double)(Params.size() - pstart);
//      //fxn = TestFxn;
//    } else
      Params[0] = 0.0;
  } else {
    fxn = MultiExponential;
    Params.resize( nexp_ * 2 );
  }
  Params.resize( pstart + (nexp_ * 2) );
  // Set initial parameters and bounds if necessary.
  std::vector<bool> bounds;
  CurveFit::Darray UB, LB, Weights;
  if (usePrefactorBounds_) {
//    bounds.assign( Params.size(), false );
//    LB.assign( Params.size(), 0.0 );
//    UB.assign( Params.size(), 0.0 );
    Weights.assign( dset_->Size(), 0.0 );
    double Dsize = (double)dset_->Size();
    for (unsigned int i = 0; i != dset_->Size(); i++) {
      Weights[i] = (Dsize - (double)i) / Dsize;
      //mprintf("DEBUG:\t\tWeight[%u]= %g\n", i, Weights[i]);
    }
  }
  for (unsigned int j = pstart; j < Params.size(); j += 2)
  {
    double prefac =  1.0;
    double expfac = -1.0;
    if (usePrefactorBounds_) {
      prefac /= (double)(Params.size() - pstart);
      expfac = -0.01;
//      bounds[j] = true;
//      bounds[j+1] = true;
//      LB[j] = 0.0;
//      UB[j] = 1.0;
//      LB[j+1] = -1E10;
//      UB[j+1] = 0.0;
    }
    Params[j  ] = prefac;
    Params[j+1] = expfac;
  }
  for (unsigned int j = 0; j != Params.size(); j++) {
    mprintf("\t\tInitial Param %u = %g", j, Params[j]);
    if (!bounds.empty() && bounds[j])
      mprintf("    Bounds: %g < P%u < %g", LB[j], UB[j]);
    mprintf("\n");
  }

  // Set up X and Y values.
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
  //int info = fit.LevenbergMarquardt( fxn, Xvals, Yvals, Params, tolerance_, maxIt_ );
  int info = fit.LevenbergMarquardt( fxn, Xvals, Yvals, Params, 
                                     bounds, LB, UB, Weights, tolerance_, maxIt_ );
  mprintf("\t%s\n", fit.Message(info));
  if (info == 0) {
    mprinterr("Error: %s\n", fit.ErrorMessage());
    return Analysis::ERR;
  } 
  for (unsigned int j = 0; j != Params.size(); j++)
    mprintf("\t\tFinal Param %u = %g\n", j, Params[j]);

  // Write out final form
  mprintf("\tFinal Eq: Y =");
  if (useConstant_)
    mprintf(" %g +", Params[0]);
  for (unsigned int p = pstart; p < Params.size(); p += 2) {
    if (p > pstart) mprintf(" +");
    mprintf(" (%g * exp(x * %g))", Params[p], Params[p+1]);
  }
  mprintf("\n");

  // TEST: Integrate the function out a ways.
/*
  double sum = 0.0;
  double xstep = Set.Dim(Dimension::X).Step();
  double lastY = fxn( 0.0, Params);// - Params[0];
  double delta = 0.0;
  double lastslice = lastY;
  unsigned int intmax = 1000000;
  for (unsigned int i = 1; i < intmax; i++) {
    double newY = fxn( (double)i * xstep, Params);// - Params[0];
    double slice = (xstep * (lastY + newY) * 0.5);
    //mprintf("DBG: %g\n", slice);
    delta = lastslice - slice;
    sum += slice;
    lastslice = slice;
  }
  mprintf("\tDEBUG: Integration from 0.0 to %g is %g (delta %g)\n", intmax*xstep, sum, delta);
*/
  // FIXME Should probably have better error handling here.
  // Construct output data.
  DataSet_Mesh& Yout = static_cast<DataSet_Mesh&>( *finalY_ );
  CurveFit::Darray::const_iterator ny = fit.FinalY().begin();
  Yout.Allocate1D( dset_->Size() );
  for (CurveFit::Darray::const_iterator x = Xvals.begin(); x != Xvals.end(); ++x, ++ny)
    Yout.AddXY( *x, *ny );
    //Yout.AddXY( *x, fxn(*x, Params) );

  // Statistics
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

  // Report decay constants from exponential parameters.
  mprintf("\tTime constants:");
  for (unsigned int p = pstart + 1; p < Params.size(); p += 2) {
    if (Params[p] != 0.0) {
      double tau = 1.0 / -Params[p];
      mprintf(" %12.6g", tau);
    } else
      mprintf("Warning: exp parameter %u is zero.\n", ((p - pstart) / 2) + 1);
  }
  mprintf("\n");

  return Analysis::OK;
}
