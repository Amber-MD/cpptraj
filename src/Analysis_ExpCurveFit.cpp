#include <cmath> // exp
#include "Analysis_ExpCurveFit.h"
#include "CpptrajStdio.h"
#include "CurveFit.h"
#include "DataSet_Mesh.h"

/** Multi exponential of form Y =  B0 + SUM[Bi * exp(X * Bi+1)] */
double MultiExponentialK(double X, CurveFit::Darray const& Params)
{
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
  return Y;
}

/** Multi exponential of form Y = SUM[Bi * exp(X * Bi+1)] */
double MultiExponential(double X, CurveFit::Darray const& Params)
{
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
  return Y;
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

  mprintf("    EXPCURVEFIT: Fitting multi-exponential curve using %i exponentials\n"
          "\tto data in set '%s'. Final Y values will be saved in set '%s'\n",
          nexp_, dset_->Legend().c_str(), finalY_->Legend().c_str());
  mprintf("\tFitting to form:");
  if (useConstant_)
    mprintf("Y = B0 + SUM[Bi * exp(X * Bi+1)\n");
  else
    mprintf("Y = SUM[Bi * exp(X * Bi+1)\n");
  mprintf("\tTolerance= %g, maximum iterations= %i\n", tolerance_, maxIt_);
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

  // Set up function to use and initial params.
  CurveFit::Darray Params;
  CurveFit::FitFunctionType fxn = 0;
  if (useConstant_) {
    fxn = MultiExponentialK; 
    Params.resize( 1 + (nexp_ * 2) );
    Params[0] = 0.0;
    for (unsigned int j = 1; j < Params.size(); j += 2) {
      Params[j  ] =  1.0;
      Params[j+1] = -1.0;
    }
  } else {
    fxn = MultiExponential;
    Params.resize( nexp_ * 2 );
    for (unsigned int j = 0; j < Params.size(); j += 2) {
      Params[j  ] =  1.0;
      Params[j+1] = -1.0;
    }
  }

  // Set up X and Y values.
  DataSet_1D& Set = static_cast<DataSet_1D&>( *dset_ );
  CurveFit::Darray Xvals, Yvals;
  Xvals.reserve( Set.Size() );
  Yvals.reserve( Set.Size() );
  for (unsigned int i = 0; i != Set.Size(); i++) {
    Xvals.push_back( Set.Xcrd(i) );
    Yvals.push_back( Set.Dval(i) );
  }

  // Perform curve fitting.
  CurveFit fit;
  fit.LevenbergMarquardt( fxn, Xvals, Yvals, Params, tolerance_, maxIt_ );

  // Write out final form
  mprintf("\tFinal Eq: Y =");
  unsigned int pstart = 0;
  if (useConstant_) {
    pstart = 1;
    mprintf(" %g +", Params[0]);
  }
  for (unsigned int p = pstart; p < Params.size(); p += 2) {
    if (p > pstart) mprintf(" +");
    mprintf(" [%g * exp(X * %g)]", Params[p], Params[p+1]);
  }
  mprintf("\n");

  // FIXME Should probably have better error handling here.
  // Construct output data.
  DataSet_Mesh& Yout = static_cast<DataSet_Mesh&>( *finalY_ );
  Yout.Allocate1D( dset_->Size() );
  for (CurveFit::Darray::const_iterator x = Xvals.begin(); x != Xvals.end(); x++)
    Yout.AddXY( *x, fxn(*x, Params) );
  
  return Analysis::OK;
}
