#include "Analysis_EvalPlateau.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"
#include "DataSet_Mesh.h"
#include "CurveFit.h"
#include "StringRoutines.h" // convertToInteger

Analysis_EvalPlateau::Analysis_EvalPlateau() :
  statsout_(0),
  tolerance_(0),
  initpct_(0),
  finalpct_(0),
  valaCut_(0),
  chisqCut_(0),
  slopeCut_(0),
  maxIt_(0),
  debug_(0)
{
}

/** Aspects for output data sets (data_[]). KEEP IN SYNC WITH OdataType */
const char* Analysis_EvalPlateau::OdataStr_[] = {
  "A0",     // A0
  "A1",     // A1
  "A2",     // A2
  "OneA1",  // ONEA1
  //"F",      // FVAL
  "corr",   // CORR
  "vala",   // VALA
  //"lchisq", // LCHISQ
  "chisq",  // CHISQ
  "pltime", // PLTIME
  "fslope", // FSLOPE
  "name",   // NAME
  "result"  // RESULT
};

/** Set type for output data sets (data_[]). KEEP IN SYNC WITH OdataType */
DataSet::DataType Analysis_EvalPlateau::OdataType_[] = {
  DataSet::DOUBLE, // A0
  DataSet::DOUBLE, // A1
  DataSet::DOUBLE, // A2
  DataSet::DOUBLE, // ONEA1
  //DataSet::DOUBLE, // FVAL
  DataSet::DOUBLE, // CORR
  DataSet::DOUBLE, // VALA
  //DataSet::DOUBLE, // LCHISQ
  DataSet::DOUBLE, // CHISQ
  DataSet::DOUBLE, // PLTIME
  DataSet::DOUBLE, // FSLOPE
  DataSet::STRING, // NAME
  DataSet::STRING  // RESULT
};

// Analysis_EvalPlateau::Help()
void Analysis_EvalPlateau::Help() const {
  mprintf("\t[name <set out name>] [tol <tol>] [valacut <valacut>]\n"
          "\t[initpct <initial pct>] [finalpct <final pct>]\n"
          "\t[chisqcut <chisqcut>] [slopecut <slopecut>] [maxit <maxit>]\n"
          "\t[out <outfile>] [resultsout <resultsfile>] [statsout <statsfile>]\n"
          "\t<input set args> ...\n"
          "  Evaluate whether the input data sets have reached a plateau after\n"
          "  fitting to a single exponential.\n");
}

// Analysis_EvalPlateau::Setup()
Analysis::RetType Analysis_EvalPlateau::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
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
  initpct_ = analyzeArgs.getKeyDouble("initpct", 1.0) / 100.0;
  if (initpct_ <= 0) {
    mprinterr("Error: Initial percent must be greater than 0.\n");
    return Analysis::ERR;
  }
  finalpct_ = analyzeArgs.getKeyDouble("finalpct", 50.0) / 100.0;
  if (finalpct_ <= 0) {
    mprinterr("Error: Final percent must be greater than 0.\n");
    return Analysis::ERR;
  }
  valaCut_ = analyzeArgs.getKeyDouble("valacut", 0.01);
  if (valaCut_ <= 0) {
    mprinterr("Error: valacut must be > 0\n");
    return Analysis::ERR;
  }
  chisqCut_ = analyzeArgs.getKeyDouble("chisqcut", 0.5);
  if (chisqCut_ <= 0) {
    mprinterr("Error: chisqcut must be > 0\n");
    return Analysis::ERR;
  }
  slopeCut_ = analyzeArgs.getKeyDouble("slopecut", 0.000001);
  if (slopeCut_ <= 0) {
    mprinterr("Error: slopecut must be > 0\n");
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
    setOut->SetupFormat() = TextFormat(TextFormat::SCIENTIFIC, 16, 8);
    outputSets_.push_back( setOut );
    if (outfile != 0) {
      outfile->AddDataSet( *it );
      outfile->AddDataSet( setOut );
    }
  }
  // Results data sets
  data_.reserve( inputSets_.size() );
  DataSet::SizeArray nData(1, inputSets_.size());
  for (int idx = 0; idx != (int)NDATA; idx++)
  {
    DataSet* ds = setup.DSL().AddSet( OdataType_[idx], MetaData( dsname_, OdataStr_[idx] ) );
    if (ds == 0) return Analysis::ERR;
    ds->Allocate( nData );
    if (resultsOut != 0)
      resultsOut->AddDataSet( ds );
    // Set formatting
    if ((OdataType)idx == A1 ||
        (OdataType)idx == VALA ||
        (OdataType)idx == FSLOPE ||
        (OdataType)idx == CHISQ //||
        //(OdataType)idx == LCHISQ
       )
      ds->SetupFormat().SetFormatType(TextFormat::GDOUBLE);
    data_.push_back( ds );
  }

  mprintf("    EVALPLATEAU: Evaluate plateau time of %zu sets.\n", inputSets_.size());
  mprintf("\tOutput set name: %s\n", dsname_.c_str());
  mprintf("\tTolerance for curve fit: %g\n", tolerance_);
  mprintf("\tWill use initial %g%% of data for initial guess of A0\n", initpct_*100.0);
  mprintf("\tWill use final %g%% of data for initial guess of A2\n", finalpct_*100.0);
  mprintf("\tMax iterations for curve fit: %i\n", maxIt_);
  if (outfile != 0)
    mprintf("\tFit curve output to '%s'\n", outfile->DataFilename().full());
  mprintf("\tStatistics output to '%s'\n", statsout_->Filename().full());
  if (resultsOut != 0)
    mprintf("\tResults output to '%s'\n", resultsOut->DataFilename().full());
  mprintf("\tCutoff for last half average vs estimated long term value: %g\n", valaCut_);
  mprintf("\tCutoff for non-linear fit chi^2: %g\n", chisqCut_);
  mprintf("\tCutoff for slope: %g\n", slopeCut_);

  return Analysis::OK;
}

/** Exponential decay with time constant A1 from an initial value A0
  * to a plateau A2.
  */
int EQ_plateau(CurveFit::Darray const& Xvals, CurveFit::Darray const& Params,
               CurveFit::Darray& Yvals)
{
  double A0 = Params[0];
  double A1 = Params[1];
  double A2 = Params[2];
  for (unsigned int n = 0; n != Xvals.size(); ++n)
    Yvals[n] = A0 + ((A2 - A0) * (1 - exp(-A1 * Xvals[n])));
  return 0;
}

/** Determine the slope of the given mesh data set and whether it satisfies the cutoff.
  * \return true if slope cutoff is satisfied.
  * \param OUT Set to calculate slope of.
  * \param pltime Set to X value when slope cutoff satisfied.
  * \param finalslope Set to final slope value.
  */
bool Analysis_EvalPlateau::CheckSlope(DataSet_Mesh const& OUT,
                                      double& pltime, double& finalslope)
{
  DataSet_1D::Darray slopeX, slopeY;
  OUT.FiniteDifference(DataSet_1D::FORWARD, slopeX, slopeY);
  double finaly = 0;
  bool slopeCutSatisfied = false;
  pltime = -1;
  finalslope = 0; 
  for (unsigned int sidx = 0; sidx != slopeY.size(); ++sidx) {
    //finalslope = slopeY[sidx];
    //mprintf("DEBUG: slope %g %g\n", slopeX[sidx], slopeY[sidx]);
    double absSlope = slopeY[sidx];
    if (absSlope < 0) absSlope = -absSlope;
    finalslope = absSlope;
    if (!slopeCutSatisfied && absSlope < slopeCut_) {
      pltime = slopeX[sidx];
      finaly = slopeY[sidx];
      slopeCutSatisfied = true;
    }
  }
  if (slopeCutSatisfied)
    statsout_->Printf("\tSlope cutoff satisfied at %g %g\n", pltime, finaly);
  else
    statsout_->Printf("\tSlope cutoff not satisfied.\n");
  return slopeCutSatisfied;
}

/** This is used to add a non-result whenever there is an error 
  * evaluating the input data.
  */
void Analysis_EvalPlateau::BlankResult(long int oidx, const char* legend)
{
  static const double ZERO = 0.0;
  data_[A0]->Add(oidx, &ZERO);
  data_[A1]->Add(oidx, &ZERO);
  data_[A2]->Add(oidx, &ZERO);
  data_[ONEA1]->Add(oidx, &ZERO);
  //data_[FVAL]->Add(oidx, &ZERO);
  data_[CORR]->Add(oidx, &ZERO);
  data_[VALA]->Add(oidx, &ZERO);
  //data_[LCHISQ]->Add(oidx, &ZERO);
  data_[CHISQ]->Add(oidx, &ZERO);
  data_[PLTIME]->Add(oidx, &ZERO);
  data_[FSLOPE]->Add(oidx, &ZERO);
  data_[NAME]->Add(oidx, legend);
  data_[RESULT]->Add(oidx, "err");
}

// Analysis_EvalPlateau::Analyze()
Analysis::RetType Analysis_EvalPlateau::Analyze() {
  std::vector<DataSet*>::const_iterator ot = outputSets_.begin();
  for (Array1D::const_iterator it = inputSets_.begin(); it != inputSets_.end(); ++it, ++ot)
  {
    mprintf("\tEvaluating: %s\n", (*it)->legend());
    if (!statsout_->IsStream())
      statsout_->Printf("# %s\n", (*it)->legend());
    DataSet_1D const& DS = static_cast<DataSet_1D const&>( *(*it) );
    // oidx will be used to index the sets in data_
    long int oidx = (it - inputSets_.begin());
    int err = 0;
    // -------------------------------------------
/*
    // First do a linear fit. TODO may not need this anymore
    statsout_->Printf("\t----- Linear Fit -----\n");
    if (DS.Size() < 2) {
      mprintf("Warning: Not enough data in '%s' to evaluate.\n", DS.legend());
      BlankResult(oidx, DS.legend());
      continue;
    }
    double slope, intercept, correl, Fval;
    err = DS.LinearRegression( slope, intercept, correl, Fval, statsout_ );
    if (err != 0) {
      mprintf("Warning: Could not perform linear regression fit for '%s'.\n", DS.legend());
      BlankResult(oidx, DS.legend());
      continue;
    }
    // Do a chi-squared for the linear fit (residual sum sq.)
    double linChiSq = 0;
    double linAvgr = 0;
    for (unsigned int n = 0; n != DS.Size(); n++) {
      double linY = (slope*DS.Xcrd(n)) + intercept;
      //mprintf("DBG: linear X= %g Y= %g\n", DS.Xcrd(n), linY);
      double diff = linY - DS.Dval(n);
      linChiSq += (diff * diff);
      linAvgr += diff;
    }
    linAvgr /= (double)DS.Size();
    statsout_->Printf("\tAvgResidual= %g\n", linAvgr);
    //statsout_->Printf("\tlinear chisq= %g\n", linChiSq);
*/

    // -------------------------------------------
    // Process expression to fit
    //if (calc.ProcessExpression( dsname_ + "=A0*(exp(A1*(1/X)))" )) {
    //  mprinterr("Error: Unable to process equation expression.\n");
    //  return Analysis::ERR;
    //}

    statsout_->Printf("\t----- Nonlinear Fit -----\n");
    // Determine general relaxation direction
    CurveFit::FitFunctionType fxn = EQ_plateau;

    fxn = EQ_plateau;

    // Set up initial X and Y values.
    double offset = DS.Xcrd(0);
    CurveFit::Darray Xvals, Yvals;
    Xvals.reserve( DS.Size() );
    Yvals.reserve( DS.Size() );
    for (unsigned int i = 0; i != DS.Size(); i++) {
      double xval = DS.Xcrd(i);
      if (xval < 0) {
        mprintf("Warning: Ignoring X value < 0: %g\n", xval);
      } else {
        Xvals.push_back( xval - offset );
        Yvals.push_back( DS.Dval(i) );
      }
    }

    // Get the average of the first initpct_% of data.
    // Will be used as the initial guess for A0, initial density.
    unsigned int initPt = (double)Yvals.size() * initpct_;
    if (initPt < 1) initPt = 1;
    double YinitialAvg = 0;
    for (unsigned int hidx = 0; hidx < initPt; hidx++)
      YinitialAvg += Yvals[hidx];
    YinitialAvg /= (double)initPt;
    statsout_->Printf("\tAvg of first %g%% of the data: %g\n", initpct_*100.0, YinitialAvg);
    // Get the average of the last finalpct_% of data.
    // Will be used as the intial guess for A2, final density.
    unsigned int finalPt = (double)Yvals.size() * finalpct_;
    if (finalPt < 1) finalPt = 1;
    unsigned int finalStart = Yvals.size() - finalPt;
    double YfinalAvg = 0;
    for (unsigned int hidx = finalStart; hidx < Yvals.size(); hidx++)
      YfinalAvg += Yvals[hidx];
    YfinalAvg /= (double)finalPt;
    statsout_->Printf("\tAvg of last %g%% of the data: %g\n", finalpct_*100.0, YfinalAvg);
/*    // Determine the average value of the last half of the data.
    unsigned int halfwayPt = (Yvals.size() / 2);
    double Yavg1half = 0;
    for (unsigned int hidx = 0; hidx < halfwayPt; hidx++)
      Yavg1half += Yvals[hidx];
    Yavg1half /= (double)halfwayPt;
    statsout_->Printf("\tFirst half <Y> = %g\n", Yavg1half);
    double Yavg2half = 0;
    for (unsigned int hidx = halfwayPt; hidx < Yvals.size(); hidx++)
      Yavg2half += Yvals[hidx];
    Yavg2half /= (double)(Yvals.size() - halfwayPt);
    statsout_->Printf("\tLast half <Y> = %g\n", Yavg2half);*/

    // Set initial guesses for parameters.
    CurveFit::Darray Params(3);
    Params[0] = YinitialAvg;
    if (Params[0] < 0) Params[0] = -Params[0];
    Params[1] = 0.1; // TODO absolute slope?
    //Params[2] = Yavg2half;
    Params[2] = YfinalAvg;

    for (CurveFit::Darray::const_iterator ip = Params.begin(); ip != Params.end(); ++ip) {
      statsout_->Printf("\tInitial Param A%li = %g\n", ip - Params.begin(), *ip);
    }

    // Perform curve fitting
    CurveFit fit;
    int info = fit.LevenbergMarquardt( fxn, Xvals, Yvals, Params, tolerance_, maxIt_ );
    mprintf("\t%s\n", fit.Message(info));
    if (info == 0) {
      mprintf("Warning: Curve fit failed for '%s'.\n"
              "Warning: %s\n", DS.legend(), fit.ErrorMessage());
      BlankResult(oidx, DS.legend());
      continue;
    }
    for (CurveFit::Darray::const_iterator ip = Params.begin(); ip != Params.end(); ++ip) {
      statsout_->Printf("\tFinal Param A%li = %g\n", ip - Params.begin(), *ip);
    }

    // -------------------------------------------
    // Params[0] = A0 = Initial value. 
    // Params[1] = A1 = Relaxation constant
    // Params[2] = A2 = Final value at long time.
    // Determine the absolute difference of the long-time estimated value
    // from the average value of the last half of the data.
    //double ValA = Yavg2half - Params[2];
    double ValA = YfinalAvg - Params[2];
    if (ValA < 0) ValA = -ValA;
    statsout_->Printf("\tValA = %g\n", ValA);

    // Create output curve
    double avgr = 0.0;
    //CpptrajFile residualout;
    //residualout.OpenWrite( "residual." + integerToString(it-inputSets_.begin()) + ".dat");
    DataSet_Mesh& OUT = static_cast<DataSet_Mesh&>( *(*ot) );
    for (unsigned int i = 0; i != Xvals.size(); i++) {
      OUT.AddXY( Xvals[i] + offset, fit.FinalY()[i] );
      double residual = fit.FinalY()[i] - Yvals[i];
      //residualout.Printf("%g\t%g\n", Xvals[i], residual);
      avgr += residual;
    }
    //residualout.CloseFile();
    avgr /= (double)Xvals.size();
    statsout_->Printf("\tAvgResidual= %g\n", avgr);

    // Calculate where slope reaches slopeCut_
/*
    DataSet_1D::Darray slopeX, slopeY;
    OUT.FiniteDifference(DataSet_1D::FORWARD, slopeX, slopeY);
    double finalx=-1, finaly=0;
    double finalslope = 0;
    bool slopeCutSatisfied = false;
    for (unsigned int sidx = 0; sidx != slopeY.size(); ++sidx) {
      //finalslope = slopeY[sidx];
      //mprintf("DEBUG: slope %g %g\n", slopeX[sidx], slopeY[sidx]);
      double absSlope = slopeY[sidx];
      if (absSlope < 0) absSlope = -absSlope;
      finalslope = absSlope;
      if (absSlope < slopeCut_) {
        finalx = slopeX[sidx];
        finaly = slopeY[sidx];
        slopeCutSatisfied = true;
        break;
      }
    }
    statsout_->Printf("\tFinal slope: %g\n", slopeY.back());
    if (slopeCutSatisfied)
      statsout_->Printf("\tSlope cutoff satisfied at %g %g\n", finalx, finaly);
    else
      statsout_->Printf("\tSlope cutoff not satisfied.\n");
*/

    // Statistics
    double corr_coeff, ChiSq, TheilU, rms_percent_error;
    err = fit.Statistics( Yvals, corr_coeff, ChiSq, TheilU, rms_percent_error);
    if (err != 0) mprintf("Warning: %s\n", fit.Message(err));
    statsout_->Printf("\tCorrelation coefficient: %g\n"
                      "\tChi squared: %g\n"
                      "\tUncertainty coefficient: %g\n"
                      "\tRMS percent error: %g\n",
                      corr_coeff, ChiSq, TheilU, rms_percent_error);

    data_[A0]->Add(oidx, &Params[0]);
    data_[A1]->Add(oidx, &Params[0] + 1);
    data_[A2]->Add(oidx, &Params[0] + 2);
    double onea1 = 1.0 / Params[1];
    data_[ONEA1]->Add(oidx, &onea1);
    //data_[FVAL]->Add(oidx, &Fval);
    data_[CORR]->Add(oidx, &corr_coeff);
    //if (ValA < 0.0) ValA = -ValA;
    data_[VALA]->Add(oidx, &ValA);
    //data_[LCHISQ]->Add(oidx, &linChiSq);
    data_[CHISQ]->Add(oidx, &ChiSq);
    //data_[PLTIME]->Add(oidx, &finalx);
    //data_[FSLOPE]->Add(oidx, &finalslope);
    data_[NAME]->Add(oidx, (*it)->legend());

    // Determine if criteria met.

    // This version checks chi squared, ValA, and final slope.
    statsout_->Printf("----- Criteria -----\n");
    // Chi^2
    bool chiSqSatisfied = false;
    if (ChiSq < chisqCut_)
      chiSqSatisfied = true;
    // ValA
    bool valaSatisfied = false;
    if (ValA < valaCut_)
      valaSatisfied = true;
    // Slope
    double pltime = -1;
    double finalslope = 0;
    bool slopeSatisfied = CheckSlope(OUT, pltime, finalslope);
    data_[FSLOPE]->Add(oidx, &finalslope);
    data_[PLTIME]->Add(oidx, &pltime);
    // Final result
    char finalResult[4];
    finalResult[3] = '\0';
    if (chiSqSatisfied && valaSatisfied && slopeSatisfied) {
      finalResult[0] = 'y';
      finalResult[1] = 'e';
      finalResult[2] = 's';
    } else {
      if (chiSqSatisfied)
        finalResult[0] = 'c';
      else
        finalResult[0] = '.';
      if (valaSatisfied)
        finalResult[1] = 'v';
      else
        finalResult[1] = '.';
      if (slopeSatisfied)
        finalResult[2] = 's';
      else
        finalResult[2] = '.';
    }
    data_[RESULT]->Add(oidx, finalResult);

/*
    // This version checks the exp fit against a linear fit.
    double finalslope = 0;
    double pltime = -1;
    bool satisfied = false;
    statsout_->Printf("----- Criteria -----\n");
    if ( (ChiSq < chisqCut_) || (linChiSq < chisqCut_) ) {
      statsout_->Printf("\tChi-squared cutoff satisfied.\n");
      // Chi-squared cutoff satisfied for one or both.
      // Which model fit better? Lower chi-squared wins. Tie goes to exp.
      if (linChiSq < ChiSq) {
        statsout_->Printf("\tLinear model has lower chi-squared value.\n");
        finalslope = slope;
        if (finalslope < 0)
          finalslope = -finalslope;
        if (finalslope < slopeCut_) {
          statsout_->Printf("\tLinear model slope satisfies the slope cut.\n");
          data_[RESULT]->Add(oidx, "lin");
          pltime = OUT.Xcrd(0);
          satisfied = true;
        } else
          statsout_->Printf("\tLinear model does not satistfy the slope cut.\n");
      } else {
        statsout_->Printf("\tExponential model has lower chi-squared value.\n");
        DataSet_1D::Darray slopeX, slopeY;
        OUT.FiniteDifference(DataSet_1D::FORWARD, slopeX, slopeY);
        double finalx=-1, finaly=0;
        bool slopeCutSatisfied = false;
        for (unsigned int sidx = 0; sidx != slopeY.size(); ++sidx) {
          //finalslope = slopeY[sidx];
          //mprintf("DEBUG: slope %g %g\n", slopeX[sidx], slopeY[sidx]);
          double absSlope = slopeY[sidx];
          if (absSlope < 0) absSlope = -absSlope;
          finalslope = absSlope;
          if (absSlope < slopeCut_) {
            finalx = slopeX[sidx];
            finaly = slopeY[sidx];
            slopeCutSatisfied = true;
            break;
          }
        }
        pltime = finalx;
        statsout_->Printf("\tFinal slope: %g\n", slopeY.back());
        if (slopeCutSatisfied) {
          statsout_->Printf("\tSlope cutoff satisfied at %g %g\n", finalx, finaly);
          data_[RESULT]->Add(oidx, "exp");
          satisfied = true;
        } else
          statsout_->Printf("\tSlope cutoff not satisfied.\n");
      }
    } else
      statsout_->Printf("\tChi-squared cutoff not satisfied.\n"); 
    data_[FSLOPE]->Add(oidx, &finalslope);
    data_[PLTIME]->Add(oidx, &pltime);
    if (!satisfied)
      data_[RESULT]->Add(oidx,"no");
*/
   

/*
    bool longAvgCutSatisfied, chiCutSatisfied;
    if (ValA < valaCut_)
      longAvgCutSatisfied = true;
    else {
      longAvgCutSatisfied = false;
      mprintf("\tLong-time average cut not satisfied.\n");
    }
    if (ChiSq < chisqCut_)
      chiCutSatisfied = true;
    else {
      chiCutSatisfied = false;
      mprintf("\tNon-linear fit chi-squared not satisfied.\n");
    }
    if ( longAvgCutSatisfied && chiCutSatisfied && slopeCutSatisfied)
      data_[RESULT]->Add(oidx, "yes");
    else
      data_[RESULT]->Add(oidx, "no");*/

    statsout_->Printf("\n");
  } // END loop over input data sets
  return Analysis::OK;
}
