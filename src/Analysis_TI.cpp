#include <cmath> // sqrt
#include "Analysis_TI.h"
#include "CpptrajStdio.h"
#include "DataSet_Mesh.h"
#include "StringRoutines.h" // integerToString
#include "Random.h"

/// CONSTRUCTOR
Analysis_TI::Analysis_TI() :
  nskip_(0),
  dAout_(0),
  dA_SD_(0),
  curveout_(0),
  masterDSL_(0),
  bootstrap_fac_(0.0),
  mode_(GAUSSIAN_QUAD),
  avgType_(AVG),
  debug_(0),
  n_bootstrap_pts_(0),
  n_bootstrap_samples_(0),
  bootstrap_seed_(-1),
  avg_increment_(0),
  avg_max_(0),
  avg_skip_(0)
{}

// Analysis_TI::Help()
void Analysis_TI::Help() const {
  mprintf("\t<dset0> [<dset1> ...] {nq <n quad pts> | xvals <x values>}\n"
          "\t[name <set name>] [out <file>] [curveout <ti curve file>]\n"
          "\t[nskip <#s to skip>]\n"
          "\t[avgincrement <#> [avgmax <#>] [avgskip <#>]]\n"
          "\t[bs_samples <samples> [bs_points <points>] [bs_seed <#>]\n"
          "\t [bs_fac <factor>]]\n"
          "  Calculate free energy using DV/DL energies from thermodynamic integration.\n"
          "  The results of integration of the DV/DL curve will be written to <file>,\n"
          "  while the curves themselves will be written to <ti curve file>.\n"
          "  Use 'nq' to specify number of Gaussian quadrature points; otherwise\n"
          "  the lambda values should be specified by 'xvals', where <x values>\n"
          "  is a comma-separated list.\n" 
          "  If 'nskip' is specified (where <# to skip> may be a comma-separated\n"
          "  list of numbers) the average DV/DL and final free energy will be\n"
          "  calculated skipping over the specified number(s) of points.\n"
          "  If 'avgincrement' is specified the average DV/DL and final free\n"
          "  energy will be calculated from incremental averages.\n"
          "  If 'bs_samples' is specified the average DV/DL will be calcuated\n"
          "  in bootstrap fashion using <points> for <samples> samples.\n"
          "  If <points> is not specified use the data size times <factor> points.\n");
}

// Analysis_TI::Setup()
Analysis::RetType Analysis_TI::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  debug_ = debugIn;
  int nq = analyzeArgs.getKeyInt("nq", 0);
  ArgList nskipArg(analyzeArgs.GetStringKey("nskip"), ","); // Comma-separated
  avg_increment_ = analyzeArgs.getKeyInt("avgincrement", -1);
  avg_max_ = analyzeArgs.getKeyInt("avgmax", -1);
  avg_skip_ = analyzeArgs.getKeyInt("avgskip", 0);
  n_bootstrap_pts_ = analyzeArgs.getKeyInt("bs_pts", -1);
  n_bootstrap_samples_ = analyzeArgs.getKeyInt("bs_samples", 0);
  bootstrap_seed_ = analyzeArgs.getKeyInt("bs_seed", -1);
  bootstrap_fac_ = analyzeArgs.getKeyDouble("bs_fac", 0.75);
  if (!nskipArg.empty()) {
    avgType_ = SKIP;
    // Specified numbers of points to skip
    nskip_.clear();
    for (int i = 0; i != nskipArg.Nargs(); i++) {
      nskip_.push_back( nskipArg.getNextInteger(0) );
      if (nskip_.back() < 0) nskip_.back() = 0;
    }
  } else if (avg_increment_ > 0)
    avgType_ = INCREMENT;
  else if (n_bootstrap_samples_ > 0)
    avgType_ = BOOTSTRAP;
  else
    avgType_ = AVG;
  masterDSL_ = setup.DslPtr();
  // Get lambda values
  ArgList xArgs(analyzeArgs.GetStringKey("xvals"), ","); // Also comma-separated
  if (!xArgs.empty()) {
    xval_.clear();
    for (int i = 0; i != xArgs.Nargs(); i++)
      xval_.push_back( xArgs.getNextDouble(0.0) );
  }
  std::string setname = analyzeArgs.GetStringKey("name");
  DataFile* outfile = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  curveout_ = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("curveout"), analyzeArgs);
  // Select datasets from remaining args
  if (input_dsets_.AddSetsFromArgs( analyzeArgs.RemainingArgs(), setup.DSL() )) {
    mprinterr("Error: Could not add data sets.\n");
    return Analysis::ERR;
  }
  if (input_dsets_.empty()) {
    mprinterr("Error: No input data sets.\n");
    return Analysis::ERR;
  }
  if (SetQuadAndWeights(nq)) return Analysis::ERR;
  // Determine integration mode
  if (nq > 0)
    mode_ = GAUSSIAN_QUAD;
  else
    mode_ = TRAPEZOID;
  // Check that # abscissas matches # data sets
  if (xval_.size() != input_dsets_.size()) {
     mprinterr("Error: Expected %zu data sets for integration, got %zu\n",
               input_dsets_.size(), xval_.size());
    return Analysis::ERR;
  }
  // Set up output data sets
  DataSet::DataType dtype = DataSet::DOUBLE;
  if (avgType_ == SKIP || avgType_ == INCREMENT)
    dtype = DataSet::XYMESH;
  dAout_ = setup.DSL().AddSet(dtype, setname, "TI");
  if (dAout_ == 0) return Analysis::ERR;
  if (outfile != 0) outfile->AddDataSet( dAout_ );
  MetaData md(dAout_->Meta().Name(), "TIcurve");
  if (avgType_ == AVG) {
    // Single curve
    curve_.push_back( setup.DSL().AddSet(DataSet::XYMESH, md) );
    if (curve_.back() == 0) return Analysis::ERR;
    if (curveout_ != 0) curveout_->AddDataSet( curve_.back() );
    if (outfile != 0) outfile->ProcessArgs("noxcol");
  } else if (avgType_ == SKIP) {
    // As many curves as skip values
    for (Iarray::const_iterator it = nskip_.begin(); it != nskip_.end(); ++it) {
      md.SetIdx( *it );
      DataSet* ds = setup.DSL().AddSet(DataSet::XYMESH, md);
      if (ds == 0) return Analysis::ERR;
      ds->SetLegend( md.Name() + "_Skip" + integerToString(*it) );
      if (curveout_ != 0) curveout_->AddDataSet( ds );
      curve_.push_back( ds );
    }
  } else if (avgType_ == BOOTSTRAP) {
    // As many curves as resamples
    for (int nsample = 0; nsample != n_bootstrap_samples_; nsample++) {
      md.SetIdx(nsample);
      DataSet* ds = setup.DSL().AddSet(DataSet::XYMESH, md);
      if (ds == 0) return Analysis::ERR;
      ds->SetLegend( md.Name() + "_Sample" + integerToString(nsample) );
      if (curveout_ != 0) curveout_->AddDataSet( ds );
      curve_.push_back( ds );
    }
    // Standard devation of avg free energy over samples
    dA_SD_ = setup.DSL().AddSet(DataSet::DOUBLE, MetaData(md.Name(), "SD"));
    if (dA_SD_ == 0) return Analysis::ERR;
    if (outfile != 0) {
      outfile->AddDataSet( dA_SD_ );
      outfile->ProcessArgs("noxcol");
    }
  }
  // NOTE: INCREMENT is set up once data set size is known 

  mprintf("    TI: Calculating TI");
  if (mode_ == GAUSSIAN_QUAD) {
    mprintf(" using Gaussian quadrature with %zu points.\n", xval_.size());
    mprintf("\t%6s %8s %8s %s\n", "Point", "Abscissa", "Weight", "SetName");
    for (unsigned int i = 0; i != xval_.size(); i++)
      mprintf("\t%6i %8.5f %8.5f %s\n", i, xval_[i], wgt_[i], input_dsets_[i]->legend());
  } else {
    mprintf(" using the trapezoid rule.\n");
    mprintf("\t%6s %8s %s\n", "Point", "Abscissa", "SetName");
    for (unsigned int i = 0; i != xval_.size(); i++)
      mprintf("\t%6i %8.5f %s\n", i, xval_[i], input_dsets_[i]->legend());
  }
  mprintf("\tResult(s) of integration(s) saved in set '%s'\n", dAout_->legend());
  if (avgType_ == AVG)
    mprintf("\tUsing all data points in <DV/DL> calc.\n");
  else if (avgType_ == SKIP) {
    mprintf("\tSkipping first");
    for (Iarray::const_iterator it = nskip_.begin(); it != nskip_.end(); ++it)
      mprintf(" %i", *it);
    mprintf(" data points for <DV/DL> calc.\n");
  } else if (avgType_ == INCREMENT) {
    mprintf("\tCalculating <DV/DL> starting from point %i, increment by %i.",
            avg_skip_, avg_increment_);
    if (avg_max_ != -1)
      mprintf(" Max %i points.", avg_max_);
    mprintf("\n");
  } else if (avgType_ == BOOTSTRAP) {
    mprintf("\tStandard devation of result stored in set '%s'\n", dA_SD_->legend());
    mprintf("\tCalculating <DV/DL> from %i bootstrap resamples.\n", n_bootstrap_samples_);
    if (n_bootstrap_pts_ > 0)
      mprintf("\tBootstrap resample size is %i data points.\n", n_bootstrap_pts_);
    else
      mprintf("\tWill use bootstrap resample size of %g%% of total points.\n",
              bootstrap_fac_*100.0);
    if (bootstrap_seed_ != -1)
      mprintf("\tBoostrap base seed is %i\n", bootstrap_seed_);
  }
  mprintf("\tTI curve(s) saved in set(s)");
  if (avgType_ != INCREMENT)
    for (DSarray::const_iterator ds = curve_.begin(); ds != curve_.end(); ++ds)
      mprintf(" '%s'", (*ds)->legend());
  else
    mprintf(" named '%s'", md.PrintName().c_str());
  mprintf("\n");
  if (outfile != 0) mprintf("\tResults written to '%s'\n", outfile->DataFilename().full());
  if (curveout_!= 0) mprintf("\tTI curve(s) written to '%s'\n", curveout_->DataFilename().full());

  return Analysis::OK;
}

// Analysis_TI::SetQuadAndWeights()
int Analysis_TI::SetQuadAndWeights(int nq) {
  if (nq < 1) return 0;
  xval_.resize(nq);
  wgt_.resize(nq);
  switch (nq) {
    case 1:
      xval_[0] = 0.5; wgt_[0] = 1.0;
      break;
    case 2:
      xval_[0] = 0.21132; wgt_[0] = wgt_[1] = 0.5;
      xval_[1] = 0.78867;
      break;
    case 3:
      xval_[0] = 0.1127;  wgt_[0] = wgt_[2] = 0.27777;
      xval_[1] = 0.5;     wgt_[1] = 0.44444;
      xval_[2] = 0.88729;
      break;
    case 5: 
      xval_[0] = 0.04691; wgt_[0] = wgt_[4] = 0.11846;
      xval_[1] = 0.23076; wgt_[1] = wgt_[3] = 0.23931;
      xval_[2] = 0.5;     wgt_[2] = 0.28444;
      xval_[3] = 0.76923;
      xval_[4] = 0.95308;
      break; 
    case 7:
      xval_[0] = 0.02544; wgt_[0] = wgt_[6] = 0.06474;
      xval_[1] = 0.12923; wgt_[1] = wgt_[5] = 0.13985;
      xval_[2] = 0.29707; wgt_[2] = wgt_[4] = 0.19091;
      xval_[3] = 0.5;     wgt_[3] = 0.20897;
      xval_[4] = 0.70292;
      xval_[5] = 0.87076;
      xval_[6] = 0.97455;
      break;
    case 9:
      xval_[0] = 0.01592; wgt_[0] = wgt_[8] = 0.04064;
      xval_[1] = 0.08198; wgt_[1] = wgt_[7] = 0.09032;
      xval_[2] = 0.19331; wgt_[2] = wgt_[6] = 0.13031;
      xval_[3] = 0.33787; wgt_[3] = wgt_[5] = 0.15617;
      xval_[4] = 0.5;     wgt_[4] = 0.16512;
      xval_[5] = 0.66213;
      xval_[6] = 0.80669;
      xval_[7] = 0.91802;
      xval_[8] = 0.98408;
      break;
    case 12:
      xval_[0] = 0.00922; wgt_[0] = wgt_[11] = 0.02359;
      xval_[1] = 0.04794; wgt_[1] = wgt_[10] = 0.05347;
      xval_[2] = 0.11505; wgt_[2] = wgt_[9]  = 0.08004;
      xval_[3] = 0.20634; wgt_[3] = wgt_[8]  = 0.10158;
      xval_[4] = 0.31608; wgt_[4] = wgt_[7]  = 0.11675;
      xval_[5] = 0.43738; wgt_[5] = wgt_[6]  = 0.12457;
      xval_[6] = 0.56262;
      xval_[7] = 0.68392;
      xval_[8] = 0.79366;
      xval_[9] = 0.88495;
      xval_[10] = 0.95206;
      xval_[11] = 0.99078;
      break;
    default:
      mprinterr("Error: Unsupported quadrature: %i\n", nq);
      return 1;
  }
  return 0;
}

/** Integrate each curve in the curve_ array using trapezoid method. */
void Analysis_TI::Integrate_Trapezoid(Darray& sum) const {
  // Integrate each curve if not doing quadrature
  for (unsigned int j = 0; j != curve_.size(); j++) {
    DataSet_Mesh const& CR = static_cast<DataSet_Mesh const&>( *(curve_[j]) );
    sum[j] = CR.Integrate_Trapezoid();
  }
}

static inline int CheckSet(DataSet_1D const& ds) {
  if (ds.Size() < 1) {
    mprinterr("Error: Set '%s' is empty.\n", ds.legend());
    return Analysis::ERR;
  }
  mprintf("\t%s (%zu points).\n", ds.legend(), ds.Size());
  return 0;
}

// Analysis_TI::Calc_Bootstrap()
int Analysis_TI::Calc_Bootstrap() {
  // sum: Hold the results of integration for each curve (bootstrap resample)
  Darray sum(n_bootstrap_samples_, 0.0);
  Random_Number RN;
  RN.rn_set( bootstrap_seed_ );
  
  // Loop over input data sets.
  for (unsigned int idx = 0; idx != input_dsets_.size(); idx++) {
    DataSet_1D const& ds = static_cast<DataSet_1D const&>( *(input_dsets_[idx]) );
    if (CheckSet(ds)) return 1;
    // True if original point already chosen this round
    std::vector<bool> chosen;
    // Hold averages for each resample
    Darray Avgs(n_bootstrap_samples_, 0.0);
    // Hold average of all resample averages
    double Mean = 0.0;
    if (n_bootstrap_pts_ < 1) {
      n_bootstrap_pts_ = (int)(bootstrap_fac_ * (double)ds.Size());
      mprintf("\tUsing %i bootstrap points.\n", n_bootstrap_pts_);
    }
    if (n_bootstrap_pts_ < 1) {
      mprinterr("Error: Not enough bootstrap points.\n");
      return 1;
    }
    if (n_bootstrap_pts_ >= (int)ds.Size()) {
      mprinterr("Error: Bootstrap sample size (%i) must be less than data set size (%zu)\n",
                n_bootstrap_pts_, ds.Size());
      return 1;
    }
    if (n_bootstrap_samples_ > (int)ds.Size() - n_bootstrap_pts_)
    mprintf("Warning: Bootstrap # resamples (%i) > data size (%zu) - sample size (%i)\n",
            n_bootstrap_samples_, ds.Size(), n_bootstrap_pts_);
    // Loop over resamples
    double d_ndata = (double)ds.Size();
    for (int nsample = 0; nsample != n_bootstrap_samples_; nsample++)
    {
      chosen.assign(ds.Size(), false);
      // Sum up sample_size randomly chosen points
      for (int npoint = 0; npoint != n_bootstrap_pts_; npoint++)
      {
        bool pointOK = false;
        unsigned int pt = 0;
        while (!pointOK)
        {
          pt = (unsigned int)(RN.rn_gen() * d_ndata);
          pointOK = (chosen[pt] == false);
        }
        chosen[pt] = true;
        Avgs[nsample] += ds.Dval( pt );
      }
      Avgs[nsample] /= (double)n_bootstrap_pts_;
      Mean += Avgs[nsample];
    }
    // Mean of all resamples
    Mean /= (double)n_bootstrap_samples_;
    mprintf("\tLambda %g: Avg of resamples= %g\n", xval_[idx], Mean);
    // Store average DV/DL for each resample 
    for (unsigned int j = 0; j != Avgs.size(); j++) {
      if (debug_ > 0)
        mprintf("\t%s Resample %u <DV/DL>= %g\n", ds.legend(), j, Avgs[j]);
      DataSet_Mesh& CR = static_cast<DataSet_Mesh&>( *(curve_[j]) );
      CR.AddXY(xval_[idx], Avgs[j]);
      if (mode_ == GAUSSIAN_QUAD)
        sum[j] += (wgt_[idx] * Avgs[j]);
    }
  } // END loop over input data sets
  if (mode_ == TRAPEZOID) Integrate_Trapezoid(sum);
  // Calculate average DA from individual TI integration values.
  double DA_avg = 0.0;
  double DA_sd = 0.0;
  for (unsigned int j = 0; j != sum.size(); j++) {
    mprintf("\tResample %u DA= %g\n", j, sum[j]);
    DA_avg += sum[j];
    DA_sd += (sum[j] * sum[j]);
  }
  DA_avg /= (double)sum.size();
  DA_sd /= (double)sum.size();
  DA_sd -= (DA_avg * DA_avg);
  if (DA_sd > 0.0)
    DA_sd = sqrt(DA_sd);
  else
    DA_sd = 0.0;
  mprintf("\tBootstrap avg +/- SD = %g +/- %g\n", DA_avg, DA_sd);
  dAout_->ModifyDim(Dimension::X).SetLabel("Bootstrap");
  dAout_->Add(0, &DA_avg);
  dA_SD_->ModifyDim(Dimension::X).SetLabel("Bootstrap");
  dA_SD_->Add(0, &DA_sd);

  return 0;
}

// Analysis_TI::Calc_Nskip()
int Analysis_TI::Calc_Nskip() {
  // sum: Hold the results of integration for each curve (skip value)
  Darray sum(nskip_.size(), 0.0);
  // lastSkipPoint: Points after which averages can be recorded
  Iarray lastSkipPoint;
  for (Iarray::const_iterator it = nskip_.begin(); it != nskip_.end(); ++it)
    lastSkipPoint.push_back( *it - 1 );
  // Loop over input data sets. 
  for (unsigned int idx = 0; idx != input_dsets_.size(); idx++) {
    DataSet_1D const& ds = static_cast<DataSet_1D const&>( *(input_dsets_[idx]) );
    if (CheckSet(ds)) return 1; 
    // Determine if skip values are valid for this set.
    Darray Npoints; // Number of points after skipping
    for (Iarray::const_iterator it = nskip_.begin(); it != nskip_.end(); ++it) {
      int np = (int)ds.Size() - *it;
      if (np < 1) {
        mprinterr("Error: Skipped too many points (set '%s' size is %zu)\n",ds.legend(),ds.Size());
        return 1;
      }
      Npoints.push_back((double)np);
    }
    // Calculate averages for each value of skip
    Darray avg(nskip_.size(), 0.0);
    for (int i = 0; i != (int)ds.Size(); i++) {
      for (unsigned int j = 0; j != nskip_.size(); j++)
        if (i > lastSkipPoint[j])
          avg[j] += ds.Dval( i );
    }
    // Store average DV/DL for each value of skip
    for (unsigned int j = 0; j != nskip_.size(); j++) {
      avg[j] /= Npoints[j];
      if (debug_ > 0)
        mprintf("\t%s Skip= %i <DV/DL>= %g\n", ds.legend(), nskip_[j], avg[j]);
      DataSet_Mesh& CR = static_cast<DataSet_Mesh&>( *(curve_[j]) );
      CR.AddXY(xval_[idx], avg[j]);
      if (mode_ == GAUSSIAN_QUAD)
        sum[j] += (wgt_[idx] * avg[j]);
    }
  } // END loop over input data sets
  if (mode_ == TRAPEZOID) Integrate_Trapezoid(sum);
  // Store final TI integration values.
  DataSet_Mesh& DA = static_cast<DataSet_Mesh&>( *dAout_ );
  DA.ModifyDim(Dimension::X).SetLabel("PtsSkipped");
  for (unsigned int j = 0; j != nskip_.size(); j++)
    DA.AddXY(nskip_[j], sum[j]);

  return 0;
}

// Analysis_TI::Calc_Increment()
int Analysis_TI::Calc_Increment() {
  // sum: Hold the results of integration for each curve (increment)
  Darray sum;
  // points: Hold point values at which each avg is being calculated
  Iarray points;
  // Loop over input data sets. 
  for (unsigned int idx = 0; idx != input_dsets_.size(); idx++) {
    DataSet_1D const& ds = static_cast<DataSet_1D const&>( *(input_dsets_[idx]) );
    if (CheckSet(ds)) return 1; 
    // Determine max pts if not given
    int maxpts = avg_max_;
    if (maxpts == -1)
      maxpts = (int)ds.Size();
    else if (maxpts > (int)ds.Size()) {
      mprintf("Warning: 'avgmax' (%i) > data size (%zu); setting to %zu\n",
              maxpts, ds.Size(), ds.Size());
      maxpts = (int)ds.Size();
    }
    if (avg_skip_ >= maxpts) {
      mprinterr("Error: 'avgskip' (%i) > max (%i).\n", avg_skip_, maxpts);
      return 1;
    }
    // Calculate averages for each increment
    Darray avg;
    Iarray increments;
    int count = 0;
    int endpt = maxpts -1;
    double currentSum = 0.0;
    if (debug_ > 0) mprintf("DEBUG: Lambda %g\n", xval_[idx]);
    for (int pt = avg_skip_; pt != maxpts; pt++)
    {
      currentSum += ds.Dval(pt);
      count++;
      if (count == avg_increment_ || pt == endpt) {
        avg.push_back( currentSum / ((double)(pt - avg_skip_ + 1)) );
        increments.push_back(pt+1);
        if (debug_ > 0)
          mprintf("DEBUG:\t\tAvg from %i to %i: %g\n", avg_skip_+1, pt+1, avg.back());
        count = 0;
      }
    }
    if (sum.empty()) {
      sum.resize(avg.size());
      points = increments;
    } else if (sum.size() != avg.size()) {
      mprinterr("Error: Different # of increments for set '%s'; got %zu, expected %zu.\n",
                ds.legend(), avg.size(), sum.size());
      return 1;
    }
    // Create increment curve data sets
    if (curve_.empty()) {
      MetaData md(dAout_->Meta().Name(), "TIcurve");
      for (unsigned int j = 0; j != avg.size(); j++) {
        md.SetIdx( increments[j] );
        DataSet* ds = masterDSL_->AddSet(DataSet::XYMESH, md);
        if (ds == 0) return Analysis::ERR;
        ds->SetLegend( md.Name() + "_Skip" + integerToString(increments[j]) );
        if (curveout_ != 0) curveout_->AddDataSet( ds );
        curve_.push_back( ds );
      }
    }
    for (unsigned int j = 0; j != avg.size(); j++) {
      DataSet_Mesh& CR = static_cast<DataSet_Mesh&>( *(curve_[j]) );
      CR.AddXY(xval_[idx], avg[j]);
      if (mode_ == GAUSSIAN_QUAD)
        sum[j] += (wgt_[idx] * avg[j]);
    }
  } // END loop over data sets
  if (mode_ == TRAPEZOID) Integrate_Trapezoid(sum);
  // Store final integration values
  DataSet_Mesh& DA = static_cast<DataSet_Mesh&>( *dAout_ );
  DA.ModifyDim(Dimension::X).SetLabel("Point");
  for (unsigned int j = 0; j != points.size(); j++)
    DA.AddXY(points[j], sum[j]);

  return 0;
}

// Analysis_TI::Calc_Avg()
int Analysis_TI::Calc_Avg() {
  // sum: Hold the results of integration single curve
  Darray sum(1, 0.0);
   // Loop over input data sets. 
  for (unsigned int idx = 0; idx != input_dsets_.size(); idx++) {
    DataSet_1D const& ds = static_cast<DataSet_1D const&>( *(input_dsets_[idx]) );
    if (CheckSet(ds)) return 1; 
    // Calculate average for this set
    double avg_dvdl = ds.Avg();
    ((DataSet_Mesh*)curve_[0])->AddXY( xval_[idx], avg_dvdl );
    // Gaussian Quad.
    if (mode_ == GAUSSIAN_QUAD)
      sum[0] += (wgt_[idx] * avg_dvdl);
  }
  if (mode_ == TRAPEZOID) Integrate_Trapezoid(sum);
  // Store final integration values
  dAout_->ModifyDim(Dimension::X).SetLabel("TI");
  dAout_->Add(0, &(sum[0]));
  return 0;
}

// Analysis_TI::Analyze()
Analysis::RetType Analysis_TI::Analyze() {
  int err = 0;
  if (avgType_ == SKIP)
    err = Calc_Nskip();
  else if (avgType_ == AVG)
    err = Calc_Avg();
  else if (avgType_ == INCREMENT)
    err = Calc_Increment();
  else if (avgType_ == BOOTSTRAP)
    err = Calc_Bootstrap();
  if (err != 0) return Analysis::ERR;

  return Analysis::OK;
}
