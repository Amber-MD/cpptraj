#include "Analysis_TI.h"
#include "CpptrajStdio.h"
#include "DataSet_Mesh.h"
#include "StringRoutines.h" // integerToString

void Analysis_TI::Help() const {
  mprintf("\t<dset0> [<dset1> ...] [nq <n quad pts>] [nskip <# to skip>]\n"
          "\t[name <set name>] [out <file>]\n" 
          "  Calculate free energy from Amber TI output.\n");
}

int Analysis_TI::SetQuadAndWeights(int nq) {
  quad_.clear();
  wgt_.clear();
  if (nq < 1) return 1;
  quad_.resize(nq);
  wgt_.resize(nq);
  switch (nq) {
    case 1:
      quad_[0] = 0.5; wgt_[0] = 1.0;
      break;
    case 2:
      quad_[0] = 0.21132; wgt_[0] = wgt_[1] = 0.5;
      quad_[1] = 0.78867;
      break;
    case 3:
      quad_[0] = 0.1127;  wgt_[0] = wgt_[2] = 0.27777;
      quad_[1] = 0.5;     wgt_[1] = 0.44444;
      quad_[2] = 0.88729;
      break;
    case 5: 
      quad_[0] = 0.04691; wgt_[0] = wgt_[4] = 0.11846;
      quad_[1] = 0.23076; wgt_[1] = wgt_[3] = 0.23931;
      quad_[2] = 0.5;     wgt_[2] = 0.28444;
      quad_[3] = 0.76923;
      quad_[4] = 0.95308;
      break; 
    case 7:
      quad_[0] = 0.02544; wgt_[0] = wgt_[6] = 0.06474;
      quad_[1] = 0.12923; wgt_[1] = wgt_[5] = 0.13985;
      quad_[2] = 0.29707; wgt_[2] = wgt_[4] = 0.19091;
      quad_[3] = 0.5;     wgt_[3] = 0.20897;
      quad_[4] = 0.70292;
      quad_[5] = 0.87076;
      quad_[6] = 0.97455;
      break;
    case 9:
      quad_[0] = 0.01592; wgt_[0] = wgt_[8] = 0.04064;
      quad_[1] = 0.08198; wgt_[1] = wgt_[7] = 0.09032;
      quad_[2] = 0.19331; wgt_[2] = wgt_[6] = 0.13031;
      quad_[3] = 0.33787; wgt_[3] = wgt_[5] = 0.15617;
      quad_[4] = 0.5;     wgt_[4] = 0.16512;
      quad_[5] = 0.66213;
      quad_[6] = 0.80669;
      quad_[7] = 0.91802;
      quad_[8] = 0.98408;
      break;
    case 12:
      quad_[0] = 0.00922; wgt_[0] = wgt_[11] = 0.02359;
      quad_[1] = 0.04794; wgt_[1] = wgt_[10] = 0.05347;
      quad_[2] = 0.11505; wgt_[2] = wgt_[9]  = 0.08004;
      quad_[3] = 0.20634; wgt_[3] = wgt_[8]  = 0.10158;
      quad_[4] = 0.31608; wgt_[4] = wgt_[7]  = 0.11675;
      quad_[5] = 0.43738; wgt_[5] = wgt_[6]  = 0.12457;
      quad_[6] = 0.56262;
      quad_[7] = 0.68392;
      quad_[8] = 0.79366;
      quad_[9] = 0.88495;
      quad_[10] = 0.95206;
      quad_[11] = 0.99078;
      break;
    default:
      mprinterr("Error: Unsupported quadrature: %i\n", nq);
      return 1;
  }
  return 0;
}

// Analysis_TI::Setup()
Analysis::RetType Analysis_TI::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  int nq = analyzeArgs.getKeyInt("nq", 0);
  ArgList nskipArg(analyzeArgs.GetStringKey("nskip"), ","); // Comma-separated
  if (nskipArg.empty())
    nskip_.resize(1, 0);
  else {
    nskip_.clear();
    for (int i = 0; i != nskipArg.Nargs(); i++) {
      nskip_.push_back( nskipArg.getNextInteger(0) );
      if (nskip_.back() < 0) nskip_.back() = 0;
    }
  }
  std::string setname = analyzeArgs.GetStringKey("name");
  DataFile* outfile = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  DataFile* curveout = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("curveout"), analyzeArgs);
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
  if (quad_.size() != input_dsets_.size()) {
    mprinterr("Error: Expected %zu data sets based on nq, got %zu\n",
              quad_.size(), input_dsets_.size());
    return Analysis::ERR;
  }
  dAout_ = setup.DSL().AddSet(DataSet::XYMESH, setname, "TI");
  if (dAout_ == 0) return Analysis::ERR;
  if (outfile != 0) outfile->AddDataSet( dAout_ );
  MetaData md(dAout_->Meta().Name(), "TIcurve");
  for (Iarray::const_iterator it = nskip_.begin(); it != nskip_.end(); ++it) {
    md.SetIdx( *it );
    DataSet* ds = setup.DSL().AddSet(DataSet::XYMESH, md);
    if (ds == 0) return Analysis::ERR;
    ds->SetLegend( md.Name() + "_Skip" + integerToString(*it) );
    if (curveout != 0) curveout->AddDataSet( ds );
    curve_.push_back( ds );
  }

  mprintf("    TI: Calculating TI using Gaussian quadrature with %zu points.\n",
          quad_.size());
  mprintf("\t%6s %8s %8s %s\n", "Point", "Abscissa", "Weight", "SetName");
  for (unsigned int i = 0; i != quad_.size(); i++)
    mprintf("\t%6i %8.5f %8.5f %s\n", i, quad_[i], wgt_[i], input_dsets_[i]->legend());
  if (nskip_.front() > 0) {
    mprintf("\tSkipping first");
    for (Iarray::const_iterator it = nskip_.begin(); it != nskip_.end(); ++it)
      mprintf(" %i", *it);
    mprintf(" data points for <DV/DL> calc.\n");
  }
  mprintf("\tResults saved in set '%s'\n", dAout_->legend());
  mprintf("\tTI curve(s) saved in set(s)");
  for (DSarray::const_iterator ds = curve_.begin(); ds != curve_.end(); ++ds)
    mprintf(" '%s'", (*ds)->legend());
  mprintf("\n");
  if (outfile != 0) mprintf("\tResults written to '%s'\n", outfile->DataFilename().full());
  if (curveout!= 0) mprintf("\tTI curve written to '%s'\n", curveout->DataFilename().full());
  return Analysis::OK;
}

Analysis::RetType Analysis_TI::Analyze() {
  Darray sum(nskip_.size(), 0.0);
  DataSet_Mesh& DA = static_cast<DataSet_Mesh&>( *dAout_ );
  Iarray lastSkipPoint; // Points after which averages can be recorded
  for (Iarray::const_iterator it = nskip_.begin(); it != nskip_.end(); ++it)
    lastSkipPoint.push_back( *it - 1 );
  // Run for multiple skip values, helps test convergences.
  for (unsigned int idx = 0; idx != input_dsets_.size(); idx++) {
    DataSet_1D const& ds = static_cast<DataSet_1D const&>( *(input_dsets_[idx]) );
    if (ds.Size() < 1) {
      mprinterr("Error: Set '%s' is empty.\n", ds.legend());
      return Analysis::ERR;
    }
    // Determine if skip values are valid for this set.
    Darray Npoints; // Number of points after skipping
    for (Iarray::const_iterator it = nskip_.begin(); it != nskip_.end(); ++it) {
      int np = (int)ds.Size() - *it;
      if (np < 1) {
        mprinterr("Error: Skipped too many points (set '%s' size is %zu)\n",ds.legend(),ds.Size());
        return Analysis::ERR;
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
      //mprintf("\t<DV/DL>=%g\n", avg);
      DataSet_Mesh& CR = static_cast<DataSet_Mesh&>( *(curve_[j]) );
      CR.AddXY(quad_[idx], avg[j]);
      sum[j] += (wgt_[idx] * avg[j]);
    }
  }
  for (unsigned int j = 0; j != nskip_.size(); j++)
    DA.AddXY(nskip_[j], sum[j]);

  return Analysis::OK;
}
