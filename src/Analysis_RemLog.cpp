#include "Analysis_RemLog.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"
#include "Analysis_Lifetime.h"
#include "StringRoutines.h" // integerToString
#include "DataSet_Mesh.h" // Regression

Analysis_RemLog::Analysis_RemLog() :
  debug_(0),
  calculateStats_(false),
  calculateLifetimes_(false),
  printIndividualTrips_(false), 
  remlog_(0),
  mode_(NONE),
  lifetimes_(0),
  statsout_(0),
  reptime_(0),
  calcRepFracSlope_(0),
  repFracSlope_(0)
{}

void Analysis_RemLog::Help() {
  mprintf("\t{<remlog dataset> | <remlog filename>} [out <filename>] [crdidx | repidx]\n"
          "\t[stats [statsout <file>] [printtrips] [reptime <file>]] [lifetime <file>]\n"
          "\t[reptimeslope <n> reptimeslopeout <file>] [acceptout <file>] [name <setname>]\n"
          "    crdidx: Print coordinate index vs exchange; output sets contain replica indices.\n"
          "    repidx: Print replica index vs exchange; output sets contain coordinate indices.\n"
          "  Analyze previously read in replica log data.\n");
}

// Analysis_RemLog::Setup()
Analysis::RetType Analysis_RemLog::Setup(ArgList& analyzeArgs, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  // Get remlog dataset
  std::string remlogName = analyzeArgs.GetStringNext();
  if (remlogName.empty()) {
    mprinterr("Error: no remlog data set or file name specified.\n");
    return Analysis::ERR;
  }
  // Check if data set exists
  remlog_ = (DataSet_RemLog*)DSL->FindSetOfType( remlogName, DataSet::REMLOG );
  if (remlog_ == 0) {
    mprinterr("Error: remlog data with name %s not found.\n", remlogName.c_str());
    return Analysis::ERR;
  }
  if (remlog_->Size() < 1 || remlog_->NumExchange() < 1) {
    mprinterr("Error: remlog data set appears to be empty.\n");
    return Analysis::ERR;
  }
  acceptout_ = DFL->AddCpptrajFile( analyzeArgs.GetStringKey("acceptout"), "replica acceptance",
                                      DataFileList::TEXT, true );
  if (acceptout_ == 0) return Analysis::ERR;
  lifetimes_ = DFL->AddCpptrajFile( analyzeArgs.GetStringKey("lifetime"), "remlog lifetimes" );
  calculateLifetimes_ = (lifetimes_ != 0);
  calculateStats_ = analyzeArgs.hasKey("stats");
  if (calculateStats_) {
    statsout_ = DFL->AddCpptrajFile( analyzeArgs.GetStringKey("statsout"), "remlog stats",
                                       DataFileList::TEXT, true );
    reptime_ = DFL->AddCpptrajFile( analyzeArgs.GetStringKey("reptime"), "replica times",
                                      DataFileList::TEXT, true );
    if (statsout_ == 0 || reptime_ == 0) return Analysis::ERR;
  }
  calcRepFracSlope_ = analyzeArgs.getKeyInt("reptimeslope", 0);
  std::string rfs_name = analyzeArgs.GetStringKey("reptimeslopeout");
  if (!calculateStats_) {
    calcRepFracSlope_ = 0;
    rfs_name.clear();
  }
  if ( (calcRepFracSlope_ > 0) != (!rfs_name.empty()) ) {
    mprinterr("Error: Both reptimeslope and reptimeslopeout must be specified.\n");
    return Analysis::ERR;
  }
  repFracSlope_ = DFL->AddCpptrajFile( rfs_name, "replica fraction slope" );
  printIndividualTrips_ = analyzeArgs.hasKey("printtrips");
  // Get mode
  if (analyzeArgs.hasKey("crdidx"))
    mode_ = CRDIDX;
  else if (analyzeArgs.hasKey("repidx"))
    mode_ = REPIDX;
  else
    mode_ = NONE;
  const char* def_name = 0;
  const char* yaxis = 0;
  if (mode_ == CRDIDX) {
    def_name = "repidx";
    yaxis = "ylabel CrdIdx";
  } else if (mode_ == REPIDX) {
    def_name = "crdidx";
    yaxis = "ylabel RepIdx";
  }
  // Set up an output set for each replica
  DataFile* dfout = 0;
  if (mode_ != NONE) {
    // Get output filename
    std::string outname = analyzeArgs.GetStringKey("out");
    if (!outname.empty()) {
      dfout = DFL->AddDataFile( outname, analyzeArgs );
      if (dfout == 0 ) return Analysis::ERR;
      if (yaxis != 0 ) dfout->ProcessArgs(yaxis);
    }
    std::string dsname = analyzeArgs.GetStringKey("name");
    if (dsname.empty())
      dsname = DSL->GenerateDefaultName(def_name);
    MetaData md(dsname);
    for (int i = 0; i < (int)remlog_->Size(); i++) {
      md.SetIdx(i+1);
      DataSet_integer* ds = (DataSet_integer*)DSL->AddSet(DataSet::INTEGER, md);
      if (ds == 0) return Analysis::ERR;
      outputDsets_.push_back( (DataSet*)ds );
      if (dfout != 0) dfout->AddDataSet( (DataSet*)ds );
      ds->Resize( remlog_->NumExchange() ); 
    }
  }
  mprintf("   REMLOG: %s, %i replicas, %i exchanges\n", remlog_->legend(),
          remlog_->Size(), remlog_->NumExchange());
  if (mode_ == CRDIDX)
    mprintf("\tGetting coordinate index vs exchange.\n");
  else if (mode_ == REPIDX)
    mprintf("\tGetting replica index vs exchange.\n");
  if (mode_ != NONE && dfout != 0)
    mprintf("\tOutput is to %s\n", dfout->DataFilename().base());
  if (calculateStats_) {
    mprintf("\tGetting replica exchange stats, output to %s\n", statsout_->Filename().full());
    if (printIndividualTrips_)
      mprintf("\tIndividual round trips will be printed.\n");
    mprintf("\tWriting time spent at each replica to %s\n", reptime_->Filename().full());
  }
  if (calculateLifetimes_)
    mprintf("\tThe lifetime of each crd at each replica will be calculated.\n");
  if (acceptout_ != 0)
    mprintf("\tOverall exchange acceptance % will be written to %s\n",
            acceptout_->Filename().full());

  return Analysis::OK;
}

// Analysis_RemLog::Analyze()
Analysis::RetType Analysis_RemLog::Analyze() {
  if (remlog_->Size() < 1) {
    mprinterr("Error: No replicas in remlog data '%s'\n", remlog_->legend());
    return Analysis::ERR;
  }
  int Ndims = remlog_->DimTypes().Ndims();
  mprintf("\t'%s' %i replicas, %i exchanges, %i dims.\n", remlog_->legend(),
         remlog_->Size(), remlog_->NumExchange(), Ndims);
  // Set up arrays for tracking replica stats.
  std::vector<RepStats> DimStats;
  std::vector<TripStats> DimTrips;
  for (int i = 0; i != Ndims; i++) {
    DimStats.push_back( RepStats(remlog_->Size()) );
    if (calculateStats_)
      DimTrips.push_back( TripStats(remlog_->Size()) );
  }
  std::vector< std::vector<int> > replicaFrac;
  if (calculateStats_) {
    replicaFrac.resize( remlog_->Size() ); // [replica][crdidx]
    for (std::vector< std::vector<int> >::iterator it = replicaFrac.begin();
                                                   it != replicaFrac.end(); ++it)
      it->resize( remlog_->Size(), 0 );
  }
  // Variables for calculating replica lifetimes
  Analysis_Lifetime Lifetime;
  Array1D dsLifetime;
  std::vector< std::vector<DataSet_integer> > series; // 2D - repidx, crdidx
  if (calculateLifetimes_) {
    mprintf("\tData size used for lifetime analysis= %zu bytes.\n",
            remlog_->Size() * remlog_->Size() * remlog_->NumExchange() * sizeof(int));
    series.resize( remlog_->Size() );
    for (unsigned int i = 0; i < remlog_->Size(); i++) {
      series[i].resize( remlog_->Size() );
      for (unsigned int j = 0; j < remlog_->Size(); j++) {
        series[i][j].Resize( remlog_->NumExchange() );
        series[i][j].SetLegend("Rep"+integerToString(i+1)+",Crd"+integerToString(j+1));
        dsLifetime.push_back( (DataSet_1D*)&(series[i][j]) );
      }
    }
    if (Lifetime.ExternalSetup( dsLifetime, lifetimes_ ) == Analysis::ERR) {
      mprinterr("Error: Could not set up remlog lifetime analysis.\n");
      return Analysis::ERR;
    }
  }

  DataSet_Mesh mesh;
  if ( calcRepFracSlope_ > 0 ) {
    mesh.CalculateMeshX( remlog_->Size(), 1, remlog_->Size() );
    repFracSlope_->Printf("%-8s", "#Exchg");
    for (int crdidx = 0; crdidx < (int)remlog_->Size(); crdidx++)
      repFracSlope_->Printf("  C%07i_slope C%07i_corel", crdidx + 1, crdidx + 1);
    repFracSlope_->Printf("\n");
  }

  ProgressBar progress( remlog_->NumExchange() );
  for (int frame = 0; frame < remlog_->NumExchange(); frame++) {
    progress.Update( frame );
    for (int replica = 0; replica < (int)remlog_->Size(); replica++) {
      DataSet_RemLog::ReplicaFrame const& frm = remlog_->RepFrame( frame, replica );
      int crdidx = frm.CoordsIdx() - 1;
      int repidx = frm.ReplicaIdx() - 1;
      int dim = frm.Dim();
      // Exchange acceptance.
      // NOTE: Because currently the direction of the attempt is not always
      //       known unless the attempt succeeds for certain remlog types,
      //       the results will be skewed if dimension size is 2 since in that
      //       case the left partner is the right partner.
      if (replica == 0) DimStats[dim].attempts_++; // Assume same # attempts for every rep in dim
      if (frm.Success()) {
        if (frm.PartnerIdx() - 1 == remlog_->ReplicaInfo()[replica][dim].RightID())
          DimStats[dim].acceptUp_[replica]++;
        else // Assume down
          DimStats[dim].acceptDown_[replica]++;
      }
      if (mode_ == CRDIDX) {
        DataSet_integer& ds = static_cast<DataSet_integer&>( *(outputDsets_[repidx]) );
        ds[frame] = frm.CoordsIdx();
      } else if (mode_ == REPIDX) {
        DataSet_integer& ds = static_cast<DataSet_integer&>( *(outputDsets_[crdidx]) );
        ds[frame] = frm.ReplicaIdx();
      }
      if (calculateLifetimes_)
        series[repidx][crdidx][frame] = 1;
      if (calculateStats_) {
        TripStats& trip = static_cast<TripStats&>( DimTrips[dim] );
        // Fraction spent at each replica
        replicaFrac[repidx][crdidx]++;
        // Replica round-trip calculation
        if (trip.status_[crdidx] == UNKNOWN) {
          if (remlog_->ReplicaInfo()[repidx][dim].Location() == DataSet_RemLog::BOTTOM) {
            trip.status_[crdidx] = HIT_BOTTOM;
            trip.bottom_[crdidx] = frame;
          }
        } else if (trip.status_[crdidx] == HIT_BOTTOM) {
          if (remlog_->ReplicaInfo()[repidx][dim].Location() == DataSet_RemLog::TOP)
            trip.status_[crdidx] = HIT_TOP;
        } else if (trip.status_[crdidx] == HIT_TOP) {
          if (remlog_->ReplicaInfo()[repidx][dim].Location() == DataSet_RemLog::BOTTOM) {
            int rtrip = frame - trip.bottom_[crdidx];
            if (printIndividualTrips_)
              statsout_->Printf("[%i] CRDIDX %i took %i exchanges to travel"
                               " up and down (exch %i to %i)\n",
                               replica, crdidx+1, rtrip, trip.bottom_[crdidx]+1, frame+1);
            trip.roundTrip_[crdidx].AddElement( rtrip );
            trip.status_[crdidx] = HIT_BOTTOM;
            trip.bottom_[crdidx] = frame;
          }
        }
      }
    } // END loop over replicas
    if (calcRepFracSlope_ > 0 && frame > 0 && (frame % calcRepFracSlope_) == 0) {
      repFracSlope_->Printf("%8i", frame+1);
      for (int crdidx = 0; crdidx < (int)remlog_->Size(); crdidx++) {
        for (int replica = 0; replica < (int)remlog_->Size(); replica++)
          mesh.SetY(replica, (double)replicaFrac[replica][crdidx] / (double)frame);
        double slope, intercept, correl;
        mesh.LinearRegression(slope, intercept, correl, true);
        repFracSlope_->Printf("  %14.7g %14.7g", slope * 100.0, correl);
                //frame+1, crdidx, slope * 100.0, intercept * 100.0, correl
      }
      repFracSlope_->Printf("\n");
    }
  } // END loop over exchanges
  // Exchange acceptance calc.
  for (int dim = 0; dim != Ndims; dim++) {
    // Assume number of exchange attempts is actually /2 since in Amber
    // attempts alternate up/down.
    acceptout_->Printf("DIMENSION %i\n", dim+1);
    if (debug_ > 0) {
    for (int replica = 0; replica != (int)remlog_->Size(); replica++)
      mprintf("Rep %i attempts %i up %i down %i\n", replica, DimStats[dim].attempts_, DimStats[dim].acceptUp_[replica], DimStats[dim].acceptDown_[replica]);
    }
    acceptout_->Printf("%-8s %8s %8s\n", "#Replica", "%UP", "%DOWN");
    double exchangeAttempts = (double)DimStats[dim].attempts_ / 2.0;
    for (int replica = 0; replica != (int)remlog_->Size(); replica++)
      acceptout_->Printf("%8i %8.3f %8.3f\n", replica+1,
            ((double)DimStats[dim].acceptUp_[replica] / exchangeAttempts) * 100.0,
            ((double)DimStats[dim].acceptDown_[replica] / exchangeAttempts) * 100.0);
  }
  if (calculateStats_) {
    statsout_->Printf("# %i replicas, %i exchanges.\n", remlog_->Size(), remlog_->NumExchange());
    for (int dim = 0; dim != Ndims; dim++) {
      if (Ndims > 1)
        statsout_->Printf("#Dim%i Round-trip stats:\n", dim+1);
      else
        statsout_->Printf("#Round-trip stats:\n");
      statsout_->Printf("#%-8s %8s %12s %12s %12s %12s\n", "CRDIDX", "RndTrips", 
                       "AvgExch.", "SD_Exch.", "Min", "Max");
      unsigned int idx = 1;
      for (DSI_array::const_iterator rt = DimTrips[dim].roundTrip_.begin();
                                     rt != DimTrips[dim].roundTrip_.end(); ++rt)
      {
        double stdev = 0.0;
        double avg = rt->Avg( stdev );
        statsout_->Printf("%-8u %8i %12.4f %12.4f %12.0f %12.0f\n", 
                          idx++, rt->Size(), avg, stdev, rt->Min(), rt->Max());
      }
    }
    reptime_->Printf("#Percent time spent at each replica:\n%-8s", "#Replica");
    for (int crd = 0; crd < (int)remlog_->Size(); crd++)
      reptime_->Printf(" CRD_%04i", crd + 1);
    reptime_->Printf("\n");
    double dframes = (double)remlog_->NumExchange();
    for (int replica = 0; replica < (int)remlog_->Size(); replica++) {
      reptime_->Printf("%8i", replica+1);
      for (int crd = 0; crd < (int)remlog_->Size(); crd++)
        reptime_->Printf(" %8.3f", ((double)replicaFrac[replica][crd] / dframes) * 100.0);
      reptime_->Printf("\n");
    }
  }
  if (calculateLifetimes_) {
    mprintf("\tCalculating remlog lifetimes:\n");
    Lifetime.Analyze();
  }
  return Analysis::OK;
}
