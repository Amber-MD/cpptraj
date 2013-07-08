#include "Analysis_RemLog.h"
#include "CpptrajStdio.h"
#include "DataSet_integer.h"
#include "DS_Math.h"

Analysis_RemLog::Analysis_RemLog() :
  calculateStats_(false), 
  remlog_(0),
  mode_(NONE)
{}

void Analysis_RemLog::Help() {
  mprintf("\t{<remlog dataset> | <remlog filename>} [out <filename>]"
          " [crdidx | repidx] [stats [statsout <file>]]\n");
}

// Analysis_RemLog::Setup()
Analysis::RetType Analysis_RemLog::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  // Get remlog dataset
  std::string remlogName = analyzeArgs.GetStringNext();
  if (remlogName.empty()) {
    mprinterr("Error: no remlog data set or file name specified.\n");
    return Analysis::ERR;
  }
  // Check if data set exists
  remlog_ = (DataSet_RemLog*)datasetlist->FindSetOfType( remlogName, DataSet::REMLOG );
  if (remlog_ == 0) {
    mprinterr("Error: remlog data with name %s not found.\n", remlogName.c_str());
    return Analysis::ERR;
  }
  if (remlog_->Size() < 1 || remlog_->NumExchange() < 1) {
    mprinterr("Error: remlog data set appears to be empty.\n");
    return Analysis::ERR;
  }
  calculateStats_ = analyzeArgs.hasKey("stats");
  if (calculateStats_) {
    if (statsout_.OpenWrite( analyzeArgs.GetStringKey("statsout") )) return Analysis::ERR;
  }
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
      dfout = DFLin->AddDataFile( outname, analyzeArgs );
      if (dfout == 0 ) return Analysis::ERR;
      if (yaxis != 0 ) dfout->ProcessArgs(yaxis);
    }
    std::string dsname = analyzeArgs.GetStringNext();
    if (dsname.empty())
      dsname = datasetlist->GenerateDefaultName(def_name);
    for (int i = 0; i < remlog_->Size(); i++) {
      DataSet_integer* ds = (DataSet_integer*)datasetlist->AddSetIdx(DataSet::INT, dsname, i+1);
      if (ds == 0) return Analysis::ERR;
      outputDsets_.push_back( (DataSet*)ds );
      if (dfout != 0) dfout->AddSet( (DataSet*)ds );
      ds->Resize( remlog_->NumExchange() ); 
    }
  }
  mprintf("   REMLOG: %s, %i replicas, %i exchanges\n", remlog_->Legend().c_str(),
          remlog_->Size(), remlog_->NumExchange());
  if (mode_ == CRDIDX)
    mprintf("\tGetting coordinate index vs exchange.\n");
  else if (mode_ == REPIDX)
    mprintf("\tGetting replica index vs exchange.\n");
  if (calculateStats_)
    mprintf("\tGetting replica exchange stats.\n");
  if (dfout != 0)
    mprintf("\tOutput is to %s\n", dfout->DataFilename().base());

  return Analysis::OK;
}

// Analysis_RemLog::Analyze()
Analysis::RetType Analysis_RemLog::Analyze() {
  enum RepStatusType { UNKNOWN = 0, HIT_BOTTOM, HIT_TOP };
  std::vector<int> replicaStatus;
  std::vector<int> replicaBottom;
  std::vector<DataSet_integer> roundTrip;
  std::vector< std::vector<int> > replicaFrac;
  if (calculateStats_) {
    replicaStatus.resize( remlog_->Size(), UNKNOWN );  
    replicaBottom.resize( remlog_->Size(), 0 );
    roundTrip.resize( remlog_->Size() );
    replicaFrac.resize( remlog_->Size() );
    for (std::vector< std::vector<int> >::iterator it = replicaFrac.begin();
                                                   it != replicaFrac.end(); ++it)
      (*it).resize( remlog_->Size(), 0 );
  }

  for (int frame = 0; frame < remlog_->NumExchange(); frame++) {
    for (int replica = 0; replica < remlog_->Size(); replica++) {
      DataSet_RemLog::ReplicaFrame const& frm = remlog_->RepFrame( frame, replica );
      int crdidx = frm.CoordsIdx();
      int repidx = frm.ReplicaIdx();
      if (mode_ == CRDIDX) {
        DataSet_integer& ds = static_cast<DataSet_integer&>( *(outputDsets_[repidx-1]) );
        ds[frame] = crdidx;
      } else if (mode_ == REPIDX) {
        DataSet_integer& ds = static_cast<DataSet_integer&>( *(outputDsets_[crdidx-1]) );
        ds[frame] = repidx;
      }
      if (calculateStats_) {
        // Fraction spent at each replica
        replicaFrac[repidx-1][crdidx-1]++;
        // Replica round-trip calculation
        if (replicaStatus[crdidx-1] == UNKNOWN) {
          if (repidx == 1) {
            replicaStatus[crdidx-1] = HIT_BOTTOM;
            replicaBottom[crdidx-1] = frame;
          }
        } else if (replicaStatus[crdidx-1] == HIT_BOTTOM) {
          if (repidx == remlog_->Size())
            replicaStatus[crdidx-1] = HIT_TOP;
        } else if (replicaStatus[crdidx-1] == HIT_TOP) {
          if (repidx == 1) {
            int rtrip = frame - replicaBottom[crdidx-1];
            mprintf("[%i] CRDIDX %i took %i exchanges to travel up and down (exch %i to %i)\n",
                    replica, crdidx, rtrip, replicaBottom[crdidx-1]+1, frame+1);
            roundTrip[crdidx-1].push_back( rtrip );
            replicaStatus[crdidx-1] = HIT_BOTTOM;
            replicaBottom[crdidx-1] = frame;
          }
        }
      }
    } // END loop over exchanges for replica
  } // END loop over replicas

  if (calculateStats_) {
    mprintf("#Round-trip stats:\n");
    for (std::vector<DataSet_integer>::iterator rt = roundTrip.begin();
                                                rt != roundTrip.end(); ++rt)
    {
      double stdev = 0.0;
      double avg = DS_Math::Avg( *rt, &stdev );
      mprintf("CRDIDX %u made %i round trips. %f +/- %f exchanges.\n", 
              rt - roundTrip.begin() + 1, (*rt).Size(), avg, stdev);
    }
   
    mprintf("#Percent time spent at each replica:\n%-8s", "#Replica");
    for (int crd = 0; crd < remlog_->Size(); crd++)
      mprintf(" CRD_%04i", crd + 1);
    mprintf("\n");
    double dframes = (double)remlog_->NumExchange();
    for (int replica = 0; replica < remlog_->Size(); replica++) {
      mprintf("%8i", replica+1);
      for (int crd = 0; crd < remlog_->Size(); crd++)
        mprintf(" %8.3f", ((double)replicaFrac[replica][crd] / dframes) * 100.0);
      mprintf("\n");
    }
  }
  return Analysis::OK;
}

