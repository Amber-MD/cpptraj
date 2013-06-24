#include "Analysis_RemLog.h"
#include "CpptrajStdio.h"

Analysis_RemLog::Analysis_RemLog() : remlog_(0) {}

void Analysis_RemLog::Help() {
  mprintf("\t{<remlog dataset> | <remlog filename>}\n");
}

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
  // Get output filename
  DataFile* dfout = 0;
  std::string outname = analyzeArgs.GetStringKey("out");
  if (!outname.empty()) {
    dfout = DFLin->AddDataFile( outname, analyzeArgs );
    if (dfout == 0 ) return Analysis::ERR;
  }
  // Set up an output set for each replica
  std::string dsname = analyzeArgs.GetStringNext();
  if (dsname.empty())
    dsname = datasetlist->GenerateDefaultName("crdidx");
  for (int i = 0; i < remlog_->Size(); i++) {
    outputDsets_.push_back( datasetlist->AddSetIdx(DataSet::INT, dsname, i+1) );
    if (outputDsets_.back() == 0) return Analysis::ERR;
    if (dfout != 0) dfout->AddSet( outputDsets_.back() );
  }
  mprintf("   REMLOG: %s, %i replicas, %i exchanges\n", remlog_->Legend().c_str(),
          remlog_->Size(), remlog_->NumExchange());

  return Analysis::OK;
}

Analysis::RetType Analysis_RemLog::Analyze() {
  // For now, just write Coordinate index vs time for each replica
  int replica = 0;
  for (DataSet_RemLog::ensemble_it rep = remlog_->begin();
                                   rep != remlog_->end(); ++rep)
  {
    int frame = 0;
    for (DataSet_RemLog::replica_it frm = (*rep).begin();
                                    frm != (*rep).end(); ++frm)
    {
      int crdidx = (*frm).CoordsIdx();
      outputDsets_[replica]->Add( frame++, &crdidx );
      //int repidx = (*frm).ReplicaIdx();
      //#outputDsets_[(*frm).CoordsIdx()-1]->Add( frame++, &repidx );
    }
    replica++;
  }
  
  return Analysis::OK;
}

