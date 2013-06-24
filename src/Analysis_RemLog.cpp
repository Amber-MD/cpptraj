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
  mprintf("   REMLOG: %s\n", remlog_->Legend().c_str());

  return Analysis::OK;
}

Analysis::RetType Analysis_RemLog::Analyze() {
  return Analysis::OK;
}

