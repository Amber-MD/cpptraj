#ifndef INC_ANALYSIS_RUNNNINGAVG_H
#define INC_ANALYSIS_RUNNNINGAVG_H
#include "Analysis.h"
class Analysis_RunningAvg : public Analysis {
  public:
    Analysis_RunningAvg();

    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_RunningAvg(); }
    static void Help();

    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    DataSetList dsets_;
    bool cumulative_;
    int window_;
    std::vector<DataSet*> outputData_;
};
#endif
