#ifndef INC_ANALYSIS_REMLOG_H
#define INC_ANALYSIS_REMLOG_H
#include "Analysis.h"
#include "DataSet_RemLog.h"
class Analysis_RemLog : public Analysis {
  public:
    Analysis_RemLog();
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_RemLog(); }
    static void Help();
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    DataSet_RemLog* remlog_;
};
#endif
