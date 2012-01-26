#ifndef INC_ANALYSISLIST_H
#define INC_ANALYSISLIST_H
#include "Analysis.h"
/// Hold all analyses to be performed.
class AnalysisList {
    std::vector<Analysis*> analysisList;
    int Nanalysis;
    int debug;
  public: 
    AnalysisList();
    ~AnalysisList();
    
    void SetDebug(int);
    int AddAnalysis(ArgList &);
    int Setup(DataSetList*,ParmFileList*);
    void Analyze(DataFileList*); 
};
#endif
