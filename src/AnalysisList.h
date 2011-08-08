#ifndef INC_ANALYSISLIST_H
#define INC_ANALYSISLIST_H
#include "Analysis.h"
class AnalysisList {
    Analysis** analysisList;
    int Nanalysis;
    int debug;
  public: 
    AnalysisList();
    ~AnalysisList();
    
    void SetDebug(int);
    int Add(ArgList *A);
    int Setup(DataSetList*); 
};
#endif
