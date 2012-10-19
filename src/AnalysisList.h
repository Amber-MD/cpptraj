#ifndef INC_ANALYSISLIST_H
#define INC_ANALYSISLIST_H
#include "Analysis.h"
/// Hold all analyses to be performed.
class AnalysisList {
  public: 
    AnalysisList();
    ~AnalysisList();
    
    void SetDebug(int);
    int AddAnalysis(ArgList &);
    int Setup(DataSetList*,TopologyList*);
    void Analyze(DataFileList*);
    void List();
  private:
    typedef std::vector<Analysis*> aListType;
    typedef std::vector<Analysis*>::iterator aListIt;
    std::vector<Analysis*> analysisList_;
    int debug_;
};
#endif
