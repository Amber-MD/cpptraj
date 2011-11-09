#ifndef INC_ANALYSIS_PTRAJANALYSIS
#define INC_ANALYSIS_PTRAJANALYSIS
#include "Analysis.h"
#include "ptraj_analyze.h"
#include "ptraj_arg.h"
// Class: PtrajAnalysis
/// Wrapper for ptraj analysis functions in ptraj_analyze.c
class PtrajAnalysis: public Analysis {
    analyzeInformation *analyzeinfo;
    argStackType *argumentStack;
  public: 
    PtrajAnalysis();
    ~PtrajAnalysis();

    int Setup(DataSetList*);
    int Analyze();
    void Print(DataFileList*);
};
#endif
