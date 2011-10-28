#ifndef INC_ANALYSIS_H
#define INC_ANALYSIS_H
/// Class: Analysis
/// Base class that all analysis routines will inherit.
/// Analysis occurs after trajectories are read and data sets populated.
/// Analysis operates on those data sets.
#include "ArgList.h"
#include "DataSetList.h"
#include "DataFileList.h"
class Analysis {
  protected:
    int debug;
    ArgList analyzeArgs;
  public:
    Analysis();
    virtual ~Analysis();

    bool noSetup;

    void SetArg(const ArgList &);
    void SetDebug(int);
    const char *AnalysisCommand();              // Print the command that calls the analysis
    const char *CmdLine();           // Print the entire argument line 
  
    virtual int Setup(DataSetList*)   {return 1;}
    virtual int Analyze()             {return 1;}
    virtual void Print(DataFileList*) {return;  }
};
#endif
