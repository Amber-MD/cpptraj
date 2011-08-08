#ifndef INC_ANALYSIS_H
#define INC_ANALYSIS_H
/// Class: Analysis
/// Base class that all analysis routines will inherit.
/// Analysis occurs after trajectories are read and data sets pupulated.
/// Analysis operates on those data sets.
#include "ArgList.h"
#include "DataSetList.h"
class Analysis {
  protected:
    int debug;
    ArgList *analyzeArg;
  public:
    Analysis();
    virtual ~Analysis();

    bool noSetup;

    void SetArg(ArgList *);
    void SetDebug(int);
    char *Name(); 
  
    virtual int Setup(DataSetList *) {return 1;}
    virtual int Analyze()            {return 1;}
    virtual void Print()             {return;  }
};
#endif
