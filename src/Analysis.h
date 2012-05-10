#ifndef INC_ANALYSIS_H
#define INC_ANALYSIS_H
#include "ArgList.h"
#include "DataSetList.h"
#include "DataFileList.h"
#include "TopologyList.h"
// Class: Analysis
/// Base class that all analysis routines will inherit.
/** Analysis occurs after trajectories are read and data sets populated.
  * Analysis operates on those data sets.
  */
class Analysis {
  public:
    Analysis();
    virtual ~Analysis();

    void SetArg(const ArgList &);
    void SetDebug(int);
    void SetParm(TopologyList*);
    const char *AnalysisCommand();   ///< Print the command that calls the analysis
    const char *CmdLine();           ///< Print the entire argument line
    void SetSetup(bool sIn) { isSetup_ = sIn; }
    bool IsSetup() { return isSetup_; }
  
    virtual int Setup(DataSetList*)   {return 1;}
    virtual int Analyze()             {return 1;}
    virtual void Print(DataFileList*) {return;  }
  protected:
    int debug_;
    ArgList analyzeArgs_;
    // NOTE: Only used for ptraj analysis
    Topology *analyzeParm_;
  private:
    bool isSetup_; ///< True if analysis could be setup
};
#endif
