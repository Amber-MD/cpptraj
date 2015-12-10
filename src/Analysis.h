#ifndef INC_ANALYSIS_H
#define INC_ANALYSIS_H
#include "DispatchObject.h"
#include "ArgList.h"
#include "AnalysisState.h"
/// The abstract base class that all other analyses inherit.
/** Analysis occurs after trajectories are read and data sets populated.
  * Analysis operates on those data sets.
  */
class Analysis : public DispatchObject {
  public:
    /// Constructor
    Analysis() : DispatchObject(ANALYSIS) {}
    /// Constructor - override ANALYSIS (e.g. HIDDEN)
    Analysis(DispatchObject::Otype o) : DispatchObject(o) {}
    /// Enumerate potential return stats from Setup and Analyze.
    enum RetType { OK = 0, ERR };
    /// Destructor - virtual since this class is inherited
    virtual ~Analysis() {}
    /// Set up Analysis 
    virtual RetType Setup(ArgList&, AnalysisSetup&, int) = 0;
    /// Execute Analysis
    virtual RetType Analyze() = 0;
};
#endif
