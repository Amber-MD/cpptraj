#ifndef INC_ANALYSISLIST_H
#define INC_ANALYSISLIST_H
#include "Analysis.h"
/// Hold all Analyses to be performed.
class AnalysisList {
  public:
    AnalysisList();
    ~AnalysisList();
    void Clear(); 
    void SetDebug(int d) { debug_ = d; }
    int Debug() const { return debug_; }
    int AddAnalysis(Analysis*, ArgList&, AnalysisSetup&);
    int DoAnalyses();
    void List() const;
    bool Empty() const { return analysisList_.empty(); }
  private:
    /// Analysis setup status
    enum AnalysisStatusType { NO_SETUP = 0, SETUP, INACTIVE };
    struct AnaHolder {
      Analysis* ptr_;             ///< Pointer to Analysis
      ArgList args_;              ///< Arguments associated with Analysis
      AnalysisStatusType status_; ///< Current Analysis status.
    };
    typedef std::vector<AnaHolder> Aarray;
    Aarray analysisList_; ///< List of Analyses
    int debug_;           ///< Default debug level for new Analyses
};
#endif
