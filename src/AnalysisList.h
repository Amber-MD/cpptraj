#ifndef INC_ANALYSISLIST_H
#define INC_ANALYSISLIST_H
#include "Analysis.h"
/// Hold all Analyses to be performed.
class AnalysisList {
  public:
    /// Analysis setup status
    enum AnalysisStatusType { NO_SETUP = 0, SETUP, INACTIVE };

    AnalysisList();
    ~AnalysisList();
    void Clear(); 
    void SetDebug(int d) { debug_ = d; }
    int Debug() const { return debug_; }
    int AddAnalysis(Analysis*, ArgList&, AnalysisSetup&);
    int DoAnalyses();
    void List() const;
    bool Empty() const { return analysisList_.empty(); }
    unsigned int size() const { return analysisList_.size(); }
#   ifdef MPI
    Analysis& Ana(int i)                   { return *(analysisList_[i].ptr_); }
    AnalysisStatusType Status(int i) const { return analysisList_[i].status_; }
    ArgList const& Args(int i)       const { return analysisList_[i].args_;   }
#   endif
  private:
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
