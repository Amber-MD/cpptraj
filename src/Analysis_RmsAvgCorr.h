#ifndef INC_ANALYSIS_RMSAVGCORR_H
#define INC_ANALYSIS_RMSAVGCORR_H
#include "Analysis.h"
#include "DataSet_Coords.h"
// Class: Analysis_RmsAvgCorr
/// Calculate rmsd using running avg structures
class Analysis_RmsAvgCorr: public Analysis {
  public:
    Analysis_RmsAvgCorr();

    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_RmsAvgCorr(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    AtomMask tgtMask_;
    CpptrajFile* separate_;
    DataSet_Coords* coords_;
    DataSet* Ct_;
    DataSet* Csd_;
    Frame refFrame_;
    int maxwindow_;
    int lagOffset_;
    bool useMass_;
    bool useFirst_; ///< If true, use first running-avgd frame as reference.
};
#endif  
