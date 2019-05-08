#ifndef INC_ANALYSIS_CRDFLUCT_H
#define INC_ANALYSIS_CRDFLUCT_H
#include "Analysis.h"
#include "DataSet_Coords.h"
class Analysis_CrdFluct : public Analysis {
  public:
    Analysis_CrdFluct();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_CrdFluct(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    void CalcBfactors( Frame const&, Frame const&, double, DataSet& );

    DataSet_Coords* coords_;
    AtomMask mask_;
    Frame sum_;  ///< Temp. space for calculating average
    Frame sum2_; ///< Temp. space for calculating average of squares
    typedef std::vector<DataSet*> SetList;
    SetList outSets_;
    bool bfactor_;
    int windowSize_;
};
#endif
