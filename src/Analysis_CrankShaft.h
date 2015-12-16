#ifndef INC_ANALYSIS_CRANKSHAFT_H
#define INC_ANALYSIS_CRANKSHAFT_H
#include "Analysis.h"
#include "DataSet_1D.h"
class Analysis_CrankShaft : public Analysis {
  public: 
    Analysis_CrankShaft();

    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_CrankShaft(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    enum CStype { ANGLE=0, DISTANCE };
    static const char* CSstring[];
    static const char* distance_ss_2D[][6];
    static const char* torsion_ss_2D[][6];

    CpptrajFile* frame_vs_bin_; ///< File for frame vs bin number
    CpptrajFile* results_;      ///< File for crankshaft results
    int debug_;
    int start_;
    int stop_;
    int offset_;
    CStype type_;
    DataSet_1D* scalar1_;
    DataSet_1D* scalar2_;
    std::string info_;
};
#endif
