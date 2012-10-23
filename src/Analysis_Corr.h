#ifndef INC_ANALYSIS_CORR_H
#define INC_ANALYSIS_CORR_H
#include "Analysis.h"
// Class: Analysis_Corr
/// Calculate autocorrelation or correlation of dataset(s).
class Analysis_Corr : public Analysis {
  public:
    Analysis_Corr();

    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Corr(); }
    static void Help();

    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,int);
    Analysis::RetType Analyze();
    void Print(DataFileList*);
  private:
    DataSet *D1_;
    DataSet *D2_;
    int lagmax_;
    DataSet* Ct_;
    std::string outfilename_;
    bool usefft_;
    bool calc_covar_;
};
#endif
