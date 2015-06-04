#ifndef INC_ANALYSIS_MULTICURVE_H
#define INC_ANALYSIS_MULTICURVE_H
#include "Analysis.h"
#include "Array1D.h"
class Analysis_Multicurve : public Analysis {
  public:
    Analysis_Multicurve() : masterDSL_(0), masterDFL_(0), debug_(0) {}
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Multicurve(); }
    static void Help();

    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    Array1D inputDsets_;
    ArgList args_;
    DataSetList* masterDSL_;
    DataFileList* masterDFL_;
    int debug_;
};
#endif
