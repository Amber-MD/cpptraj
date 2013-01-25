#ifndef INC_ANALYSIS_MATRIX_H
#define INC_ANALYSIS_MATRIX_H
#include "Analysis.h"
#include "DataSet_Matrix.h"
#include "DataSet_Modes.h"
class Analysis_Matrix : public Analysis {
  public:
    Analysis_Matrix();

    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Matrix(); }
    static void Help();


    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    DataSet_Matrix* matrix_;
    DataSet_Modes* modes_;
    std::string outfilename_;
    std::string outthermo_;
    int nevec_;
    bool thermopt_;
    bool reduce_;
    bool eigenvaluesOnly_;
};
#endif
