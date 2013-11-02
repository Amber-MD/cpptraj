#ifndef INC_ANALYSIS_FFT_H
#define INC_ANALYSIS_FFT_H
#include "Analysis.h"
#include "DataSet_1D.h"
/// Calculate FFT of dataset(s)
class Analysis_FFT : public Analysis {
  public:
    Analysis_FFT();

    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_FFT(); }
    static void Help();

    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    DataFile* outfile_;
    DataSetList input_dsets_;
    std::vector<DataSet_1D*> output_dsets_;
    size_t maxsize_;
    double dt_;
};
#endif
