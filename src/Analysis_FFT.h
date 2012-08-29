#ifndef INC_ANALYSIS_FFT_H
#define INC_ANALYSIS_FFT_H
#include "Analysis.h"
/// Calculate FFT of dataset(s)
class Analysis_FFT : public Analysis {
  public:
    Analysis_FFT();

    int Setup(DataSetList*);
    int Analyze();
    void Print(DataFileList*);
  private:
    DataSetList input_dsets_;
    std::string outfilename_;
    std::string setname_;
    std::vector<DataSet*> output_dsets_;
    int maxsize_;
    double dt_;
    double f0_; 
};
#endif
