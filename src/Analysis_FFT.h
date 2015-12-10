#ifndef INC_ANALYSIS_FFT_H
#define INC_ANALYSIS_FFT_H
#include "Analysis.h"
#include "Array1D.h"
/// Calculate FFT of dataset(s)
class Analysis_FFT : public Analysis {
  public:
    Analysis_FFT();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_FFT(); }
    void Help() const;
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    Array1D input_dsets_;
    Array1D output_dsets_;
    double dt_; ///< Sampling interval (timestep)
};
#endif
