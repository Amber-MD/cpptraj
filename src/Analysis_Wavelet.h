#ifndef INC_ANALYSIS_WAVELET_H
#define INC_ANALYSIS_WAVELET_H
#include "Analysis.h"
#include "ComplexArray.h"
/// Perform wavelet analysis
class Analysis_Wavelet : public Analysis {
  public:
    Analysis_Wavelet();
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Wavelet(); }
    static void Help();
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    enum WaveletType { W_MORLET = 0, W_PAUL, W_NONE };
    // Wavelet functions
    ComplexArray F_Morlet(std::vector<int> const&, double) const;
    ComplexArray F_Paul(std::vector<int> const&, double) const;
    /** Function prototype for wavelet: Fxn( complexOutput, x ) */
    //typedef void (*WaveletFxnType)(ComplexArray&, std::vector<int> const&, double);

    struct WaveletToken { const char* key_; const char* description_; };
    static const WaveletToken Tokens_[];

    AtomMask mask_;
    DataSet_Coords* coords_;
    DataSet* output_;
    double S0_;
    double ds_;
    double correction_;
    double chival_;
    WaveletType wavelet_type_;
    int nb_;
};
#endif
