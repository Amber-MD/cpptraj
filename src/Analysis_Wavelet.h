#ifndef INC_ANALYSIS_WAVELET_H
#define INC_ANALYSIS_WAVELET_H
#include "Analysis.h"
#include "ComplexArray.h"
/// Perform wavelet analysis
/** \author Original code: Zahra Heidari
  * \author Implemented in CPPTRAJ by Dan Roe
  */
class Analysis_Wavelet : public Analysis {
  public:
    Analysis_Wavelet();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Wavelet(); }
    void Help() const;
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
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
