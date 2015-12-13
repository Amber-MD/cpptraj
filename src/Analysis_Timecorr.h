#ifndef INC_ANALYSIS_TIMECORR_H
#define INC_ANALYSIS_TIMECORR_H
#include "Analysis.h"
#include "DataSet_Vector.h"
#include "Corr.h"
/** \author Original Code by Alrun N. Koller & H. Gohlke
  * \author Adapted by DRR
  */
class Analysis_Timecorr : public Analysis {
  public:
    Analysis_Timecorr();

    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Timecorr(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    struct AvgResults {
      double avgr_;
      double rave_;
      double r3iave_;
      double r6iave_;
    };
    enum ModeType { AUTOCORR = 0, CROSSCORR };
    static const char* ModeString[];

    std::vector<double> CalculateAverages(DataSet_Vector const&, AvgResults&);
    void CalcCorr(int);
    void Normalize( DataSet*, int, double );

    double tstep_;
    double tcorr_;
    int order_;
    ModeType mode_;
    bool dplr_;
    bool norm_;
    bool drct_;
    bool ptrajformat_;
    ComplexArray data1_;
    ComplexArray data2_;
    DataSet_Vector* vinfo1_;
    DataSet_Vector* vinfo2_;
    DataSet* tc_c_;
    DataSet* tc_p_;
    DataSet* tc_r3r3_;
    CpptrajFile* outfile_;         ///< Timecorr/dipolar output
    static const char* Plegend_[]; ///< <P0>, <P1>, <P2>
    CorrF_FFT pubfft_;
    CorrF_Direct corfdir_;
};
#endif
