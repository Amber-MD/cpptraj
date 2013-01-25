#ifndef INC_ANALYSIS_TIMECORR_H
#define INC_ANALYSIS_TIMECORR_H
#include "Analysis.h"
#include "DataSet_Vector.h"
#include "ComplexArray.h"
/** \author Original Code by Alrun N. Koller & H. Gohlke
  * \author Adapted by DRR
  */
class Analysis_Timecorr : public Analysis {
  public:
    Analysis_Timecorr();

    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Timecorr(); }
    static void Help();

    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    struct AvgResults {
      double avgr_;
      double rave_;
      double r3iave_;
      double r6iave_;
    };
    enum timecorrMode { AUTO = 0, CROSS };
    static const char* ModeString[];

    double tstep_;
    double tcorr_;
    int order_;
    timecorrMode mode_;
    bool dplr_;
    bool norm_;
    bool drct_;
    ComplexArray data1_;
    ComplexArray data2_;
    DataSet_Vector* vinfo1_;
    DataSet_Vector* vinfo2_;
    std::string filename_;
    CorrF_FFT pubfft_;
    CorrF_Direct corfdir_;
    
    std::vector<double> CalculateAverages(DataSet_Vector&, AvgResults&);
    void CalcCorr(int);
};
#endif
