#ifndef INC_ANALYSIS_TIMECORR_H
#define INC_ANALYSIS_TIMECORR_H
#include "Analysis.h"
#include "DataSet_Vector.h"
#include "PubFFT.h"
/** \author Original Code by Alrun N. Koller & H. Gohlke
  * \author Adapted by DRR
  */
class Analysis_Timecorr : public Analysis {
  public:
    Analysis_Timecorr();

    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Timecorr(); }
    static void Help();

    ~Analysis_Timecorr();

    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,int);
    Analysis::RetType Analyze();
    void Print(DataFileList*) { return; }
  private:
    enum timecorrMode { AUTO = 0, CROSS };
    //static const double factor_;
    // 4/5*PI due to spherical harmonics addition theorem
    // FOURFIFTHSPI in Constants.h
    static const char ModeString[][6];

    double tstep_;
    double tcorr_;
    int order_;
    timecorrMode mode_;
    bool dplr_;
    bool norm_;
    bool drct_;
    double* table_;
    double* data1_;
    double* data2_;
    double* cf_;
    double* p2cf_;
    double* rcf_;
    DataSet_Vector* vinfo1_;
    DataSet_Vector* vinfo2_;
    std::string filename_;
    PubFFT pubfft_;
    
    void CalcCorr(int,int,int);
};
#endif
