#ifndef INC_ANALYSIS_TIMECORR_H
#define INC_ANALYSIS_TIMECORR_H
#include "Analysis.h"
#include "VectorType.h"
/** \author Original Code by Alrun N. Koller & H. Gohlke
  * \author Adapted by DRR
  */
class Analysis_Timecorr : public Analysis {
  public:
    Analysis_Timecorr();
    ~Analysis_Timecorr();

    int Setup(DataSetList*);
    int Analyze();
    //void Print(DataFileList*);
  private:
    enum timecorrMode { M_UNKNOWN, AUTO, CROSS };
    enum timecorrType { T_UNKNOWN, IRED, NORMAL };
    //static const double factor_;
    // 4/5*PI due to spherical harmonics addition theorem
    // FOURFIFTHSPI in Constants.h

    static const char ModeString[][6];
    timecorrMode mode_;
    timecorrType type_;
    bool dplr_;
    bool norm_;
    bool drct_;
    bool relax_;
    double tstep_;
    double tcorr_;
    double distnh_;
    double freq_;
    std::string filename_;
    std::string noeFilename_;
    int ndata_;
    double* table_;
    double* data1_;
    double* data2_;
    double* cf_;
    double* cfinf_;
    double* p2cf_;
    double* rcf_;
    // IRED only
    double* cf_cjt_;
    // IRED relax only
    double* taum_;
    
    VectorType* vinfo1_;
    VectorType* vinfo2_;
};
#endif
