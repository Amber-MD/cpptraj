#ifndef INC_ANALYSIS_IRED_H
#define INC_ANALYSIS_IRED_H
#include "Analysis.h"
#include "DataSet_Vector.h"
#include "DataSet_Modes.h"
#include "PubFFT.h"
/** \author Original Code by Alrun N. Koller & H. Gohlke
  * \author Adapted by DRR
  */
class Analysis_IRED : public Analysis {
  public:
    Analysis_IRED();
    ~Analysis_IRED();

    int Setup(DataSetList*);
    int Analyze();
    //void Print(DataFileList*);
  private:
    //static const double factor_;
    // 4/5*PI due to spherical harmonics addition theorem
    // FOURFIFTHSPI in Constants.h
    double freq_;
    double tstep_;
    double tcorr_;
    double distnh_;
    int order_;
    bool relax_;
    bool norm_;
    bool drct_;
    double* data1_;
    double* table_;
    double* cf_;
    double* cf_cjt_;
    double* cfinf_;
    double* taum_;
    std::string noeFilename_;
    std::string filename_;
    std::string orderparamfile_;
    PubFFT pubfft_;
    DataSet_Modes* modinfo_;
    std::vector<DataSet_Vector*> IredVectors_;

    double calc_spectral_density(int, double);
};
#endif
