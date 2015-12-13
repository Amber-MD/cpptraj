#ifndef INC_ANALYSIS_IRED_H
#define INC_ANALYSIS_IRED_H
#include "Analysis.h"
#include "DataSet_Vector.h"
#include "DataSet_Modes.h"
/** \author Original Code by Alrun N. Koller & H. Gohlke
  * \author Adapted by DRR
  */
class Analysis_IRED : public Analysis {
  public:
    Analysis_IRED();

    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_IRED(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    double Jw(int, double, std::vector<double>) const;

    double freq_;           ///< Frequency for calculation of relaxation parameters.
    double tstep_;          ///< Time step
    double tcorr_;          ///< Total correlation time.
    double distnh_;         ///< N-H distance for relaxation calc in Ang.
    int order_;             ///< Order of spherical harmonics for calculation autocorrelation fns
    int debug_;
    bool relax_;            ///< If true calculate relaxation and NOEs
    bool norm_;             ///< If true output normalized time correlation functions
    bool drct_;             ///< If true use direct calculation of autocorrelation instead of FFT
    DataFile* cmtfile_;     ///< File to write mode autocorrelation functions to
    DataFile* cjtfile_;     ///< File to write reconstructed vector autocorrelation functions to
    std::string dsname_;    ///< Data set name
    DataSet* data_s2_;      ///< Order parameters (1 per vector).
    DataSet* data_plateau_; ///< Cm(t) plateau values, i.e. Cm(t->T), 1 per vector
    DataSet* data_tauM_;    ///< Cm(t) relaxation values, 1 per vector.
    DataSet* data_noe_;     ///< NOEs (1 per vector)
    DataSet* data_t1_;      ///< T1 relaxation (1 per vector).
    DataSet* data_t2_;      ///< T2 relaxation (1 per vector).
    DataSet* data_ds2_mat_; ///< delta * S^2 matrix(j,m)
    typedef std::vector<DataSet*> DataListType;
    DataListType CmtArray_; ///< Cm(t) for each mode
    DataListType CjtArray_; ///< Cj(t) for each vector
    DataSetList* masterDSL_; 
    
    DataSet_Modes* modinfo_; ///< Modes data from prior diagonalization of IRED matrix
    std::vector<DataSet_Vector*> IredVectors_; ///< IRED vectors used to calc IRED matrix
};
#endif
