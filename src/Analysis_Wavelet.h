#ifndef INC_ANALYSIS_WAVELET_H
#define INC_ANALYSIS_WAVELET_H
#include "Analysis.h"
#include "ComplexArray.h"
#include "DataSet_MatrixFlt.h"
#include "ClusterMap.h"
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
    Frame currentFrame_;
    DataSet_Coords* coords_;
    DataSet* output_;
    double S0_;
    double ds_;
    double correction_;
    double chival_;
    WaveletType wavelet_type_;
    int nb_;

    // Wavelet map clustering --------------------
    int WAFEX(DataSet_MatrixFlt const&);
    void ComputeKdist(int, DataSet_2D const&) const;

#   ifdef _OPENMP
    //std::vector<Iarray> thread_neighbors_; ///< RegionQuery neighbors for each thread.
    int numthreads_;                       ///< Total number of OpenMP threads.
#   endif
    ClusterMap CMAP_;
    DataSet* clustermap_; ///< Output cluster map
    DataSet* c_points_;
    DataSet* c_minatm_;
    DataSet* c_maxatm_;
    DataSet* c_minfrm_;
    DataSet* c_maxfrm_;
    DataSet* c_avgval_;
    std::string cprefix_; ///< Output cluster traj prefix
    std::string overlayName_;
    std::string overlayParm_;
    bool doClustering_; ///< Perform clustering on wavelet map
    bool cmap_square_;  ///< If true write cluster map by min/max rows/cols.
    bool doKdist_;      ///< If true calculate Kdist plot
#   ifdef TIMER
    Timer t_overall_;
#   endif
};
#endif
