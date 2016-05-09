#ifndef INC_ANALYSIS_WAVELET_H
#define INC_ANALYSIS_WAVELET_H
#include "Analysis.h"
#include "ComplexArray.h"
#include "DataSet_MatrixFlt.h"
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

    // Wavelet map clustering --------------------
    typedef std::vector<int> Iarray;

    int ClusterMap(DataSet_MatrixFlt const&);
    void RegionQuery(Iarray&, double, int, DataSet_2D const&);
    void AddCluster(Iarray const&, DataSet_2D const&);

    class Cluster;
    typedef std::vector<Cluster> Carray;
    Carray clusters_;                      ///< Hold all clusters.
#   ifdef _OPENMP
    std::vector<Iarray> thread_neighbors_; ///< RegionQuery neighbors for each thread.
    int mythread_;                         ///< Current OpenMP thread.
#   endif
    DataSet* clustermap_; ///< Output cluster map
    DataSet* c_points_;
    DataSet* c_minatm_;
    DataSet* c_maxatm_;
    DataSet* c_minfrm_;
    DataSet* c_maxfrm_;
    DataSet* c_avgval_;
    double epsilon_;  ///< Distance to search for neighboring points within.
    double epsilon2_; ///< Epsilon squared.
    double Avg_;      ///< Average value of map, used as cutoff (points below are noise).
    int minPoints_;   ///< Minimum number of points within epsilon to qualify as cluster.
    int nClusters_;   ///< Current number of clusters.
    bool doClustering_; ///< Perform clustering on wavelet map
#   ifdef TIMER
    Timer t_overall_;
    Timer t_query1_;
    Timer t_query2_;
#   endif
};
// -----------------------------------------------
class Analysis_Wavelet::Cluster {
  public:
    Cluster() : avg_(0.0), cnum_(-1), min_col_(-1), max_col_(-1), min_row_(-1), max_row_(-1) {}
    Cluster(Iarray const& p, double A, int C, int minc, int maxc, int minr, int maxr) :
      points_(p), avg_(A), cnum_(C), min_col_(minc), max_col_(maxc), min_row_(minr), max_row_(maxr)
      {}
    void SetCnum(int c) { cnum_ = c; }
    Iarray const& Points() const { return points_;  }
    double Avg()           const { return avg_;     }
    int Cnum()             const { return cnum_;    }
    int MinCol()           const { return min_col_; }
    int MaxCol()           const { return max_col_; }
    int MinRow()           const { return min_row_; }
    int MaxRow()           const { return max_row_; }
    // Use > since we give higher priority to larger clusters
    bool operator <(Cluster const& rhs) const {
      return ( points_.size() > rhs.points_.size() );
    }
  private:
    Iarray points_;
    double avg_;
    int cnum_;
    int min_col_;
    int max_col_;
    int min_row_;
    int max_row_;
};
#endif
