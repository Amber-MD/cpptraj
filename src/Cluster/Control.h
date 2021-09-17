#ifndef INC_CLUSTER_CONTROL_H
#define INC_CLUSTER_CONTROL_H
#include "Algorithm.h" // Algorithm::AType
#include "BestReps.h"  // BestReps::RepMethodType
#include "DrawGraph.h" // GraphType
#include "List.h"
#include "MetricArray.h"
#include "Node.h" // Node::CnormType
#include "Sieve.h"
#include "../Timer.h"
class DataSet_Coords;
class DataSet_PairwiseCache;
class DataFile;
class DataFileList;
class DataSetList;
namespace Cpptraj {
namespace Cluster {
class Algorithm;
class Results;
/// Hold clusters, algorithm, and pairwise matrix.
class Control {
  public:
    Control();
    ~Control();

    /// Set debug level
    void SetDebug(int d) { debug_ = d; }

    /// For determining how any sieved frames should be restored.
    enum SieveRestoreType { NO_RESTORE = 0, CLOSEST_CENTROID, EPSILON_CENTROID, EPSILON_FRAME };
    /// For determining how frames to cluster will be determined.
    enum FrameSelectType { UNSPECIFIED = 0, FROM_CACHE };

    /// Setup clustering
    int SetupClustering(DataSetList const&, DataSet*, ArgList&, DataSetList&, DataFileList&, int);
    /// Provide information on how clustering calculation is currently set up.
    void Info() const;
    /// Perform clustering.
    int Run();
    /// Do any clustering output. TODO make const?
    int Output(DataSetList&);
    /// Print timing data
    void Timing(double) const;
    /// Print help text to STDOUT
    static void Help();
  private:
    // Help keywords
    static const char* AlgorithmArgs_;
    static const char* CoordsDataSetArgs_;
    static const char* SieveArgs1_;
    static const char* SieveArgs2_;
    static const char* BestRepArgs_;
    static const char* OutputArgs1_;
    static const char* OutputArgs2_;
    static const char* OutputArgs3_;
    static const char* OutputArgs4_;
    static const char* GraphArgs_;

    /// \return Algorithm of given type
    static Algorithm* AllocateAlgorithm(Algorithm::AType);
    /// Allocate algorithm from keywords
    int AllocateAlgorithm(ArgList&);    
    /// Initialize clusters from Info file
    int ReadInfo(std::string const&);
    /// Initialize clusters from cluster number vs time DataSet
    int InitClustersFromSet(DataSet*);

    List clusters_;                ///< Hold cluster results.
    Sieve frameSieve_;             ///< Hold frames to cluster, frames to "sieve" out.
    MetricArray metrics_;          ///< Hold the distance metrics and any cached distances.
    Algorithm* algorithm_;         ///< Hold the clustering algorithm.
    Results* results_;             ///< Hold output routines specific to data being clustered.
    std::string dsname_;           ///< Name for output data sets.
    int verbose_;

    FrameSelectType frameSelect_;    ///< How frames to cluster should be determined.
    int sieve_;                      ///< Sieve value
    int sieveSeed_;                  ///< Seed if doing random sieve
    SieveRestoreType sieveRestore_;  ///< Specify if/how sieved frames should be restored.
    double restoreEpsilon_;          ///< Cutoff to use if restoring by epsilon to centroids.
    bool includeSieveInCalc_;        ///< If true include sieved frames in later calculations.
    bool includeSieveCdist_;         ///< If true include sieved frames in cluster distance calc.

    BestReps::RepMethodType bestRep_; ///< How to determine best rep frames.
    int nRepsToSave_;                 ///< How many rep frames to save.

    bool suppressInfo_;               ///< If true do not write cluster info
    std::string clusterinfo_;         ///< Cluster info file name.
    std::string summaryfile_;         ///< Cluster summary file name.
    std::string sil_file_;            ///< File prefix for writing silhouette data

    DataSet* cnumvtime_;              ///< Cluster number vs time data set.
    bool grace_color_;                ///< If true change cluster number to grace color

    DataSet* clustersVtime_;          ///< Number of unique clusters vs time
    int windowSize_;                  ///< Window size for determining number unique clusters v time

    DataFile* cpopvtimefile_;         ///< Cluster population vs time file.
    Node::CnormType norm_pop_;        ///< Cluster pop vs time normalization type

    bool calc_lifetimes_;             ///< If true create cluster lifetime data sets

    std::string splitfile_;           ///< Output file for splitting cluster results
    Cframes splitFrames_;             ///< Frames at which to split
    // TODO make class vars in DrawGraph.cpp
    GraphType drawGraph_;             ///< Indicate whether a cluster graph should be drawn
    double draw_tol_;                 ///< Graph draw tolerance for min
    int draw_maxit_;                  ///< Graph draw max iterations for min

    int debug_;                       ///< Cluster debug level

    CpptrajFile* metricContribFile_;  ///< If not null, determine how much each metric contributes to total distance

    // Timers
    Timer timer_setup_;          ///< Run - metric, frames to cluster setup 
    Timer timer_pairwise_;       ///< Run - pairwise caching
    Timer timer_cluster_;        ///< Run - clustering
    Timer timer_post_;           ///< Run - post-clustering calcs (rep frames etc)
    Timer timer_post_renumber_;  ///< Run Post - cluster sort/renumber
    Timer timer_post_bestrep_;   ///< Run Post - best rep frame calc
    Timer timer_run_;            ///< Total Run time
    Timer timer_output_info_;    ///< Output - info file write
    Timer timer_output_summary_; ///< Output - summary write
    Timer timer_output_results_; ///< Output - results (e.g. coords writes).
    Timer timer_output_;         ///< Total output time
};

} /** END namespace Cluster. */
} /** END namespace Cpptraj. */
#endif
