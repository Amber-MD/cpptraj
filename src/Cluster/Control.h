#ifndef INC_CLUSTER_CONTROL_H
#define INC_CLUSTER_CONTROL_H
#include "List.h"
#include "PairwiseMatrix.h"
#include "Algorithm.h"
#include "Metric.h"
#include "BestReps.h"
#include "../Timer.h"
#include "../DataSetList.h"
#include "../DataFileList.h"
#include "../DataSet_Coords.h"
#include "../DataSet_PairwiseCache.h"
namespace Cpptraj {
namespace Cluster {

/// Hold clusters, algorithm, and pairwise matrix.
class Control {
  public:
    Control();
    ~Control();

    static const char* PairwiseArgs;
    static const char* AlgorithmArgs;

    enum SieveRestoreType { NO_RESTORE = 0, CLOSEST_CENTROID, EPSILON_CENTROID, EPSILON_FRAME };

    int SetupForCoordsDataSet(DataSet_Coords*, std::string const&, ArgList&, DataSetList&,
                              DataFileList&, int);
    /// Provide information on how clustering calculation is currently set up.
    void Info() const;
    /// Perform clustering.
    int Run();
    /// Do any clustering output. TODO
    //int Output();
    /// Print timing data
    void Timing(double) const;

    List const& Clusters()     const { return clusters_; }
    Metric const& DistMetric() const { return *metric_;  }
  private:
    int AllocatePairwise(ArgList&, DataSetList&, DataFileList&);

    static Metric* AllocateMetric(Metric::Type);

    static Algorithm* AllocateAlgorithm(Algorithm::AType);
    int AllocateAlgorithm(ArgList&);    

    int Common(ArgList&, DataSetList&, DataFileList&);

    static const char* DEFAULT_PAIRDIST_NAME_;

    static DataFile::DataFormatType DEFAULT_PAIRDIST_TYPE_;

    List clusters_;                ///< Hold cluster results.
    Metric* metric_;               ///< Hold the distance metric.
    DataSet_PairwiseCache* cache_; ///< Hold any cached pairwise distances.
    PairwiseMatrix pmatrix_;       ///< Encapsulates the metric and any cached distances.
    Algorithm* algorithm_;         ///< Hold the clustering algorithm.
    int verbose_;

    int sieve_;                     ///< Sieve value
    int sieveSeed_;                 ///< Seed if doing random sieve
    SieveRestoreType sieveRestore_; ///< Specify if/how sieved frames should be restored.
    double restoreEpsilon_;         ///< Cutoff to use if restoring by epsilon to centroids.
    bool includeSieveInCalc_;       ///< If true include sieved frames in later calculations.

    BestReps::RepMethodType bestRep_; ///< How to determine best rep frames.
    int nRepsToSave_;                 ///< How many rep frames to save.

    bool suppressInfo_;               ///< If true do not write cluster info
    std::string clusterinfo_;         ///< Cluster info file name.
    std::string summaryfile_;         ///< Cluster summary file name.
    std::string sil_file_;            ///< File prefix for writing silhouette data

    DataFile* cpopvtimefile_;
    Node::CnormType norm_pop_;

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
    Timer timer_output_;         ///< Total output time
};

} /** END namespace Cluster. */
} /** END namespace Cpptraj. */
#endif
