#ifndef INC_ANALYSIS_CLUSTERING_H
#define INC_ANALYSIS_CLUSTERING_H
#include "Analysis.h"
#include "ClusterList.h"
#include "TrajectoryFile.h"
#include "DataSet_Coords.h"
// Class: Analysis_Clustering
/// Used to perform clustering of frames, currently by RMSD only.
class Analysis_Clustering: public Analysis {
  public:
    Analysis_Clustering();
    ~Analysis_Clustering();

    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Clustering(); }
    static void Help();
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    DataSetList* masterDSL_;    ///< For Cluster pop v time DataSets.
    DataSet_Coords* coords_;    ///< Hold coordinates of frames being clustered.
    ClusterList* CList_;        ///< Hold specified clustering algorithm.
    std::string maskexpr_;      ///< If RMSD, Atoms to cluster on
    int sieve_;                 ///< If > 1, frames to skip on initial clustering pass.
    int splitFrame_;            ///< Frame to split at when comparing 1st to 2nd half.
    DataSet* cnumvtime_;        ///< Cluster vs time dataset.
    DataFile* cpopvtimefile_;   ///< Cluster pop v time file.
    std::string summaryfile_;   ///< Summary file name
    std::string halffile_;      ///< 1st/2nd half summary file name
    std::string clusterfile_;   ///< Cluster trajectory base filename.
    std::string singlerepfile_; ///< Cluster all rep single trajectory filename.
    std::string reptrajfile_;   ///< Cluster rep to separate trajectory filename.
    std::string clusterinfo_;   ///< Name for Ptraj-like cluster output file.
    std::string pairdistfile_;  ///< Name of pairwise-distances file.
    bool nofitrms_;             ///< If true do not best-fit when calc RMSD.
    bool usedme_;
    bool useMass_;
    bool grace_color_;          ///< If true print grace colors instead of cluster number
    bool norm_pop_;             ///< If true cluster pops v time will be normalized to 1.0
    bool load_pair_;            ///< If true, previously calcd pair dist file will be used if found
    ClusterDist::DsArray cluster_dataset_;  ///< DataSet(s) to use for clustering.
    /// Cluster trajectory format.
    TrajectoryFile::TrajFormatType clusterfmt_;
    /// Cluster all rep single trajectory format.
    TrajectoryFile::TrajFormatType singlerepfmt_;
    /// Cluster rep to separate trajectory format.
    TrajectoryFile::TrajFormatType reptrajfmt_;
    int debug_;
    static const char* PAIRDISTFILE;

    void CreateCnumvtime( ClusterList const& );
    void CreateCpopvtime( ClusterList const& );
    void WriteClusterTraj( ClusterList const& );
    void WriteSingleRepTraj( ClusterList const& );
    void WriteRepTraj( ClusterList const& );
};
#endif
