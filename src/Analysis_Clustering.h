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

    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Clustering(); }
    static void Help();
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,int);
    Analysis::RetType Analyze();
    void Print(DataFileList*);
  private:
    DataSet_Coords* coords_;
    std::string maskexpr_;      ///< If RMSD, Atoms to cluster on
    double epsilon_;            ///< Once the min distance is > epsilon, stop clustering
    int targetNclusters_;       ///< Once there are targetNclusters, stop clustering
    int sieve_;
    DataSet* cnumvtime_;        ///< Cluster vs time dataset.
    std::string cnumvtimefile_; ///< Cluster vs time filename.
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
    bool load_pair_;            ///< If true, previously calcd pair dist file will be used if found
    ClusterDist::DsArray cluster_dataset_;  ///< DataSet(s) to use for clustering.
    /// Cluster linkage type
    ClusterList::LINKAGETYPE Linkage_;
    /// Clustering algorithm to use.
    ClusterList::ClusterAlgorithm mode_;
    /// Cluster trajectory format.
    TrajectoryFile::TrajFormatType clusterfmt_;
    /// Cluster all rep single trajectory format.
    TrajectoryFile::TrajFormatType singlerepfmt_;
    /// Cluster rep to separate trajectory format.
    TrajectoryFile::TrajFormatType reptrajfmt_;
    int debug_;
    static const char* PAIRDISTFILE; // TODO: Make this a user option

    void CreateCnumvtime( ClusterList & );
    void WriteClusterTraj( ClusterList & );
    void WriteSingleRepTraj( ClusterList & );
    void WriteRepTraj( ClusterList & );
};
#endif
