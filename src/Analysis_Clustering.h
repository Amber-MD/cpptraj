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

    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Clustering(); }
    void Help() const;
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    inline void GetClusterTrajArgs(ArgList&, const char*, const char*, std::string&,
                                   TrajectoryFile::TrajFormatType&) const;
    void CreateCnumvtime( ClusterList const&, unsigned int );
    void CreateCpopvtime( ClusterList const&, unsigned int );
    void ClusterLifetimes( ClusterList const&, unsigned int );
    void NclustersObserved(ClusterList const&, unsigned int);
    void WriteClusterTraj( ClusterList const& );
    void WriteAvgStruct( ClusterList const& );
    void WriteSingleRepTraj( ClusterList const& );
    void WriteRepTraj( ClusterList const& );

    DataSetList* masterDSL_;    ///< For Cluster pop v time DataSets.
    DataSet_Coords* coords_;    ///< Hold coordinates of frames being clustered.
    ClusterList* CList_;        ///< Hold specified clustering algorithm.
    std::string maskexpr_;      ///< If RMSD, Atoms to cluster on
    int sieve_;                 ///< If > 1, frames to skip on initial clustering pass.
    int sieveSeed_;             ///< Used to seed random number gen for sieve
    int windowSize_;            ///< Window size for # clusters seen vs time.
    int drawGraph_;
    int draw_maxit_;
    double draw_tol_;
    std::vector<int> splitFrames_; ///< Frames to split at when comparing parts.
    DataSet* cnumvtime_;        ///< Cluster vs time dataset.
    DataSet* clustersVtime_;    ///< # clusters seen vs time dataset.
    DataSet* pw_dist_;          ///< Cluster pairwise distance matrix dataset
    DataFile* cpopvtimefile_;   ///< Cluster pop v time file.
    DataFile* pwd_file_;        ///< Data file to write pairwise distance matrix to.
    std::string summaryfile_;   ///< Summary file name
    std::string halffile_;      ///< 1st/2nd half summary file name
    std::string clusterfile_;   ///< Cluster trajectory base filename.
    std::string singlerepfile_; ///< Cluster all rep single trajectory filename.
    std::string reptrajfile_;   ///< Cluster rep to separate trajectory filename.
    std::string avgfile_;       ///< Cluster traj average structure filename.
    std::string clusterinfo_;   ///< Name for Ptraj-like cluster output file.
    std::string sil_file_;      ///< Prefix name of file for cluster silhouette.
    bool nofitrms_;             ///< If true do not best-fit when calc RMSD.
    ClusterList::DistMetricType metric_;
    bool useMass_;
    bool grace_color_;          ///< If true print grace colors instead of cluster number
    enum normPopType { NONE=0, CLUSTERPOP, FRAME };
    normPopType norm_pop_;      ///< If set cluster pops v time will be normalized 
    bool calc_lifetimes_;       ///< If true calculate DataSets for use in lifetime analysis.
    bool writeRepFrameNum_;     ///< If true frame #s will be in rep file names.
    bool suppressInfo_;         ///< If true, do not print cluster info to STDOUT
    ClusterDist::DsArray cluster_dataset_;        ///< DataSet(s) to use for clustering.
    TrajectoryFile::TrajFormatType clusterfmt_;   ///< Cluster trajectory format.
    TrajectoryFile::TrajFormatType singlerepfmt_; ///< Cluster all rep single trajectory format.
    TrajectoryFile::TrajFormatType reptrajfmt_;   ///< Cluster rep to separate trajectory format.
    TrajectoryFile::TrajFormatType avgfmt_;       ///< Cluster traj average structure file format.
    static const TrajectoryFile::TrajFormatType DEF_TRAJ_FMT_;
    int debug_;
    static const char* PAIRDISTFILE_;              ///< Default pairwise dist file name.
    static DataFile::DataFormatType PAIRDISTTYPE_; ///< Default pairwise dist file type.
};
#endif
