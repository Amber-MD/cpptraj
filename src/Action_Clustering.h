#ifndef INC_ACTION_CLUSTERING_H
#define INC_ACTION_CLUSTERING_H
#include "Action.h"
#include "TriangleMatrix.h"
#include "ClusterList.h"
// Class: Clustering
/// Used to perform clustering of frames, currently by RMSD only.
class Clustering: public Action {
    ClusterList::LINKAGETYPE Linkage;
    FrameList ReferenceFrames; ///< Hold frames from all trajin stmts
    AtomMask Mask0;            ///< Target atom mask
    double epsilon;            ///< Once the min distance is > epsilon, stop clustering
    int targetNclusters;       ///< Once there are targetNclusters, stop clustering
    DataSet *cnumvtime;        ///< Cluster vs time dataset
    char *summaryfile;         ///< Summary file name
    char *halffile;            ///< 1st/2nd half summary file name
    char *clusterfile;         ///< Cluster trajectory base filename.
    FileFormat clusterfmt;     ///< Cluster trajectory format.
    char *singlerepfile;       ///< Cluster all rep single trajectory filename.
    FileFormat singlerepfmt;   ///< Cluster all rep single trajectory format.
    char *repfile;             ///< Cluster rep to separate trajectory filename.
    FileFormat repfmt;         ///< Cluster rep to separate trajectory format.
    char *clusterinfo;         ///< Ptraj-like cluster output
    bool nofitrms;             ///< If true do not best-fit when calc RMSD.
    bool grace_color;          ///< If true print grace colors instead of cluster number
    bool load_pair;            ///< If true, previously calcd pair dist file will be used if found
    DataSet *cluster_dataset;
    static const char PAIRDISTFILE[];

    int calcDistFromRmsd( TriangleMatrix *);
    int ClusterHierAgglo( TriangleMatrix *, ClusterList*);
    void CreateCnumvtime( ClusterList * );
    void WriteClusterTraj( ClusterList * );
    void WriteSingleRepTraj( ClusterList * );
    void WriteRepTraj( ClusterList * );
    void calcDistFromDataset( TriangleMatrix & );
  public:
    Clustering();

    int init();
    int action();
    void print();
};
#endif
