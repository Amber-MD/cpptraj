#ifndef INC_ACTION_CLUSTERING_H
#define INC_ACTION_CLUSTERING_H
/// Class: Clustering
#include "Action.h"
#include "TriangleMatrix.h"
#include "ClusterList.h"
class Clustering: public Action {
    ClusterList::LINKAGETYPE Linkage;
    FrameList ReferenceFrames; // Hold frames from all trajin stmts
    AtomMask Mask0;            // Target atom mask
    double epsilon;            // Once the min distance is > epsilon, stop clustering
    int targetNclusters;       // Once there are targetNclusters, stop clustering
    DataSet *cnumvtime;        // Cluster vs time dataset
    char *summaryfile;         // Summary file name
    char *clusterfile;         // Cluster trajectory base filename.
    FileFormat clusterfmt;     // Cluster trajectory format.

    int calcDistFromRmsd( TriangleMatrix *);
    int ClusterHierAgglo( TriangleMatrix *, ClusterList*);
    void CreateCnumvtime( ClusterList * );
    void WriteClusterTraj( ClusterList * );
  public:
    Clustering();
    ~Clustering();

    int init();
    int action();
    void print();
};
#endif
