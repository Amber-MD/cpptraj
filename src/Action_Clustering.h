#ifndef INC_ACTION_CLUSTERING_H
#define INC_ACTION_CLUSTERING_H
/// Class: Clustering
#include "Action.h"
#include "TriangleMatrix.h"
#include "ClusterList.h"
class Clustering: public Action {
    enum LINKAGETYPE {SINGLELINK, AVERAGELINK};
    LINKAGETYPE Linkage;
    FrameList ReferenceFrames; // Hold frames from all trajin stmts
    AtomMask Mask0;        // Target atom mask
    double epsilon; // Once the min distance is > epsilon, stop clustering
    int targetNclusters; // Once there are targetNclusters, stop clustering
    DataSet *cnumvtime;
    char *summaryfile; // SUmmary file name

    int calcDistFromRmsd( TriangleMatrix *);
    int ClusterHierAgglo( TriangleMatrix *, ClusterList*);
  public:
    Clustering();
    ~Clustering();

    int init();
    int setup();
    int action();
    void print();
};
#endif
