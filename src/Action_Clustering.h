#ifndef INC_ACTION_CLUSTERING_H
#define INC_ACTION_CLUSTERING_H
/// Class: Clustering
#include "Action.h"
#include "TriangleMatrix.h"
class Clustering: public Action {
    enum LINKAGETYPE {SINGLELINK, AVERAGELINK};
    LINKAGETYPE Linkage;
    FrameList ReferenceFrames; // Hold frames from all trajin stmts
    AtomMask Mask0;        // Target atom mask
    double epsilon; // Once the min distance is > epsilon, stop clustering
    int targetNclusters; // Once there are targetNclusters, stop clustering

    int calcDistFromRmsd( TriangleMatrix *);
    int ClusterHierAgglo( TriangleMatrix *);
  public:
    Clustering();
    ~Clustering();

    int init();
    int setup();
    int action();
    void print();
};
#endif
