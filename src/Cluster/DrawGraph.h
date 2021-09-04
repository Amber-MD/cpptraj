#ifndef INC_CLUSTER_DRAWGRAPH_H
#define INC_CLUSTER_DRAWGRAPH_H
class DataSet;
namespace Cpptraj {
namespace Cluster {
class Cframes;
class PairwiseMatrix;

enum GraphType { TWOD = 0, THREED };

void DrawGraph(Cframes const&, PairwiseMatrix const&, GraphType, DataSet*, double, int, int);

}
}
#endif
