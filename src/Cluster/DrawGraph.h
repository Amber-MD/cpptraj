#ifndef INC_CLUSTER_DRAWGRAPH_H
#define INC_CLUSTER_DRAWGRAPH_H
class DataSet;
namespace Cpptraj {
namespace Cluster {
class Cframes;
class PairwiseMatrix;
void DrawGraph(Cframes const&, PairwiseMatrix const&, bool, DataSet*, double, int, int);

}
}
#endif
