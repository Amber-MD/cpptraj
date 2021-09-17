#ifndef INC_CLUSTER_DRAWGRAPH_H
#define INC_CLUSTER_DRAWGRAPH_H
class DataSet;
namespace Cpptraj {
namespace Cluster {
class Cframes;
class MetricArray; 

enum GraphType { NO_DRAWGRAPH = 0, TWOD, THREED };

void DrawGraph(Cframes const&, MetricArray&, GraphType, DataSet*, double, int, int);

}
}
#endif
