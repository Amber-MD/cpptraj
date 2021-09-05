#ifndef INC_CLUSTER_ASSIGNREFSTOCLUSTERS_H
#define INC_CLUSTER_ASSIGNREFSTOCLUSTERS_H
#include <string>
class DataSet_Coords;
class DataSetList;
namespace Cpptraj {
namespace Cluster {
class List;

int AssignRefsToClusters(DataSetList const&, std::string const&, double, bool, 
                         DataSet_Coords&, List& CList);

}
}
#endif
