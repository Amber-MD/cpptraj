#ifndef INC_CLUSTER_DBI_H
#define INC_CLUSTER_DBI_H
namespace Cpptraj {
namespace Cluster {
class List;
class MetricArray;
/// Calculate the Davies-Bouldin index and average distance to centroid for each cluster
double ComputeDBI(List const&, std::vector<double>&, MetricArray&);
}
}
#endif
