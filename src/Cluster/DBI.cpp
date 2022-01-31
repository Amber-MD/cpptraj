#include <list>
#include <string>
#include <vector>
#include "DBI.h"
#include "List.h"
#include "Node.h"
#include "MetricArray.h"

/** The Davies-Bouldin Index (DBI) measures the average similarity between each
  * cluster and its most similar one; the smaller the DBI, the better. The DBI 
  * is defined as the average, for all clusters X, of fred, where fred(X) = max,
  * across other clusters Y, of (Cx + Cy)/dXY. Here Cx is the average distance
  * from points in X to the centroid, similarly Cy, and dXY is the distance 
  * between cluster centroids.
  * NOTE: To use this, cluster centroids should be fully up-to-date.
  */
double Cpptraj::Cluster::ComputeDBI(List const& clusters, std::vector<double>& averageDist, MetricArray& metricIn)
{
  averageDist.clear();
  averageDist.reserve( clusters.Nclusters() );
  for (List::cluster_iterator C1 = clusters.begincluster(); C1 != clusters.endcluster(); ++C1) {
    // Calculate average distance to centroid for this cluster
    averageDist.push_back( C1->CalcAvgToCentroid( metricIn ) );
    //if (outfile.IsOpen())
    //  outfile.Printf("#Cluster %i has average-distance-to-centroid %f\n", 
    //                 C1->Num(), averageDist.back());
  }
  double DBITotal = 0.0;
  unsigned int nc1 = 0;
  for (List::cluster_iterator c1 = clusters.begincluster(); c1 != clusters.endcluster(); ++c1, ++nc1) {
    double MaxFred = 0;
    unsigned int nc2 = 0;
    for (List::cluster_iterator c2 = clusters.begincluster(); c2 != clusters.endcluster(); ++c2, ++nc2) {
      if (c1 != c2) {
        double Fred = averageDist[nc1] + averageDist[nc2];
        Fred /= metricIn.CentroidDist( c1->Cent(), c2->Cent() );
        if (Fred > MaxFred)
          MaxFred = Fred;
      }
    }
    DBITotal += MaxFred;
  }
  DBITotal /= (double)clusters.Nclusters();
  //if (outfile.IsOpen()) outfile.Printf("#DBI: %f\n", DBITotal);
  return DBITotal;
}
