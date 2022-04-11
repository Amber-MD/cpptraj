#include <list>
#include "PseudoF.h"
#include "List.h"
#include "MetricArray.h"
#include "Node.h"
#include "../Constants.h"
#include "../CpptrajStdio.h"

/** The pseudo-F statistic is another measure of clustering goodness. It is 
  * intended to capture the 'tightness' of clusters, and is in essence a ratio
  * of the mean sum of squares between groups to the mean sum of squares within
  * group (Lattin et al., 2003: 291) High values are good. Generally, one 
  * selects a cluster-count that gives a peak in the pseudo-f statistic (or 
  * pSF, for short).
  * Formula: A/B, where A = (T - P)/(G-1), and B = P / (n-G). Here n is the 
  * number of points, G is the number of clusters, T is the total distance from
  * the all-data centroid, and P is the sum (for all clusters) of the distances
  * from the cluster centroid.
  * This will also calculate SSR (sum of squares regression, between sum of
  * squares) over the total sum of squares (SST). The SSR/SST ratio value lies
  * between 0 and 1 and gives the percentage of explained variance by the data,
  * and is similar to the R 2 value in regression analysis. As the ratio
  * inherently rises with cluster count, one looks for an “elbow” in the curve
  * where adding another cluster does not add much new information.
  * NOTE: To use this, cluster centroids should be fully up-to-date.
  * NOTE: This calc differs slightly from PTRAJ in that real centroids are used
  *       instead of representative structures.
  * \param clusters List of clusters
  * \param SSRSST sum of squares regression (SSR or between sum of squares) over total sum of squares (SST).
  */
double Cpptraj::Cluster::ComputePseudoF(List const& clusters, double& SSRSST, MetricArray& metricIn, int debugIn)
{
  // Calculation makes no sense with fewer than 2 clusters.
  if (clusters.Nclusters() < 2) {
    mprintf("Warning: Fewer than 2 clusters. Not calculating pseudo-F.\n");
    return 0.0;
  }

  // Form a cluster with all points to get a centroid. Use only frames that
  // are in clusters, i.e. ignore noise. Assumes all cluster centroids are
  // up to date.
  Node c_all;
  for (List::cluster_iterator C1 = clusters.begincluster(); C1 != clusters.endcluster(); ++C1)
  {
    for (Node::frame_iterator f1 = C1->beginframe(); f1 != C1->endframe(); ++f1)
      c_all.AddFrameToCluster( *f1 );
  }
  // Pseudo-F makes no sense if # clusters == # frames
  if (clusters.Nclusters() == c_all.Nframes()) {
    mprintf("Warning: Each frame is in a separate cluster. Not calculating pseudo-F.\n");
    return 0.0;
  }
  c_all.SortFrameList();
  c_all.CalculateCentroid( metricIn );

  // Loop over all clusters
  double gss = 0.0; // between-group sum of squares
  double wss = 0.0; // within-group sum of squares
  for (List::cluster_iterator C1 = clusters.begincluster(); C1 != clusters.endcluster(); ++C1)
  {
    for (Node::frame_iterator f1 = C1->beginframe(); f1 != C1->endframe(); ++f1)
    {
      double dist = metricIn.FrameCentroidDist(*f1, c_all.Cent());
      gss += (dist * dist);
      dist = metricIn.FrameCentroidDist(*f1, C1->Cent());
      wss += (dist * dist);
    }
  }
  double d_nclusters = (double)clusters.Nclusters();
  double d_ntotal = (double)c_all.Nframes();
  double num = (gss - wss) / (d_nclusters - 1.0);
  double den = wss / (d_ntotal - d_nclusters);
  if (den < Constants::SMALL)
    den = Constants::SMALL;
  double pseudof = num / den;
  if (debugIn > 0)
    mprintf("Pseudo-f: Total distance to centroid is %.4f\n"
            "Pseudo-f: Cluster distance to centroid is %.4f\n"
            "Pseudo-f: Numerator %.4f over denominator %.4f gives %.4f\n", 
            gss, wss, num, den, pseudof);
  // This calculation taken directly from ptraj
  SSRSST = pseudof*(d_nclusters-1)/(d_ntotal-d_nclusters+pseudof*(d_nclusters-1));

  return pseudof;
}
