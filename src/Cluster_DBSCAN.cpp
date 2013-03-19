#include "Cluster_DBSCAN.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"

Cluster_DBSCAN::Cluster_DBSCAN() :
  minPoints_(-1),
  epsilon_(-1.0)
{}

void Cluster_DBSCAN::Help() {
  mprintf("\t[dbscan minpoints <n> epsilon <e>]\n");
}

int Cluster_DBSCAN::SetupCluster(ArgList& analyzeArgs) {
  minPoints_ = analyzeArgs.getKeyInt("minpoints", -1);
  if (minPoints_ < 1) {
    mprinterr("Error: DBSCAN requires minimum # of points to be set and >= 1\n");
    mprinterr("Error: Use 'minpoints <N>'\n");
    return 1;
  }
  epsilon_ = analyzeArgs.getKeyDouble("epsilon", -1.0);
  if (epsilon_ <= 0.0) {
    mprinterr("Error: DBSCAN requires epsilon to be set and > 0.0\n");
    mprinterr("Error: Use 'epsilon <e>'\n");
    return 1;
  }
  return 0;
}

void Cluster_DBSCAN::ClusteringInfo() {
  mprintf("\tDBSCAN:\n");
  mprintf("\t\tMinimum pts to form cluster= %i\n", minPoints_);
  mprintf("\t\tCluster distance criterion= %.3f\n", epsilon_);
}

void Cluster_DBSCAN::RegionQuery(std::vector<int>& NeighborPts,
                                 std::vector<int> const& FramesToCluster,
                                 int point)
{
  NeighborPts.clear();
  for (std::vector<int>::const_iterator otherpoint = FramesToCluster.begin();
                                        otherpoint != FramesToCluster.end();
                                        ++otherpoint)
  {
    if (point == *otherpoint) continue;
    if ( FrameDistances_.GetElement(point, *otherpoint) < epsilon_ )
      NeighborPts.push_back( *otherpoint );
  }
}

/** Ester, Kriegel, Sander, Xu; Proceedings of 2nd International Conference
  * on Knowledge Discovery and Data Mining (KDD-96); pp 226-231.
  */
int Cluster_DBSCAN::Cluster() {
  std::vector<int> NeighborPts;
  std::vector<int> Npts2; // Will hold neighbors of a neighbor
  std::vector<int> FramesToCluster;
  ClusterDist::Cframes cluster_frames;
  // First determine which frames are being clustered.
  for (int frame = 0; frame < (int)FrameDistances_.Nframes(); ++frame)
    if (!FrameDistances_.IgnoringRow( frame ))
      FramesToCluster.push_back( frame );
  // Set up array to keep track of points that have been visited.
  // Make it the size of FrameDistances so we can index into it. May
  // waste memory during sieving but makes code easier.
  std::vector<bool> Visited( FrameDistances_.Nframes(), false );
  // Set up array to keep track of whether points are noise or in a cluster.
  static char UNASSIGNED = 'U';
  static char NOISE = 'N';
  static char INCLUSTER = 'C';
  std::vector<char> Status( FrameDistances_.Nframes(), UNASSIGNED);
  mprintf("\tStarting DBSCAN Clustering:\n");
  ProgressBar cluster_progress(FramesToCluster.size());
  int iteration = 0;
  for (std::vector<int>::iterator point = FramesToCluster.begin();
                                  point != FramesToCluster.end(); ++point)
  {
    if (!Visited[*point]) {
      // Mark this point as visited
      Visited[*point] = true;
      // Determine how many other points are near this point
      RegionQuery( NeighborPts, FramesToCluster, *point );
      if (debug_ > 0) {
        mprintf("\tPoint %i\n", *point + 1);
        mprintf("\t\t%u neighbors:", NeighborPts.size());
      }
      // If # of neighbors less than cutoff, noise; otherwise cluster
      if ((int)NeighborPts.size() < minPoints_) {
        if (debug_ > 0) mprintf(" NOISE\n");
        Status[*point] = NOISE;
      } else {
        // Expand cluster
        cluster_frames.clear();
        cluster_frames.push_back( *point );
        // NOTE: Use index instead of iterator since NeighborPts may be
        //       modified inside this loop.
        unsigned int endidx = NeighborPts.size();
        for (unsigned int idx = 0; idx < endidx; ++idx) {
          int neighbor_pt = NeighborPts[idx];
          if (!Visited[neighbor_pt]) {
            if (debug_ > 0) mprintf(" %i", neighbor_pt + 1);
            // Mark this neighbor as visited
            Visited[neighbor_pt] = true;
            // Determine how many other points are near this neighbor
            RegionQuery( Npts2, FramesToCluster, neighbor_pt );
            if ((int)Npts2.size() >= minPoints_) {
              // Add other points to current neighbor list
              NeighborPts.insert( NeighborPts.end(), Npts2.begin(), Npts2.end() );
              endidx = NeighborPts.size();
            }
          }
          // If neighbor is not yet part of a cluster, add it to this one.
          if (Status[neighbor_pt] != INCLUSTER) {
            cluster_frames.push_back( neighbor_pt );
            Status[neighbor_pt] = INCLUSTER;
          }
        }
        // Remove duplicate frames
        // TODO: Take care of this in Renumber?
        cluster_frames.sort();
        cluster_frames.unique();
        // Add cluster to the list
        AddCluster( cluster_frames );
        if (debug_ > 0) {
          mprintf("\n");
          PrintClusters();
        }
      }
    }
    cluster_progress.Update(iteration++);
  } // END loop over FramesToCluster
  // Count the number of noise points
  mprintf("\tNOISE FRAMES:");
  unsigned int frame = 1;
  for (std::vector<char>::iterator stat = Status.begin();
                                   stat != Status.end(); ++stat, ++frame)
  {
    if ( *stat == NOISE )
      mprintf(" %i", frame);
  }
  mprintf("\n");
  // Calculate the distances between each cluster based on centroids
  ClusterDistances_.SetupMatrix( clusters_.size() );
  // Make sure centroid for clusters are up to date
  for (cluster_it C1 = clusters_.begin(); C1 != clusters_.end(); ++C1)
    (*C1).CalculateCentroid( Cdist_ );
  // Calculate distances between each cluster centroid
  cluster_it C1end = clusters_.end();
  for (cluster_it C1 = clusters_.begin(); C1 != C1end; ++C1) {
    cluster_it C2 = C1;
    ++C2;
    for (; C2 != clusters_.end(); ++C2)
      ClusterDistances_.AddElement( Cdist_->CentroidDist( (*C1).Cent(), (*C2).Cent() ) );
  }

  return 0;
}
