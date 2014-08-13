#include "Cluster_Kmeans.h"
#include "CpptrajStdio.h"
#include "Random.h"

Cluster_Kmeans::Cluster_Kmeans() :
  nclusters_(0),
  kseed_(-1),
  maxIt_(100),
  mode_(SEQUENTIAL)
{}

void Cluster_Kmeans::Help() {
  mprintf("\t[kmeans nclusters <n> [randompoint [kseed <seed>]] [maxit <iterations>]\n");
}

int Cluster_Kmeans::SetupCluster(ArgList& analyzeArgs) {
  nclusters_ = analyzeArgs.getKeyInt("nclusters", -1);
  if (nclusters_ < 2) {
    mprinterr("Error: Specify number of clusters > 1 for K-means algorithm.\n");
    return 1;
  }
  if (analyzeArgs.hasKey("randompoint"))
    mode_ = RANDOM;
  else
    mode_ = SEQUENTIAL;
  kseed_ = analyzeArgs.getKeyInt("kseed", -1);
  maxIt_ = analyzeArgs.getKeyInt("maxit", 100);
  return 0;
}

void Cluster_Kmeans::ClusteringInfo() {
  mprintf("\tK-MEANS: Looking for %i clusters.\n", nclusters_);
  if (mode_ == SEQUENTIAL)
    mprintf("\t\tSequentially modify each point.\n");
  else
    mprintf("\t\tRandomly pick point for modification.\n");
  if (kseed_ != -1 && mode_ == RANDOM)
    mprintf("\t\tSeed for random number generator: %i\n", kseed_);
}


int Cluster_Kmeans::Cluster() {
  // First determine which frames are being clustered.
  // FIXME: Can this just be the sieved array?
  for (int frame = 0; frame < (int)FrameDistances_.Nframes(); ++frame)
    if (!FrameDistances_.IgnoringRow( frame ))
      FramesToCluster_.push_back( frame );

  // Determine seeds
  FindKmeansSeeds();

  Random_Number RN;
  if (mode_ == RANDOM)
    RN.rn_set( kseed_ );

  int pointCount = (int)FramesToCluster_.size();
  std::vector<bool> FinishedPoints( pointCount, false );

  // Add the seed clusters
  for (Iarray::const_iterator seedIdx = SeedIndices_.begin();
                              seedIdx != SeedIndices_.end(); ++seedIdx)
  {
    int seedFrame = FramesToCluster_[ *seedIdx ];
    AddCluster( ClusterDist::Cframes(1, seedFrame) );
    clusters_.back().CalculateCentroid( Cdist_ );
    FinishedPoints[ *seedIdx ] = true;
    mprintf("Put frame %i in cluster %i.\n", seedFrame, clusters_.back().Num());
  }
  int unprocessedPointCount = pointCount - nclusters_;
  int oldClusterIdx = -1;
  // Assign points in 3 passes. If a point looked like it belonged to cluster A
  // at first, but then we added many other points and altered our cluster 
  // shapes, its possible that we will want to reassign it to cluster B.
  for (int iteration = 0; iteration != maxIt_; iteration++)
  {
    // Add each point to an existing cluster, and recompute centroid
    mprintf("Round %i\n", iteration);
    if (iteration != 0) {
      FinishedPoints.assign( pointCount, false );
      unprocessedPointCount = pointCount;
    }
    int changed = 0;
    for (int processIdx = 0; processIdx != pointCount; processIdx++)
    {
      if (mode_ == SEQUENTIAL)
        pointIdx = processIdx;
      else
        pointIdx = ChooseNextPoint(FinishedPoints, pointCount, unprocessedPointCount - processIdx);
      if (pointIdx != -1 && 
           (iteration != 0 || mode != SEQUENTIAL || !FinishedPoints[pointIdx]))
      {
        int pointFrame = FramesToCluster_[ pointIdx ];
        bool yanked = true;
        if (iteration > 0) {
          // Yank this point out of its cluster, recompute the centroid
          for (cluster_it C1 = clusters_.begin(); C1 != clusters_.end(); ++C1)
          {
            if (C1->HasFrame( pointFrame )) 
            {
              // If this point is alone in its cluster its in the right place
              if (C1->Nframes() ==1) {
                yanked = false;
                //continue;
                break; // break bc if pointFrame is here it wont be anywhere else
              }
              C1->FindCentroidFrame();
              oldBestRep = C1->CentroidFrame(); 
              oldCluster = C1;
              C1->RemoveFrameFromCluster( pointFrame );
              C1->CalculateCentroid( Cdist_ );
              C1->FindCentroidFrame();
              newBestRep = C1->CentroidFrame();
              mprintf("Remove Frame %i from cluster %i\n", pointFrame, C1->Num());
              // if C2Ccentroid
              // In ptraj realign against best rep here.
            }
          }
        }
        if (Yanked) {
          double closestDist = -1.0;
          cluster_it closestCluster = clusters_.begin();
          for (cluster_it C1 = clusters_.begin(); C1 != clusters_.end(); ++C1)
          {
            double dist = Cdist_->FrameCentroidDist(pointFrame, C1->Cent());
            if (closestDist < 0.0 || dist < closestDist)
            {
              closestDist = dist;
              closestCluster = C1;
            }
          }
          closestCluster->AddFrameToCluster( pointFrame );
          closestCluster->CalculateCentroid( Cdist_ );
        }
      }
    } // END loop over points to cluster 
  } // END k-means iterations
  return 0;
}

/** Find some seed-points for K-means clustering. Take the first point as an 
  * arbitrary first choice.  Then, at each iteration, add the point whose total
  * distance from our set of seeds is as large as possible.
  */
int Cluster_Kmeans::FindKmeansSeeds() {
  // SeedIndices will hold indices into FramesToCluster_
  SeedIndices_.resize( nclusters_, 1 ); // 1 used to be consistent with ptraj

  for (int seedIdx = 1; seedIdx < nclusters_; seedIdx++)
  {
    double bestDistance = 0.0;
    int bestIdx = 0;
    for (int candidateIdx = 0; candidateIdx != (int)FramesToCluster_.size(); candidateIdx++)
    {
      // Make sure this candidate isnt already a seed
      bool skipPoint = false;
      for (int checkIdx = 0; checkIdx < seedIdx; checkIdx++)
      {
        if (SeedIndices_[checkIdx] == candidateIdx) {
          skipPoint = true;
          break;
        }
      }
      if (!skipPoint) {
        // Get the closest distance from this candidate to a current seed
        double nearestDistance = -1.0;
        int candidateFrame = FramesToCluster_[candidateIdx];
        for (int checkIdx = 0; checkIdx < seedIdx; checkIdx++)
        {
          int seedFrame = FramesToCluster_[ SeedIndices_[checkIdx] ];
          double dist = FrameDistances_.GetFdist( seedFrame, candidateFrame );
          if ( dist < nearestDistance || nearestDistance < 0.0 )
            nearestDistance = dist;
        }
        // Is this the best so far
        if (nearestDistance > bestDistance)
        {
          bestDistance = nearestDistance;
          bestIdx = candidateIdx;
        }
      }
    }
    SeedIndices_[seedIdx] = bestIdx;
  }

  return 0;
}
