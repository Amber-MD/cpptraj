#include "Algorithm_Kmeans.h"
#include "../CpptrajStdio.h"
#include "../ProgressBar.h"

Cpptraj::Cluster::Algorithm_Kmeans::Algorithm_Kmeans() :
  Algorithm(KMEANS),
  nclusters_(0),
  kseed_(-1),
  maxIt_(100),
  mode_(SEQUENTIAL),
  clusterToClusterCentroid_(false)
{}

void Cpptraj::Cluster::Algorithm_Kmeans::Help() {
  mprintf("\t[kmeans clusters <n> [randompoint [kseed <seed>]] [maxit <iterations>]\n");
}

// SetupCluster()
int Cpptraj::Cluster::Algorithm_Kmeans::Setup(ArgList& analyzeArgs) {
  nclusters_ = analyzeArgs.getKeyInt("clusters", -1);
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

// ClusteringInfo()
void Cpptraj::Cluster::Algorithm_Kmeans::Info() const {
  mprintf("\tK-MEANS: Looking for %i clusters.\n", nclusters_);
  if (mode_ == SEQUENTIAL)
    mprintf("\t\tSequentially modify each point.\n");
  else
    mprintf("\t\tRandomly pick points for modification.\n");
  if (kseed_ != -1 && mode_ == RANDOM)
    mprintf("\t\tSeed for random number generator: %i\n", kseed_);
  mprintf("\tCluster to cluster distance will be based on");
  if (clusterToClusterCentroid_)
    mprintf(" best representative from cluster.\n");
  else
    mprintf(" cluster centroids.\n");
}

// ClusterResults()
void Cpptraj::Cluster::Algorithm_Kmeans::Results(CpptrajFile& outfile) const {
  outfile.Printf("#Algorithm: Kmeans nclusters %i maxit %i\n", nclusters_, maxIt_);
}

// Cpptraj::Cluster::Algorithm_Kmeans::Cluster()
int Cpptraj::Cluster::Algorithm_Kmeans::DoClustering(List& clusters,
                                                     Cframes const& framesToCluster,
                                                     PairwiseMatrix const& pmatrix)
{
  if (mode_ == RANDOM)
    RN_.rn_set( kseed_ );

  int pointCount = (int)framesToCluster.size();

  // This array will hold the indices of the points to process each iteration.
  // If sequential this is just 0 -> pointCount. If random this will be 
  // reassigned each iteration.
  Iarray PointIndices;
  PointIndices.reserve( pointCount );
  for (int processIdx = 0; processIdx != pointCount; processIdx++)
    PointIndices.push_back( processIdx );

  if (clusters.empty()) {
    // Determine seeds
    FindKmeansSeeds( framesToCluster, pmatrix );
    // Add the seed clusters
    for (Iarray::const_iterator seedIdx = SeedIndices_.begin();
                                seedIdx != SeedIndices_.end(); ++seedIdx)
    {
      int seedFrame = framesToCluster[ *seedIdx ];
      // A centroid is created for new clusters.
      clusters.AddCluster( Node(pmatrix.MetricPtr(), Cframes(1, seedFrame), clusters.Nclusters()) );
      // NOTE: No need to calc best rep frame, only 1 frame.
      if (debug_ > 0)
        mprintf("Put frame %i in cluster %i (seed index=%i).\n", 
                seedFrame, clusters.back().Num(), *seedIdx);
    }
  } else {
    // Ensure centroids are up to date.
    clusters.UpdateCentroids( pmatrix.MetricPtr() );
  }

  // Assign points in 3 passes. If a point looked like it belonged to cluster A
  // at first, but then we added many other points and altered our cluster 
  // shapes, its possible that we will want to reassign it to cluster B.
  for (int iteration = 0; iteration != maxIt_; iteration++)
  {
    if (mode_ == RANDOM)
      ShufflePoints( PointIndices );
    // Add each point to an existing cluster, and recompute centroid
    mprintf("\tRound %i: ", iteration);
    ProgressBar progress( PointIndices.size() );
    int Nchanged = 0;
    int prog = 0;
    for (Iarray::const_iterator pointIdx = PointIndices.begin();
                                pointIdx != PointIndices.end(); ++pointIdx, ++prog)
    {
      if (debug_ < 1) progress.Update( prog );
      int oldClusterIdx = -1;
//      if ( iteration != 0 || mode_ != SEQUENTIAL) // FIXME: Should this really happen for RANDOM
//      {
        int pointFrame = framesToCluster[ *pointIdx ];
        if (debug_ > 0)
          mprintf("DEBUG: Processing frame %i (index %i)\n", pointFrame, *pointIdx);
        bool pointWasYanked = true;
        if (iteration > 0) {
          // Yank this point out of its cluster, recompute the centroid
          for (List::cluster_it C1 = clusters.begin(); C1 != clusters.end(); ++C1)
          {
            if (C1->HasFrame( pointFrame )) 
            {
              // If this point is alone in its cluster its in the right place
              if (C1->Nframes() == 1) {
                pointWasYanked = false;
                continue; // FIXME: should this be a break?
              }
              //oldBestRep = C1->BestRepFrame(); 
              oldClusterIdx = C1->Num();
              C1->RemoveFrameUpdateCentroid( pmatrix.MetricPtr(), pointFrame ); // TEST
//              C1->RemoveFrameFromCluster( pointFrame );
              //newBestRep = C1->FindBestRepFrame();
//              C1->CalculateCentroid( pmatrix.MetricPtr() );
              if (debug_ > 0)
                mprintf("Remove Frame %i from cluster %i\n", pointFrame, C1->Num());
              //if (clusterToClusterCentroid_) {
              //  if (oldBestRep != NewBestRep)
              //    C1->AlignToBestRep( pmatrix.MetricPtr() ); // FIXME: Only relevant for COORDS dist?
              //  C1->CalculateCentroid( pmatrix.MetricPtr() ); // FIXME: Seems unnessecary to align prior
              //} 
            }
          }
        } else {
          // First iteration. If this point is already in a cluster it is a seed.
          for (List::cluster_it C1 = clusters.begin(); C1 != clusters.end(); ++C1)
          {
            if (C1->HasFrame( pointFrame )) {
              pointWasYanked = false;
              if (debug_ > 0)
                mprintf("Frame %i was already used to seed cluster %i\n", 
                        pointFrame, C1->Num());
              continue; // FIXME break?
            }
          }
        }
        if (pointWasYanked) {
          // Find out what cluster this point is now closest to.
          double closestDist = -1.0;
          List::cluster_it closestCluster = clusters.begin();
          for (List::cluster_it C1 = clusters.begin(); C1 != clusters.end(); ++C1)
          {
            double dist = pmatrix.MetricPtr()->FrameCentroidDist(pointFrame, C1->Cent());
            if (closestDist < 0.0 || dist < closestDist)
            {
              closestDist = dist;
              closestCluster = C1;
            }
          }
          //oldBestRep = closestCluster->BestRepFrame();
          closestCluster->AddFrameUpdateCentroid( pmatrix.MetricPtr(), pointFrame ); // TEST
//          closestCluster->AddFrameToCluster( pointFrame );
          //newBestRep = closestCluster->FindBestFrameFrame();
//          closestCluster->CalculateCentroid( pmatrix.MetricPtr() );
          if (closestCluster->Num() != oldClusterIdx)
          {
            Nchanged++;
            if (debug_ > 0)
              mprintf("Remove Frame %i from cluster %i, but add to cluster %i (dist= %f).\n",
                      pointFrame, oldClusterIdx, closestCluster->Num(), closestDist);
          } else {
            if (debug_ > 0)
              mprintf("Frame %i staying in cluster %i\n", pointFrame, closestCluster->Num());
          }
          if (clusterToClusterCentroid_) {
            //if (oldBestRep != NewBestRep) {
            //    C1->AlignToBestRep( pmatrix.MetricPtr() ); // FIXME: Only relevant for COORDS dist?
            //  C1->CalculateCentroid( pmatrix.MetricPtr() ); // FIXME: Seems unnessecary to align prior
            //}
          }
        }
//      }
    } // END loop over points to cluster
    if (Nchanged == 0) {
      mprintf("\tK-means round %i: No change. Skipping the rest of the iterations.\n", iteration);
      break;
    } else
      mprintf("\tK-means round %i: %i points changed cluster assignment.\n", iteration, Nchanged);
  } // END k-means iterations
  // Remove any empty clusters
  // FIXME: Will there ever be empty clusters?
  clusters.RemoveEmptyClusters();
  // NOTE in PTRAJ here align all frames to best rep 
  return 0;
}

// FindKmeansSeeds()
/** Find some seed-points for K-means clustering. Take the first point as an 
  * arbitrary first choice.  Then, at each iteration, add the point whose total
  * distance from our set of seeds is as large as possible.
  */
int Cpptraj::Cluster::Algorithm_Kmeans::FindKmeansSeeds(Cframes const& FramesToCluster,
                                                        PairwiseMatrix const& pmatrix)
{
  // SeedIndices will hold indices into FramesToCluster
  SeedIndices_.resize( nclusters_, 1 ); // 1 used to be consistent with ptraj

  double bestDistance = 0.0;
  int frameCount = (int)FramesToCluster.size();
  for (int frameIdx = 0; frameIdx != frameCount; frameIdx++)
  {
    int seedFrame = FramesToCluster[ frameIdx ];
    for (int candidateIdx = frameIdx; candidateIdx < frameCount; candidateIdx++)
    {
      int candidateFrame = FramesToCluster[ candidateIdx ];
      double dist = pmatrix.Frame_Distance( seedFrame, candidateFrame );
      if (dist > bestDistance)
      {
        bestDistance = dist;
        SeedIndices_[0] = frameIdx;
        SeedIndices_[1] = candidateIdx;
      }
    }
  }

  for (int seedIdx = 2; seedIdx != nclusters_; seedIdx++)
  {
    bestDistance = 0.0;
    int bestIdx = 0;
    for (int candidateIdx = 0; candidateIdx < frameCount; candidateIdx++)
    {
      // Make sure this candidate isnt already a seed
      bool skipCandidate = false;
      for (int checkIdx = 0; checkIdx != seedIdx; checkIdx++)
      {
        if (SeedIndices_[checkIdx] == candidateIdx) {
          skipCandidate = true;
          break;
        }
      }
      if (!skipCandidate) {
        // Get the closest distance from this candidate to a current seed
        int candidateFrame = FramesToCluster[ candidateIdx ];
        double nearestDist = -1.0;
        for (int checkIdx = 0; checkIdx != seedIdx; checkIdx++)
        {
          int seedFrame = FramesToCluster[ SeedIndices_[checkIdx] ];
          double dist = pmatrix.Frame_Distance( candidateFrame, seedFrame );
          if (dist < nearestDist || nearestDist < 0.0)
            nearestDist = dist;
        }
        // Is this the best so far?
        if (nearestDist > bestDistance)
        {
          bestDistance = nearestDist;
          bestIdx = candidateIdx;
        }
      }
    }
    SeedIndices_[seedIdx] = bestIdx;
  }
  if (debug_ > 0)
    for (unsigned int si = 0; si != SeedIndices_.size(); si++)
      mprintf("DEBUG:\t\tSeedIndices[%u]= %i\n", si, SeedIndices_[si]);
  return 0;
}

/** Use modern version of the Fisher-Yates shuffle to randomly reorder the
  * given points.
  */
void Cpptraj::Cluster::Algorithm_Kmeans::ShufflePoints( Iarray& PointIndices ) {
  for (unsigned int i = PointIndices.size() - 1; i != 1; i--)
  { // 0 <= j <= i
    unsigned int j = (unsigned int)(RN_.rn_gen() * (double)i);
    int temp = PointIndices[j];
    PointIndices[j] = PointIndices[i];
    PointIndices[i] = temp;
  }
  if (debug_ > 0) { 
    mprintf("DEBUG: Shuffled points:");
    for (Iarray::const_iterator it = PointIndices.begin();
                                it != PointIndices.end(); ++it)
      mprintf(" %i", *it);
    mprintf("\n");
  }
}
