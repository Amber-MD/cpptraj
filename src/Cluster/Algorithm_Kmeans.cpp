#include "Algorithm_Kmeans.h"
#include "Cframes.h"
#include "List.h"
#include "MetricArray.h"
#include "Node.h"
#include "../ArgList.h"
#include "../CpptrajFile.h"
#include "../CpptrajStdio.h"
#include "../ProgressBar.h"

/** CONSTRUCTOR */
Cpptraj::Cluster::Algorithm_Kmeans::Algorithm_Kmeans() :
  Algorithm(KMEANS),
  nclusters_(0),
  kseed_(-1),
  maxIt_(100),
  mode_(SEQUENTIAL),
  clusterToClusterCentroid_(false)
{}

/** Print help to stdout. */
void Cpptraj::Cluster::Algorithm_Kmeans::Help() {
  mprintf("\t[kmeans clusters <n> [randompoint [kseed <seed>]] [maxit <iterations>]]\n");
}

// SetupCluster()
/** Set up kmeans clustering. */
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
/** Print kmeans clustering info. */
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
/** Write kemans info to info file. */
void Cpptraj::Cluster::Algorithm_Kmeans::Results(CpptrajFile& outfile) const {
  outfile.Printf("#Algorithm: Kmeans nclusters %i maxit %i\n", nclusters_, maxIt_);
}

// Cpptraj::Cluster::Algorithm_Kmeans::Cluster()
/** Perform kmeans clustering. */
int Cpptraj::Cluster::Algorithm_Kmeans::DoClustering(List& clusters,
                                                     Cframes const& framesToCluster,
                                                     MetricArray& pmatrix)
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

  int startIteration = 0;
  if (clusters.empty()) {
    // Determine seeds TODO have FindKmeansSeeds return Iarray with frame #s?
    FindKmeansSeeds( framesToCluster, pmatrix );
    // Add the seed clusters
    for (Iarray::const_iterator seedIdx = SeedIndices_.begin();
                                seedIdx != SeedIndices_.end(); ++seedIdx)
    {
      int seedFrame = framesToCluster[ *seedIdx ];
      // A centroid is created for new clusters.
      clusters.AddCluster( Node(pmatrix, Cframes(1, seedFrame), clusters.Nclusters()) );
      // NOTE: No need to calc best rep frame, only 1 frame.
      if (debug_ > 0)
        mprintf("Put frame %i in cluster %i (seed index=%i).\n", 
                seedFrame, clusters.back().Num(), *seedIdx);
    }
  } else {
    // Clusters already exist.
    mprintf("\t%i existing clusters.\n", clusters.Nclusters());
    if (clusters.Nclusters() > nclusters_) {
      // We currently have more clusters than target clusters.
      // Keep the top nclusters_ clusters and break up the rest.
      mprintf("\tNumber of input clusters %i larger than target number of clusters %i.\n"
              "\tRemoving low-population clusters.\n", clusters.Nclusters(),
              nclusters_);
      clusters.Sort();
      while (clusters.Nclusters() > nclusters_) {
        List::cluster_it lastCluster = clusters.end();
        --lastCluster;
        clusters.RemoveCluster( lastCluster );
      }
      mprintf("\tNow %i existing clusters.\n", clusters.Nclusters());
    } else if (clusters.Nclusters() < nclusters_) {
      // We have fewer clusters than target clusters.
      // Try to find new seeds to make up the difference.
      mprintf("\tNumber of input clusters %i smaller than target number of clusters %i.\n",
              clusters.Nclusters(), nclusters_);
      mprintf("\tWill attempt to find seeds from existing clusters.\n");
      Iarray Seeds = FindSeedsFromClusters(clusters, pmatrix);
      if (Seeds.empty()) {
        mprinterr("Error: Finding seeds from existing clusters failed.\n");
        return 1;
      }
      // Add the seed clusters
      for (Iarray::const_iterator seed = Seeds.begin(); seed != Seeds.end(); ++seed)
      {
        // A centroid is created for new clusters.
        clusters.AddCluster( Node(pmatrix, Cframes(1, *seed), clusters.Nclusters()) );
        // NOTE: No need to calc best rep frame, only 1 frame.
        if (debug_ > 0)
          mprintf("Put frame %i in cluster %i.\n", *seed, clusters.back().Num());
      }
    }

    // Ensure centroids are up to date.
    clusters.UpdateCentroids( pmatrix );
    // Since clusters already exist, go beyond the initial pass.
    startIteration = 1;
  }

  // Assign points in 3 passes. If a point looked like it belonged to cluster A
  // at first, but then we added many other points and altered our cluster 
  // shapes, its possible that we will want to reassign it to cluster B.
  for (int iteration = startIteration; iteration != maxIt_; iteration++)
  {
    if (mode_ == RANDOM) {
      RN_.ShufflePoints( PointIndices );
      if (debug_ > 0) { 
        mprintf("DEBUG: Shuffled points:");
        for (Iarray::const_iterator it = PointIndices.begin();
                                    it != PointIndices.end(); ++it)
          mprintf(" %i", *it);
        mprintf("\n");
      }
    }

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
//      if ( iteration != 0 || mode_ != SEQUENTIAL)
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
                break;
              }
              //oldBestRep = C1->BestRepFrame(); 
              oldClusterIdx = C1->Num();
              C1->RemoveFrameUpdateCentroid( pmatrix, pointFrame ); // TEST
//              C1->RemoveFrameFromCluster( pointFrame );
              //newBestRep = C1->FindBestRepFrame();
//              C1->CalculateCentroid( pmatrix.MetricPtr() );
              if (debug_ > 0)
                mprintf("Remove Frame %i from cluster %i\n", pointFrame, C1->Num());
              //if (clusterToClusterCentroid_) {
              //  if (oldBestRep != NewBestRep)
              //    C1->AlignToBestRep( pmatrix.MetricPtr() ); // Only relevant for COORDS dist?
              //  C1->CalculateCentroid( pmatrix.MetricPtr() ); // TODO Seems unnessecary to align prior
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
              break;
            }
          }
        }
        if (pointWasYanked) {
          // Find out what cluster this point is now closest to.
          double closestDist = -1.0;
          List::cluster_it closestCluster = clusters.begin();
          for (List::cluster_it C1 = clusters.begin(); C1 != clusters.end(); ++C1)
          {
            double dist = pmatrix.FrameCentroidDist(pointFrame, C1->Cent());
            if (closestDist < 0.0 || dist < closestDist)
            {
              closestDist = dist;
              closestCluster = C1;
            }
          }
          //oldBestRep = closestCluster->BestRepFrame();
          closestCluster->AddFrameUpdateCentroid( pmatrix, pointFrame ); // TEST
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
              mprintf("Frame %i staying in cluster %i (dist= %f)\n",
                      pointFrame, closestCluster->Num(), closestDist);
          }
          //if (clusterToClusterCentroid_) {
            //if (oldBestRep != NewBestRep) {
            //    C1->AlignToBestRep( pmatrix.MetricPtr() ); // Only relevant for COORDS dist?
            //  C1->CalculateCentroid( pmatrix.MetricPtr() ); // TODO Seems unnessecary to align prior
            //}
          //}
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
  clusters.RemoveEmptyClusters();
  // NOTE in PTRAJ here align all frames to best rep 
  return 0;
}

// FindSeedsFromClusters
/** Given current clusters, find points that are far away from centroids
  * that can be used as new cluster seeds.
  * \return Array containing frame numbers to be used as seeds.
  */
Cpptraj::Cluster::Algorithm_Kmeans::Iarray
  Cpptraj::Cluster::Algorithm_Kmeans::FindSeedsFromClusters(List& clusters,
                                                            MetricArray& pmatrix)
const
{
  // TODO change algorithm to do normal seed choice, then pick seeds that are
  // far away from existing centroids.
  int nSeedsToFind = nclusters_ - clusters.Nclusters();
  mprintf("\tTarget # seeds= %i\n", nSeedsToFind);
  if (nSeedsToFind < 1) return Iarray();
  // Ensure centroids are up to date
  clusters.UpdateCentroids( pmatrix );
  Iarray Seeds( nSeedsToFind, -1 );
  int nSeedsFound = 0;
  while (nSeedsFound < nSeedsToFind)
  {
    // Find the farthest point in each cluster.
    std::vector<double> MaxDst;
    std::vector<int> MaxFrame;
    for (List::cluster_it node = clusters.begin(); node != clusters.end(); ++node)
    {
      double maxdist = 0;
      int maxframe = -1;

      for (Node::frame_iterator frm = node->beginframe(); frm != node->endframe(); ++frm)
      {
        double dist = pmatrix.FrameCentroidDist(*frm, node->Cent());
        if (dist > maxdist) {
          maxdist = dist;
          maxframe = *frm;
        }
      }
      if (debug_ > 0)
        mprintf("DEBUG: cluster %i seed %i maxdist= %f maxframe= %i\n",
                node->Num(), nSeedsFound, maxdist, maxframe);
      MaxDst.push_back( maxdist );
      MaxFrame.push_back( maxframe );
    }
    if (nSeedsFound == 0) {
      // For the first seed, just choose the one with the largest distance to centroid.
      int maxi = -1;
      double maxd = 0;
      List::cluster_it maxclust = clusters.end();
      int idx = 0;
      for (List::cluster_it node = clusters.begin(); node != clusters.end(); ++node, idx++)
      {
        if (MaxDst[idx] > maxd) {
          maxd = MaxDst[idx];
          maxi = idx;
          maxclust = node;
        }
      }
      maxclust->RemoveFrameUpdateCentroid( pmatrix, MaxFrame[maxi] );
      Seeds[nSeedsFound++] = MaxFrame[maxi];
      if (debug_ > 0)
        mprintf("DEBUG: Frame %i from cluster %i (%f) chosen as first seed.\n",
                MaxFrame[maxi], maxi, MaxDst[maxi]);
    } else {
      // Out of the list of farthest points, choose the one with the largest
      // cumulative distance to existing seeds.
      std::vector<double> cumulativeDist;
      int maxi = -1;
      double maxd = 0;
      List::cluster_it maxclust = clusters.end();
      int idx = 0;
      for (List::cluster_it node = clusters.begin(); node != clusters.end(); ++node, idx++)
      {
        double cdist = 0;
        int frm1 = MaxFrame[idx];
        for (int is = 0; is < nSeedsFound; is++)
          cdist += pmatrix.Frame_Distance(frm1, Seeds[is]);
        if (cdist > maxd) {
          maxd = cdist;
          maxi = idx;
          maxclust = node;
        }
      }
      if (debug_ > 0)
        mprintf("DEBUG: Frame %i from cluster %i (%f, %f) chosen as seed %i.\n",
                MaxFrame[maxi], maxi, MaxDst[maxi], maxd, nSeedsFound);
      maxclust->RemoveFrameUpdateCentroid( pmatrix, MaxFrame[maxi] );
      Seeds[nSeedsFound++] = MaxFrame[maxi];
    }

  }

  return Seeds;
}

// FindKmeansSeeds()
/** Find some seed-points for K-means clustering. Take the first point as an 
  * arbitrary first choice.  Then, at each iteration, add the point whose total
  * distance from our set of seeds is as large as possible.
  */
int Cpptraj::Cluster::Algorithm_Kmeans::FindKmeansSeeds(Cframes const& FramesToCluster,
                                                        MetricArray& pmatrix)
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
