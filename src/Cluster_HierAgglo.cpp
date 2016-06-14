#include <cfloat> // DBL_MAX
#include "Cluster_HierAgglo.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"

Cluster_HierAgglo::Cluster_HierAgglo() :
  nclusters_(-1),
  epsilon_(-1.0),
  linkage_(AVERAGELINK)
{}

void Cluster_HierAgglo::Help() {
  mprintf("\t[hieragglo [epsilon <e>] [clusters <n>] [linkage|averagelinkage|complete]\n"
          "\t  [epsilonplot <file>]]\n");
}

static const char* LinkageString[] = {
  "single-linkage", "average-linkage", "complete-linkage"
};

int Cluster_HierAgglo::SetupCluster(ArgList& analyzeArgs) {
  nclusters_ = analyzeArgs.getKeyInt("clusters", -1);
  epsilon_ = analyzeArgs.getKeyDouble("epsilon", -1.0);
  if (analyzeArgs.hasKey("linkage"))             linkage_ = SINGLELINK;
  else if (analyzeArgs.hasKey("averagelinkage")) linkage_ = AVERAGELINK;
  else if (analyzeArgs.hasKey("complete"))       linkage_ = COMPLETELINK;
  else linkage_ = AVERAGELINK; // DEFAULT linkage
  std::string epsilonPlot = analyzeArgs.GetStringKey("epsilonplot");
  if (!epsilonPlot.empty()) {
    if (eps_v_n_.OpenWrite( epsilonPlot )) return 1;
    eps_v_n_.Printf("%-12s %12s\n", "#Epsilon", "Nclusters");
  }
  // Determine finish criteria. If nothing specified default to 10 clusters.
  if (nclusters_==-1 && epsilon_==-1.0) {
    mprintf("Warning: cluster: Neither target # of clusters nor epsilon given.\n");
    nclusters_ = 10;
    mprintf("Warning: cluster: Defaulting to %i clusters.\n", nclusters_);
  }
  return 0;
}

void Cluster_HierAgglo::ClusteringInfo() const {
  mprintf("\tHierarchical Agglomerative:");
  if (nclusters_ != -1)
    mprintf(" %i clusters,",nclusters_);
  if (epsilon_ != -1.0)
    mprintf(" epsilon %.3f,",epsilon_);
  mprintf(" %s.\n", LinkageString[linkage_]);
  if (eps_v_n_.IsOpen())
    mprintf("\tWriting epsilon vs # clusters to '%s'\n", eps_v_n_.Filename().full());
}

/** Set up the initial distances between clusters. Should be called before 
  * any clustering is performed. 
  */
void Cluster_HierAgglo::InitializeClusterDistances() {
  // Sets up matrix and ignore array
  ClusterDistances_.SetupMatrix( clusters_.size() );
  // Build initial cluster distances
  if (linkage_==AVERAGELINK) {
#   ifdef NEWCODE
    // Set up matrix to hold sums of distances to clusters.
    SumDistToCluster_.resize(0, clusters_.size());
    for (cluster_iterator C1 = clusters_.begin(); C1 != clusters_.end(); ++C1) {
      for (cluster_iterator C2 = C1; C2 != clusters_.end(); ++C2) {
        if (C1 != C2) {
          double sum = 0.0;
          for (ClusterNode::frame_iterator F1 = C1->beginframe();
                                           F1 != C1->endframe(); ++F1)
            for (ClusterNode::frame_iterator F2 = C2->beginframe();
                                             F2 != C2->endframe(); ++F2)
              sum += FrameDistances().GetFdist( *F1, *F2 );
          SumDistToCluster_.setElement( C1->Num(), C2->Num(), sum );
          double total = (double)(C1->Nframes() * C2->Nframes());
          ClusterDistances_.SetElement( C1->Num(), C2->Num(), sum / total );
        }
      }
    }
#   else
    for (cluster_it C1_it = clusters_.begin();
                    C1_it != clusters_.end(); C1_it++)
      calcAvgDist(C1_it);
#   endif
  } else if (linkage_==SINGLELINK) {
    for (cluster_it C1_it = clusters_.begin();
                    C1_it != clusters_.end(); C1_it++)
      calcMinDist(C1_it);
  } else if (linkage_==COMPLETELINK) {
    for (cluster_it C1_it = clusters_.begin();
                    C1_it != clusters_.end(); C1_it++)
      calcMaxDist(C1_it);
  }
  if (debug_ > 1) {
    mprintf("CLUSTER: INITIAL CLUSTER DISTANCES:\n");
    ClusterDistances_.PrintElements();
  }
}

/** Cluster using a hierarchical agglomerative (bottom-up) approach. All frames
  * start in their own cluster. The closest two clusters are merged, and 
  * distances between the newly merged cluster and all remaining clusters are
  * recalculated according to one of the following metrics:
  * - single-linkage  : The minimum distance between frames in clusters are used.
  * - average-linkage : The average distance between frames in clusters are used.
  * - complete-linkage: The maximum distance between frames in clusters are used.
  */
int Cluster_HierAgglo::Cluster() {
  // If epsilon not given make it huge
  if (epsilon_ == -1.0) epsilon_ = DBL_MAX;
  // If target clusters not given make it 1
  if (nclusters_ == -1) nclusters_ = 1;
  mprintf("\tStarting Hierarchical Agglomerative Clustering:\n");
  ProgressBar cluster_progress(-10);
  // Build initial clusters.
  for (int frame = 0; frame < (int)FrameDistances().OriginalNframes(); frame++) {
    if (!FrameDistances().FrameWasSieved( frame ))
      AddCluster( ClusterDist::Cframes(1, frame) );
  }
  mprintf("\t%i initial clusters.\n", Nclusters());
  // Build initial cluster distance matrix.
  InitializeClusterDistances();
  // DEBUG - print initial clusters
  if (debug_ > 1)
    PrintClusters();
  bool clusteringComplete = false;
  int iterations = 0;
  while (!clusteringComplete) {
    // Merge 2 closest clusters. Clustering complete if closest dist > epsilon.
    if (MergeClosest()) break;
    // If the target number of clusters is reached we are done
    if (Nclusters() <= nclusters_) {
      mprintf("\n\tTarget # of clusters (%i) met (%u), clustering complete.\n", nclusters_,
              Nclusters());
      break;
    }
    if (Nclusters() == 1) clusteringComplete = true; // Sanity check
    cluster_progress.Update( iterations++ );
  }
  mprintf("\tCompleted after %i iterations, %u clusters.\n",iterations,
          Nclusters());
  return 0;
}

void Cluster_HierAgglo::ClusterResults(CpptrajFile& outfile) const {
  outfile.Printf("#Algorithm: HierAgglo linkage %s nclusters %i epsilon %g\n",
                 LinkageString[linkage_], nclusters_, epsilon_);
}

#ifdef TIMER
void Cluster_HierAgglo::Timing(double total) const {
  time_findMin_.WriteTiming(2, "Find min distance", total);
  time_mergeFrames_.WriteTiming(2, "Merge cluster frames", total);
  time_calcLinkage_.WriteTiming(2, "Calculate new linkage", total);
}
#endif

/** Find and merge the two closest clusters. */
int Cluster_HierAgglo::MergeClosest() {
  int C1, C2;
  // Find the minimum distance between clusters. C1 will be lower than C2.
# ifdef TIMER
  time_findMin_.Start();
# endif
  double min = ClusterDistances_.FindMin(C1, C2);
# ifdef TIMER
  time_findMin_.Stop();
# endif
  if (eps_v_n_.IsOpen())
    eps_v_n_.Printf("%12g %12i\n", min, Nclusters());
  if (debug_>0)
    mprintf("\tMinimum found between clusters %i and %i (%f)\n",C1,C2,min);
  // If the minimum distance is greater than epsilon we are done
  if (min > epsilon_) {
    mprintf("\n\tMinimum distance (%f) is greater than epsilon (%f), clustering complete.\n",
            min, epsilon_);
    return 1;
  }

  // Find C1, the number of the cluster to be merged into.
  cluster_it C1_it = clusters_.begin();
  for (; C1_it != clusters_.end(); ++C1_it)
  {
    if ( C1_it->Num() == C1 ) break;
  }
  if (C1_it == clusters_.end()) {
    mprinterr("Error: MergeClosest: C1 (%i) not found.\n",C1);
    return 1;
  }
  // Find C2 - start from C1 since C1 < C2
  cluster_it C2_it = C1_it;
  for (; C2_it != clusters_.end(); ++C2_it) {
    if ( C2_it->Num() == C2 ) break;
  }
  if (C2_it == clusters_.end()) {
    mprinterr("Error: MergeClosest: C2 (%i) not found.\n",C2);
    return 1;
  }

  // Merge the closest clusters, C2 -> C1, remove C2
# ifdef TIMER
  time_mergeFrames_.Start();
# endif
  C1_it->MergeFrames( *C2_it );
  clusters_.erase( C2_it );
# ifdef TIMER
  time_mergeFrames_.Stop();
# endif
  // DEBUG
  if (debug_>1) {
    mprintf("\nAFTER MERGE of %i and %i:\n",C1,C2);
    PrintClusters();
  }
  // Remove all distances having to do with C2
  ClusterDistances_.Ignore(C2);
# ifdef TIMER
  time_calcLinkage_.Start();
# endif
# ifdef NEWCODE
  // Recalculate distances between C1 and all other clusters
  if (linkage_ == AVERAGELINK) { // TODO: Const
    // Update sums and average distances from C1 to other clusters, 
    // excluding any that have already been merged. 
    for (cluster_it C = clusters_.begin(); C != clusters_.end(); ++C) {
      if (!ClusterDistances_.IgnoringRow(C->Num()) &&
           C->Num() != C1 )
      {
        SumDistToCluster_.element( C1, C->Num() ) += SumDistToCluster_.element( C2, C->Num() );
        double nDist = (double)(C1_it->Nframes() * C->Nframes());
        ClusterDistances_.SetElement( C1, C->Num(), 
                                      SumDistToCluster_.element(C1, C->Num()) / nDist );
      }
    }
  } else if (linkage_ == SINGLELINK)
    calcMinDist(C1_it);
  else
    calcMaxDist(C1_it);
# else
  switch (linkage_) {
    case AVERAGELINK : calcAvgDist(C1_it); break;
    case SINGLELINK  : calcMinDist(C1_it); break;
    case COMPLETELINK: calcMaxDist(C1_it); break;
  }
# endif
# ifdef TIMER
  time_calcLinkage_.Stop();
# endif
  if (debug_>2) {
    mprintf("NEW CLUSTER DISTANCES:\n");
    ClusterDistances_.PrintElements();
  }

  return 0;
}

/** Calculate the minimum distance between frames in cluster specified by
  * iterator C1 and frames in all other clusters.
  */
void Cluster_HierAgglo::calcMinDist(cluster_it& C1_it)
{
  // All cluster distances to C1 must be recalcd.
  for (cluster_it C2_it = clusters_.begin();
                  C2_it != clusters_.end(); ++C2_it)
  {
    if (C2_it == C1_it) continue;
    //mprintf("\t\tRecalc distance between %i and %i:\n",C1,newc2);
    // Pick the minimum distance between newc2 and C1
    double min = DBL_MAX;
    for (ClusterNode::frame_iterator c1frames = C1_it->beginframe();
                                     c1frames != C1_it->endframe();
                                     ++c1frames)
    {
      for (ClusterNode::frame_iterator c2frames = C2_it->beginframe();
                                       c2frames != C2_it->endframe();
                                       ++c2frames)
      {
        double Dist = FrameDistances().GetFdist(*c1frames, *c2frames);
        //mprintf("\t\t\tFrame %i to frame %i = %f\n",*c1frames,*c2frames,Dist);
        if ( Dist < min ) min = Dist;
      }
    }
    //mprintf("\t\tMin distance between %i and %i: %f\n",C1,newc2,min);
    ClusterDistances_.SetCdist( C1_it->Num(), C2_it->Num(), min );
  }
}

/** Calculate the maximum distance between frames in cluster specified by
  * iterator C1 and frames in all other clusters.
  */
void Cluster_HierAgglo::calcMaxDist(cluster_it& C1_it)
{
  // All cluster distances to C1 must be recalcd.
  for (cluster_it C2_it = clusters_.begin();
                  C2_it != clusters_.end(); ++C2_it)
  {
    if (C2_it == C1_it) continue;
    //mprintf("\t\tRecalc distance between %i and %i:\n",C1,newc2);
    // Pick the maximum distance between newc2 and C1
    double max = -1.0;
    for (ClusterNode::frame_iterator c1frames = C1_it->beginframe();
                                     c1frames != C1_it->endframe();
                                     ++c1frames)
    {
      for (ClusterNode::frame_iterator c2frames = C2_it->beginframe();
                                       c2frames != C2_it->endframe();
                                       ++c2frames)
      {
        double Dist = FrameDistances().GetFdist(*c1frames, *c2frames);
        //mprintf("\t\t\tFrame %i to frame %i = %f\n",*c1frames,*c2frames,Dist);
        if ( Dist > max ) max = Dist;
      }
    }
    //mprintf("\t\tMax distance between %i and %i: %f\n",C1,newc2,max);
    ClusterDistances_.SetCdist( C1_it->Num(), C2_it->Num(), max );
  }
}

/** Calculate the average distance between frames in cluster specified by
  * iterator C1 and frames in all other clusters.
  */
void Cluster_HierAgglo::calcAvgDist(cluster_it& C1_it)
{
  // All cluster distances to C1 must be recalcd.
  for (cluster_it C2_it = clusters_.begin();
                  C2_it != clusters_.end(); ++C2_it)
  {
    if (C2_it == C1_it) continue;
    //mprintf("\t\tRecalc distance between %i and %i:\n",(*C1_it).Num(),(*C2_it).Num());
    // Pick the minimum distance between newc2 and C1
    double sumDist = 0;
    double N = 0;
    for (ClusterNode::frame_iterator c1frames = C1_it->beginframe();
                                     c1frames != C1_it->endframe();
                                     ++c1frames)
    {
      for (ClusterNode::frame_iterator c2frames = C2_it->beginframe();
                                       c2frames != C2_it->endframe();
                                       ++c2frames)
      {
        double Dist = FrameDistances().GetFdist(*c1frames, *c2frames);
        //mprintf("\t\t\tFrame %i to frame %i = %f\n",*c1frames,*c2frames,Dist);
        sumDist += Dist;
        N++;
      }
    }
    double Dist = sumDist / N;
    //mprintf("\t\tAvg distance between %i and %i: %f\n",(*C1_it).Num(),(*C2_it).Num(),Dist);
    ClusterDistances_.SetCdist( C1_it->Num(), C2_it->Num(), Dist );
  }
}
