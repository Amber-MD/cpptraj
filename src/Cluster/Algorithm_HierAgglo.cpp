#include <limits> // double max
#include "Algorithm_HierAgglo.h"
#include "../CpptrajStdio.h"
#include "../ProgressBar.h"

Cpptraj::Cluster::Algorithm_HierAgglo::Algorithm_HierAgglo() :
  Algorithm(HIERAGGLO),
  nclusters_(-1),
  epsilon_(-1.0),
  linkage_(AVERAGELINK)//,
//  includeSievedFrames_(false)
{}

void Cpptraj::Cluster::Algorithm_HierAgglo::Help() {
  mprintf("\t[hieragglo [epsilon <e>] [clusters <n>] [linkage|averagelinkage|complete]\n"
          "\t  [epsilonplot <file>] [includesieved_cdist]]\n");
}

static const char* LinkageString[] = {
  "single-linkage", "average-linkage", "complete-linkage"
};

int Cpptraj::Cluster::Algorithm_HierAgglo::Setup(ArgList& analyzeArgs) {
  nclusters_ = analyzeArgs.getKeyInt("clusters", -1);
  epsilon_ = analyzeArgs.getKeyDouble("epsilon", -1.0);
  if (analyzeArgs.hasKey("linkage"))             linkage_ = SINGLELINK;
  else if (analyzeArgs.hasKey("averagelinkage")) linkage_ = AVERAGELINK;
  else if (analyzeArgs.hasKey("complete"))       linkage_ = COMPLETELINK;
  else linkage_ = AVERAGELINK; // DEFAULT linkage
//  includeSievedFrames_ = analyzeArgs.hasKey("includesieved_cdist");
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

void Cpptraj::Cluster::Algorithm_HierAgglo::Info() const {
    mprintf("\tHierarchical Agglomerative:");
  if (nclusters_ != -1)
    mprintf(" %i clusters,",nclusters_);
  if (epsilon_ != -1.0)
    mprintf(" epsilon %.3f,",epsilon_);
  mprintf(" %s.\n", LinkageString[linkage_]);
  if (eps_v_n_.IsOpen())
    mprintf("\tWriting epsilon vs # clusters to '%s'\n", eps_v_n_.Filename().full());
  /*if (includeSievedFrames_)
    mprintf("\tSieved frames will be included in final cluster distance calculation.\n"
            "Warning: 'includesieved_cdist' may be very slow.\n");
  else
    mprintf("\tSieved frames will not be included in final cluster distance calculation.\n");*/
}

void Cpptraj::Cluster::Algorithm_HierAgglo::Results(CpptrajFile& outfile) const {
  outfile.Printf("#Algorithm: HierAgglo linkage %s nclusters %i epsilon %g\n",
                 LinkageString[linkage_], nclusters_, epsilon_);
}

void Cpptraj::Cluster::Algorithm_HierAgglo::Timing(double total) const {
# ifdef TIMER
  time_findMin_.WriteTiming(2, "Find min distance", total);
  time_mergeFrames_.WriteTiming(2, "Merge cluster frames", total);
  time_calcLinkage_.WriteTiming(2, "Calculate new linkage", total);
# endif
}

/** Default: put all frames in their own starting cluster. */
void Cpptraj::Cluster::Algorithm_HierAgglo::buildInitialClusters(List& clusters,
                                                                 Cframes const& framesToCluster,
                                                                 Metric* metric)
{
  int num = 0;
  for (Cframes_it frm = framesToCluster.begin(); frm != framesToCluster.end(); ++frm)
  {
    Cframes oneframe(1, *frm);
    clusters.AddCluster( Node(metric, oneframe, num++) );
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
int Cpptraj::Cluster::Algorithm_HierAgglo::DoClustering(List& clusters,
                                                        Cframes const& framesToCluster,
                                                        PairwiseMatrix const& pmatrix)
{
  // If epsilon not given make it huge
  if (epsilon_ == -1.0) epsilon_ = std::numeric_limits<double>::max();
  // If target clusters not given make it 1
  if (nclusters_ == -1) nclusters_ = 1;
  mprintf("\tStarting Hierarchical Agglomerative Clustering:\n");
  ProgressBar cluster_progress(-10);
  // Build initial clusters.
  if (clusters.empty())
    buildInitialClusters(clusters, framesToCluster, pmatrix.MetricPtr());
  mprintf("\t%i initial clusters.\n", clusters.Nclusters());
  // Build initial cluster distance matrix.
  ClusterDistances_.SetupMatrix( clusters.Nclusters() );
  double dval = -1.0;
  for (List::cluster_it C1_it = clusters.begin(); C1_it != clusters.end(); C1_it++)
  {
    List::cluster_it C2_it = C1_it;
    C2_it++;
    for(; C2_it != clusters.end(); C2_it++)
    {
      switch (linkage_) {
        case SINGLELINK   : dval = minDist(*C1_it, *C2_it, pmatrix); break;
        case COMPLETELINK : dval = maxDist(*C1_it, *C2_it, pmatrix); break;
        case AVERAGELINK  : dval = avgDist(*C1_it, *C2_it, pmatrix); break;
      }
      ClusterDistances_.SetCdist( C1_it->Num(), C2_it->Num(), dval );
    }
  }
  //InitializeClusterDistances();
  // DEBUG - print initial clusters
  if (debug_ > 1)
    clusters.PrintClusters();
  bool clusteringComplete = false;
  int iterations = 0;
  while (!clusteringComplete) {
    // Merge 2 closest clusters. Clustering complete if closest dist > epsilon.
    if (MergeClosest(clusters, pmatrix)) break;
    // If the target number of clusters is reached we are done
    if (clusters.Nclusters() <= nclusters_) {
      mprintf("\n\tTarget # of clusters (%i) met (%u), clustering complete.\n", nclusters_,
              clusters.Nclusters());
      break;
    }
    if (clusters.Nclusters() == 1) clusteringComplete = true; // Sanity check
    cluster_progress.Update( iterations++ );
  }
  mprintf("\tCompleted after %i iterations, %u clusters.\n",iterations,
          clusters.Nclusters());
  return 0;
}

/** Find and merge the two closest clusters.
  * \return 1 if clustering is complete, 0 otherwise.
  */
int Cpptraj::Cluster::Algorithm_HierAgglo::MergeClosest(List& clusters, PairwiseMatrix const& pmatrix)
{
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
    eps_v_n_.Printf("%12g %12i\n", min, clusters.Nclusters());
  if (debug_>0)
    mprintf("\tMinimum found between clusters %i and %i (%f)\n",C1,C2,min);
  // If the minimum distance is greater than epsilon we are done
  if (min > epsilon_) {
    mprintf("\n\tMinimum distance (%f) is greater than epsilon (%f), clustering complete.\n",
            min, epsilon_);
    return 1;
  }

  // Find C1, the number of the cluster to be merged into.
  List::cluster_it C1_it = clusters.begin();
  for (; C1_it != clusters.end(); ++C1_it)
  {
    if ( C1_it->Num() == C1 ) break;
  }
  if (C1_it == clusters.end()) {
    mprinterr("Error: MergeClosest: C1 (%i) not found.\n",C1);
    return 1;
  }
  // Find C2 - start from C1 since C1 < C2
  List::cluster_it C2_it = C1_it;
  for (; C2_it != clusters.end(); ++C2_it) {
    if ( C2_it->Num() == C2 ) break;
  }
  if (C2_it == clusters.end()) {
    mprinterr("Error: MergeClosest: C2 (%i) not found.\n",C2);
    return 1;
  }

  // Merge the closest clusters, C2 -> C1, remove C2
# ifdef TIMER
  time_mergeFrames_.Start();
# endif
  C1_it->MergeFrames( *C2_it );
  clusters.RemoveCluster( C2_it );
# ifdef TIMER
  time_mergeFrames_.Stop();
# endif
  // DEBUG
  if (debug_>1) {
    mprintf("\nAFTER MERGE of %i and %i:\n",C1,C2);
    clusters.PrintClusters();
  }
  // Remove all distances having to do with C2
  ClusterDistances_.Ignore(C2);
# ifdef TIMER
  time_calcLinkage_.Start();
# endif
  switch (linkage_) {
    case AVERAGELINK : calcAvgDist(C1_it, clusters, pmatrix); break;
    case SINGLELINK  : calcMinDist(C1_it, clusters, pmatrix); break;
    case COMPLETELINK: calcMaxDist(C1_it, clusters, pmatrix); break;
  }
# ifdef TIMER
  time_calcLinkage_.Stop();
# endif
  if (debug_>2) {
    mprintf("NEW CLUSTER DISTANCES:\n");
    ClusterDistances_.PrintElements();
  }

  return 0;
}

/** \return The shortest distance between any two points in C1 and C2. */
double Cpptraj::Cluster::Algorithm_HierAgglo::minDist(Node const& C1,
                                                      Node const& C2,
                                                      PairwiseMatrix const& pmatrix)
{
  double min = std::numeric_limits<double>::max();
  for (Node::frame_iterator c1frames = C1.beginframe(); c1frames != C1.endframe(); ++c1frames)
  {
    for (Node::frame_iterator c2frames = C2.beginframe(); c2frames != C2.endframe(); ++c2frames)
    {
      double Dist = pmatrix.GetFdist(*c1frames, *c2frames);
      //mprintf("\t\t\tFrame %i to frame %i = %f\n",*c1frames,*c2frames,Dist);
      if ( Dist < min ) min = Dist;
    }
  }
  return min;
}

/** \return The longest distance between any two points in C1 and C2. */
double Cpptraj::Cluster::Algorithm_HierAgglo::maxDist(Node const& C1,
                                                      Node const& C2,
                                                      PairwiseMatrix const& pmatrix)
{
  double max = -1.0; 
  for (Node::frame_iterator c1frames = C1.beginframe(); c1frames != C1.endframe(); ++c1frames)
  {
    for (Node::frame_iterator c2frames = C2.beginframe(); c2frames != C2.endframe(); ++c2frames)
    {
      double Dist = pmatrix.GetFdist(*c1frames, *c2frames);
      //mprintf("\t\t\tFrame %i to frame %i = %f\n",*c1frames,*c2frames,Dist);
      if ( Dist > max ) max = Dist;
    }
  }
  return max;
}

/** \return The average distance between points in C1 and C2. */
double Cpptraj::Cluster::Algorithm_HierAgglo::avgDist(Node const& C1,
                                                      Node const& C2,
                                                      PairwiseMatrix const& pmatrix)
{
  double sum = 0.0; 
  for (Node::frame_iterator c1frames = C1.beginframe(); c1frames != C1.endframe(); ++c1frames)
  {
    for (Node::frame_iterator c2frames = C2.beginframe(); c2frames != C2.endframe(); ++c2frames)
    {
      double Dist = pmatrix.GetFdist(*c1frames, *c2frames);
      //mprintf("\t\t\tFrame %i to frame %i = %f\n",*c1frames,*c2frames,Dist);
      sum += Dist;
    }
  }
  return sum / (double)(C1.Nframes() * C2.Nframes());
}

/** Calculate the minimum distance between frames in cluster specified by
  * iterator C1 and frames in all other clusters.
  */
void Cpptraj::Cluster::Algorithm_HierAgglo::calcMinDist(List::cluster_it& C1_it, List& clusters,
                                                        PairwiseMatrix const& pmatrix)
{
  // All cluster distances to C1 must be recalcd.
  for (List::cluster_it C2_it = clusters.begin();
                        C2_it != clusters.end(); ++C2_it)
  {
    if (C2_it != C1_it) {
      double min =  minDist(*C1_it, *C2_it, pmatrix);
      //mprintf("\t\tMin distance between %i and %i: %f\n",C1,newc2,min);
      ClusterDistances_.SetCdist( C1_it->Num(), C2_it->Num(), min );
    }
  }
}

/** Calculate the maximum distance between frames in cluster specified by
  * iterator C1 and frames in all other clusters.
  */
void Cpptraj::Cluster::Algorithm_HierAgglo::calcMaxDist(List::cluster_it& C1_it, List& clusters,
                                                        PairwiseMatrix const& pmatrix)
{
  // All cluster distances to C1 must be recalcd.
  for (List::cluster_it C2_it = clusters.begin();
                        C2_it != clusters.end(); ++C2_it)
  {
    if (C2_it != C1_it) {
      double max = maxDist( *C1_it, *C2_it, pmatrix );
      //mprintf("\t\tMax distance between %i and %i: %f\n",C1,newc2,max);
      ClusterDistances_.SetCdist( C1_it->Num(), C2_it->Num(), max );
    }
  }
}

/** Calculate the average distance between frames in cluster specified by
  * iterator C1 and frames in all other clusters.
  */
void Cpptraj::Cluster::Algorithm_HierAgglo::calcAvgDist(List::cluster_it& C1_it, List& clusters,
                                                        PairwiseMatrix const& pmatrix)
{
  // All cluster distances to C1 must be recalcd.
  for (List::cluster_it C2_it = clusters.begin();
                        C2_it != clusters.end(); ++C2_it)
  {
    if (C2_it != C1_it) {
      double Dist = avgDist( *C1_it, *C2_it, pmatrix );
      //mprintf("\t\tAvg distance between %i and %i: %f\n",(*C1_it).Num(),(*C2_it).Num(),Dist);
      ClusterDistances_.SetCdist( C1_it->Num(), C2_it->Num(), Dist );
    }
  }
}

double Cpptraj::Cluster::Algorithm_HierAgglo::ClusterDistance(Node const& C1, Node const& C2,
                                                              PairwiseMatrix const& pmatrix,
                                                              bool includeSieved,
                                                              Cframes const& sievedOut)
const
{
  double dval = -1.0;
  switch (linkage_) {
    case SINGLELINK   : dval = std::numeric_limits<double>::max(); break;
    case COMPLETELINK : dval = -1.0; break;
    case AVERAGELINK  : dval = 0.0; break;
  }
  unsigned int nvals = 0;

  if (includeSieved) {
    // Include sieved frames
    for (Node::frame_iterator f1 = C1.beginframe(); f1 != C1.endframe(); ++f1)
    {
      for (Node::frame_iterator f2 = C2.beginframe(); f2 != C2.endframe(); ++f2)
      {
        double Dist = pmatrix.Frame_Distance(*f1, *f2);
        //mprintf("\t\t\tFrame %i to frame %i = %f\n",*c1frames,*c2frames,Dist);
        switch (linkage_) {
          case SINGLELINK   : if ( Dist < dval ) dval = Dist; break;
          case COMPLETELINK : if ( Dist > dval ) dval = Dist; break;
          case AVERAGELINK  : dval += Dist; nvals++; break;
        }
      }
    }
  } else {
    // No sieved frames included.
    for (Node::frame_iterator f1 = C1.beginframe(); f1 != C1.endframe(); ++f1)
    {
      if (!sievedOut.HasFrame(*f1)) {
        for (Node::frame_iterator f2 = C2.beginframe(); f2 != C2.endframe(); ++f2)
        {
          if (!sievedOut.HasFrame(*f2)) {
            double Dist = pmatrix.GetFdist(*f1, *f2);
            //mprintf("\t\t\tFrame %i to frame %i = %f\n",*c1frames,*c2frames,Dist);
            switch (linkage_) {
              case SINGLELINK   : if ( Dist < dval ) dval = Dist; break;
              case COMPLETELINK : if ( Dist > dval ) dval = Dist; break;
              case AVERAGELINK  : dval += Dist; nvals++; break;
            }
          }
        }
      }
    }
  }
  if (linkage_ == AVERAGELINK)
    dval /= (double)nvals;

  return dval;
}
