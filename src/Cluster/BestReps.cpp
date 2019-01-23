#include "BestReps.h"
#include "../CpptrajStdio.h"

/// Save up to maxSize of the best (lowest) representative scores/frames.
void Cpptraj::Cluster::BestReps::SaveBestRep(RepMap& reps, RepPair const& Dist_Num,
                                             unsigned int maxSize)
{
  if (reps.size() < maxSize)
    reps.insert( Dist_Num );
  else {
    RepMap::reverse_iterator end = reps.rbegin();
    if (Dist_Num.first < end->first) {
      reps.insert( Dist_Num );
      if (reps.size() > maxSize) {
        RepMap::iterator it = reps.end();
        --it;
        reps.erase( it );
      }
    }
  }
}

/// Set given cluster node with best representative frames/scores in reps
void Cpptraj::Cluster::BestReps::SetBestRepFrame(Node& node, RepMap const& reps)
{
  if (!reps.empty()) {
    node.BestReps().clear();
    for (RepMap::const_iterator it = reps.begin(); it != reps.end(); ++it) {
      node.BestReps().push_back( Node::RepPair(it->second, (float)it->first) );
    }
  }
}

/** Find best representative frames for each cluster. */
int Cpptraj::Cluster::BestReps::FindBestRepFrames(RepMethodType type, int nToSave,
                                                  List& clusters, PairwiseMatrix const& pmatrix,
                                                  int debug)
{
  int err = 0;
  switch (type) {
    case CUMULATIVE:
      err = FindBestRepFrames_CumulativeDist(nToSave, clusters, pmatrix); break;
    case CENTROID:
      err = FindBestRepFrames_Centroid(nToSave, clusters, pmatrix); break;
    default:
      mprinterr("Internal Error: BestReps::FindBestRepFrames: Unhandled type.\n");
      err = 1;
  }

  // DEBUG
  if (debug > 0) {
    for (List::cluster_iterator node = clusters.begin(); node != clusters.end(); ++node)
    {
      mprintf("DEBUG: Cluster %i best reps:\n", node->Num());
      for (Node::RepPairArray::const_iterator it = node->BestReps().begin();
                                              it != node->BestReps().end(); ++it)
        mprintf("\t%i (%g)\n", it->second, it->first);
    }
  }

  return err;
}

// ClusterList::FindBestRepFrames_CumulativeDist()
/** Find the frame in each cluster that is the best representative by
  * having the lowest cumulative distance to every other point in the cluster.
  */
int Cpptraj::Cluster::BestReps::FindBestRepFrames_CumulativeDist(int nToSave, List& clusters,
                                                                 PairwiseMatrix const& pmatrix)
{
  int err = 0;
  for (List::cluster_it node = clusters.begin(); node != clusters.end(); ++node) {
    //node->Cent()->Print("centroid." + integerToString(node->Num())); // DEBUG
    //CpptrajFile tmp; // DEBUG
    //tmp.OpenWrite("c"+integerToString(node->Num())+".bestRep.dat"); // DEBUG
    RepMap bestReps;
    for (Node::frame_iterator f1 = node->beginframe(); f1 != node->endframe(); ++f1)
    {
      double cdist = 0.0;
      for (Node::frame_iterator f2 = node->beginframe(); f2 != node->endframe(); ++f2)
      {
        if (f1 != f2)
          cdist += pmatrix.Frame_Distance(*f1, *f2);
      }
      SaveBestRep(bestReps, RepPair(cdist, *f1), nToSave);
      //tmp.Printf("%i %g %g\n", *f1+1, cdist, Cdist_->FrameCentroidDist(*f1, node->Cent()));
    }
    //tmp.CloseFile();
    if (bestReps.empty()) {
      mprinterr("Error: Could not determine represenative frame for cluster %i\n",
                node->Num());
      err++;
    }
    SetBestRepFrame( *node, bestReps );
  }
  return err;
}

/** Find the frame in the cluster that is the best representative by
  * having the lowest distance to the cluster centroid.
  */
int Cpptraj::Cluster::BestReps::FindBestRepFrames_Centroid(int nToSave, List& clusters,
                                                           PairwiseMatrix const& pmatrix)
{
  int err = 0;
  for (List::cluster_it node = clusters.begin(); node != clusters.end(); ++node) {
    //mprintf("DEBUG: FindBestRepFrames_Centroid: Cluster %i\n", node->Num());
    RepMap bestReps;
    //node->Cent()->Print("centroid." + integerToString(node->Num())); // DEBUG
    for (Node::frame_iterator f1 = node->beginframe(); f1 != node->endframe(); ++f1)
    {
      double dist = pmatrix.MetricPtr()->FrameCentroidDist(*f1, node->Cent());
      //mprintf("\t%8i %10.4g %10.4g %i\n", *f1+1, dist, mindist, minframe+1);
      SaveBestRep(bestReps, RepPair(dist, *f1), nToSave);
    }
    if (bestReps.empty()) {
      mprinterr("Error: Could not determine represenative frame for cluster %i\n",
                node->Num());
      err++;
    }
    SetBestRepFrame( *node, bestReps );
  }
  return err;
}
