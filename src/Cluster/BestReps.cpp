#include "BestReps.h"
#include "List.h"
#include "MetricArray.h"
#include "Node.h"
#include "../CpptrajStdio.h"

/** CONSTRUCTOR */
Cpptraj::Cluster::BestReps::BestReps() :
  debug_(0),
  nToSave_(0),
  type_(NO_REPS)
{}

/** Save up to maxSize of the best (lowest) representative scores/frames. */
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

/** Set given cluster node with best representative frames/scores in reps */
void Cpptraj::Cluster::BestReps::SetBestRepFrame(Node& node, RepMap const& reps)
{
  if (!reps.empty()) {
    node.BestReps().clear();
    for (RepMap::const_iterator it = reps.begin(); it != reps.end(); ++it) {
      node.BestReps().push_back( Node::RepPair(it->second, it->first) );
    }
  }
}

/** Initialize best rep find. */
int Cpptraj::Cluster::BestReps::InitBestReps(RepMethodType typeIn, int nToSaveIn, int debugIn)
{
  // TODO some error checking
  type_ = typeIn;
  nToSave_ = nToSaveIn;
  debug_ = debugIn;
  return 0;
}

/** Print best reps to stdout. */
void Cpptraj::Cluster::BestReps::PrintBestReps(Node const& node) {
  mprintf("DEBUG: Cluster %i best reps:\n", node.Num());
  for (Node::RepPairArray::const_iterator it = node.BestReps().begin();
                                          it != node.BestReps().end(); ++it)
    mprintf("\t%i (%g)\n", it->first, it->second);
}

/** Find best representative frames for each cluster. */
int Cpptraj::Cluster::BestReps::FindBestRepFrames(List& clusters, MetricArray& pmatrix,
                                                  Cframes const& sievedFrames)
const
{
  int err = 0;
  switch (type_) {
    case CUMULATIVE:
      err = FindBestRepFrames_CumulativeDist(clusters, pmatrix);
      break;
    case CENTROID:
      err = FindBestRepFrames_Centroid(clusters, pmatrix);
      break;
    case CUMULATIVE_NOSIEVE:
      err = FindBestRepFrames_NoSieve_CumulativeDist(clusters, pmatrix, sievedFrames);
      break;
    case NO_REPS:
      mprintf("Warning: Skipping best representative frame calc.\n");
      break;
    default:
      mprinterr("Internal Error: BestReps::FindBestRepFrames: Unhandled type.\n");
      err = 1;
  }

  // DEBUG
  if (debug_ > 0) {
    for (List::cluster_iterator node = clusters.begin(); node != clusters.end(); ++node)
      PrintBestReps(*node);
  }

  return err;
}

/** Find best representative frames for given node. */
int Cpptraj::Cluster::BestReps::FindBestRepFrames(Node& node, MetricArray& pmatrix,
                                                  Cframes const& sievedFrames)
const
{
  int err = 0;
  switch (type_) {
    case CUMULATIVE:
      err = FindBestRepFrames_CumulativeDist(node, pmatrix);
      break;
    case CENTROID:
      err = FindBestRepFrames_Centroid(node, pmatrix);
      break;
    case CUMULATIVE_NOSIEVE:
      err = FindBestRepFrames_NoSieve_CumulativeDist(node, pmatrix, sievedFrames);
      break;
    case NO_REPS:
      mprintf("Warning: Skipping best representative frame calc.\n");
      break;
    default:
      mprinterr("Internal Error: BestReps::FindBestRepFrames: Unhandled type.\n");
      err = 1;
  }

  // DEBUG
  if (debug_ > 0)
    PrintBestReps(node);

  return err;
}

/** Find the frame(s) in the node that is best representative by
  * having the lowest cumulative distance to every other point in the cluster.
  */
int Cpptraj::Cluster::BestReps::FindBestRepFrames_CumulativeDist(Node& node,
                                                                 MetricArray& pmatrix)
const
{
  int err = 0;
  //node.Cent()->Print("centroid." + integerToString(node.Num())); // DEBUG
  //CpptrajFile tmp; // DEBUG
  //tmp.OpenWrite("c"+integerToString(node.Num())+".bestRep.dat"); // DEBUG
  // Handle special cases
  if (node.Nframes() == 1) {
    node.BestReps().clear();
    if (debug_ > 0)
      mprintf("DEBUG: Only 1 frame, best rep: %i\n", node.ClusterFrame(0));
    // Only one frame. That is the best rep.
    node.BestReps().push_back( RepPair(node.ClusterFrame(0), 0.0) );
  } else if (node.Nframes() == 2) {
    node.BestReps().clear();
    if (debug_ > 0)
      mprintf("DEBUG: 2 frames %i and %i, using former as best rep.\n",
              node.ClusterFrame(0), node.ClusterFrame(1));
    // Two frames, distance from f1 to f2 same as f2 to f1. Take f1 by convention.
    node.BestReps().push_back(RepPair(node.ClusterFrame(0),
                              pmatrix.Frame_Distance(node.ClusterFrame(0), node.ClusterFrame(1))));
  } else {
    // Find cumulative distance of each frame to all other frames.
    RepMap bestReps;
    for (Node::frame_iterator f1 = node.beginframe(); f1 != node.endframe(); ++f1)
    {
      double cdist = 0.0;
      for (Node::frame_iterator f2 = node.beginframe(); f2 != node.endframe(); ++f2)
      {
        if (f1 != f2)
          cdist += pmatrix.Frame_Distance(*f1, *f2);
      }
      SaveBestRep(bestReps, RepPair(cdist, *f1), nToSave_);
      //tmp.Printf("%i %g %g\n", *f1+1, cdist, Cdist_->FrameCentroidDist(*f1, node.Cent()));
    }
    //tmp.CloseFile();
    if (bestReps.empty()) {
      mprinterr("Error: Could not determine represenative frame for cluster %i\n",
                node.Num());
      err++;
    }
    SetBestRepFrame( node, bestReps );
    // DEBUG
    if (debug_ > 0) {
      mprintf("DEBUG: Best reps:\n");
      for (RepMap::const_iterator it = bestReps.begin(); it != bestReps.end(); ++it)
        mprintf("\t%i (%20.10E)\n", it->second, it->first);
    }
  }
  return err;
}

// ClusterList::FindBestRepFrames_CumulativeDist()
/** Find the frame in each cluster that is the best representative by
  * having the lowest cumulative distance to every other point in the cluster.
  */
int Cpptraj::Cluster::BestReps::FindBestRepFrames_CumulativeDist(List& clusters,
                                                                 MetricArray& pmatrix)
const
{
  int err = 0; 
  for (List::cluster_it node = clusters.begin(); node != clusters.end(); ++node) {
    err += FindBestRepFrames_CumulativeDist(*node, pmatrix);
  }
  return err;
/*
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
*/
}

/** Find the frame in Node that is the best representative by
  * having the lowest cumulative distance to every other point in the cluster,
  * ignoring sieved frames.
  */
int Cpptraj::Cluster::BestReps::
    FindBestRepFrames_NoSieve_CumulativeDist(Node& node, MetricArray& pmatrix,
                                             Cframes const& sievedFrames)
const
{
  int err = 0;
  RepMap bestReps;
  for (Node::frame_iterator f1 = node.beginframe(); f1 != node.endframe(); ++f1)
  {
    if (!sievedFrames.HasFrame( *f1 )) {
      double cdist = 0.0;
      for (Node::frame_iterator f2 = node.beginframe(); f2 != node.endframe(); ++f2)
      {
        if (f1 != f2 && !sievedFrames.HasFrame( *f2 ))
          //cdist += pmatrix.Cache().CachedDistance(*f1, *f2); // TODO benchmark the two ways
          cdist += pmatrix.Frame_Distance(*f1, *f2);
      }
      SaveBestRep(bestReps, RepPair(cdist, *f1), nToSave_);
    }
  }
  if (bestReps.empty()) {
    mprinterr("Error: Could not determine represenative frame for cluster %i\n",
              node.Num());
    err++;
  }
  SetBestRepFrame( node, bestReps );
  return err;
}

/** Find the frame in each cluster that is the best representative by
  * having the lowest cumulative distance to every other point in the cluster,
  * ignoring sieved frames.
  */
int Cpptraj::Cluster::BestReps::
    FindBestRepFrames_NoSieve_CumulativeDist(List& clusters, MetricArray& pmatrix,
                                             Cframes const& sievedFrames)
const
{
  if (sievedFrames.size() > 0)
    mprintf("Warning: Ignoring sieved frames while looking for best representative.\n");
  int err = 0;
  for (List::cluster_it node = clusters.begin(); node != clusters.end(); ++node) {
    err += FindBestRepFrames_NoSieve_CumulativeDist(*node, pmatrix, sievedFrames);
  }
  return err;
}

/** Find the frame in the node that is the best representative by
  * having the lowest distance to the cluster centroid.
  */
int Cpptraj::Cluster::BestReps::FindBestRepFrames_Centroid(Node& node,
                                                           MetricArray& pmatrix)
const
{
  int err = 0;
  //mprintf("DEBUG: FindBestRepFrames_Centroid: Cluster %i\n", node.Num());
  RepMap bestReps;
  //node.Cent()->Print("centroid." + integerToString(node.Num())); // DEBUG
  for (Node::frame_iterator f1 = node.beginframe(); f1 != node.endframe(); ++f1)
  {
    double dist = pmatrix.FrameCentroidDist(*f1, node.Cent());
    //mprintf("\t%8i %10.4g %10.4g %i\n", *f1+1, dist, mindist, minframe+1);
    SaveBestRep(bestReps, RepPair(dist, *f1), nToSave_);
  }
  if (bestReps.empty()) {
    mprinterr("Error: Could not determine represenative frame for cluster %i\n",
              node.Num());
    err++;
  }
  SetBestRepFrame( node, bestReps );
  return err;
}

/** Find the frame in each cluster that is the best representative by
  * having the lowest distance to the cluster centroid.
  */
int Cpptraj::Cluster::BestReps::FindBestRepFrames_Centroid(List& clusters,
                                                           MetricArray& pmatrix)
const
{
  int err = 0;
  for (List::cluster_it node = clusters.begin(); node != clusters.end(); ++node) {
    err += FindBestRepFrames_Centroid(*node, pmatrix);
  }
  return err;
}
