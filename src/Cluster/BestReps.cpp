#include "BestReps.h"

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

