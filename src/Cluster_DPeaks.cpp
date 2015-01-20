#include "Cluster_DPeaks.h"
#include "CpptrajStdio.h"

Cluster_DPeaks::Cluster_DPeaks() : epsilon_(-1.0) {}

void Cluster_DPeaks::Help() {
  mprintf("\t[dpeaks epsilon <e>]\n");
}

int Cluster_DPeaks::SetupCluster(ArgList& analyzeArgs) {
  epsilon_ = analyzeArgs.getKeyDouble("epsilon", -1.0);
  if (epsilon_ <= 0.0) {
    mprinterr("Error: DPeaks requires epsilon to be set and > 0.0\n"
              "Error: Use 'epsilon <e>'\n");
    return 1;
  }
  return 0;
}

void Cluster_DPeaks::ClusteringInfo() {
  mprintf("\tDPeaks: Cutoff (epsilon) for determining local density is %g\n", epsilon_);
}

int Cluster_DPeaks::Cluster() {
  Iarray FramesToCluster;
  Iarray local_density; // For each point, # other points within epsilon
  Darray distances; // For each point, minimum distance to point with higher density
  // First determine which frames are being clustered.
  for (int frame = 0; frame < (int)FrameDistances_.Nframes(); ++frame)
    if (!FrameDistances_.IgnoringRow( frame ))
      FramesToCluster.push_back( frame );
  // Sanity check.
  if (FramesToCluster.size() < 2) {
    mprinterr("Error: Only 1 frame in initial clustering.\n");
    return 1;
  }
  // For each point, determine how many others are within epsilon
  local_density.reserve( FramesToCluster.size() );
  for (Iarray::const_iterator point0 = FramesToCluster.begin();
                              point0 != FramesToCluster.end(); ++point0)
  {
    int density = 0;
    for (Iarray::const_iterator point1 = FramesToCluster.begin();
                                point1 != FramesToCluster.end(); ++point1)
    {
      if (point0 != point1) {
        if ( FrameDistances_.GetFdist(*point0, *point1) < epsilon_ )
          ++density;
      }
    }
    local_density.push_back( density );
  }
  // For each point, find the closest point that has higher density.
  distances.reserve( FramesToCluster.size() );
  Iarray::const_iterator p0_density = local_density.begin();
  for (Iarray::const_iterator point0 = FramesToCluster.begin();
                              point0 != FramesToCluster.end();
                            ++point0, ++p0_density)
  {
    double min_dist = -1.0;
    double max_dist = -1.0;
    Iarray::const_iterator p1_density = local_density.begin();
    for (Iarray::const_iterator point1 = FramesToCluster.begin();
                                point1 != FramesToCluster.end();
                              ++point1, ++p1_density)
    {
      if (point0 != point1) {
        max_dist = std::max(max_dist, FrameDistances_.GetFdist(*point0, *point1)); 
        if (*p1_density > *p0_density)
        {
          if (min_dist < 0.0)
            min_dist = FrameDistances_.GetFdist(*point0, *point1);
          else
            min_dist = std::min(min_dist, FrameDistances_.GetFdist(*point0, *point1));
        }
      }
    }
    // If min_dist is -1 at this point there is no point with higher density
    // i.e. this point has the highest density. Assign it the maximum observed
    // distance.
    if (min_dist < 0.0)
      distances.push_back( max_dist );
    else
      distances.push_back( min_dist );
  }
  // DEBUG - Plot density vs distance for each point.
  CpptrajFile output;
  output.OpenWrite("dpeaks.dat");
  for (unsigned int p = 0; p != FramesToCluster.size(); ++p)
    output.Printf("%i %g \"%i\"\n", local_density[p], distances[p], FramesToCluster[p]+1);
  output.CloseFile();
  // Choose points for which the min distance to point with higher density is
  // anomalously high.
  return 0;
}

void Cluster_DPeaks::ClusterResults(CpptrajFile& outfile) const {
   outfile.Printf("#Algorithm: DPeaks epsilon %g\n", epsilon_);
}

void Cluster_DPeaks::AddSievedFrames() {
  mprintf("FIXME: Adding sieved frames not yet supported.\n");
}
