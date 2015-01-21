#include <cmath> // fabs
#include <algorithm> // sort
#include "Cluster_DPeaks.h"
#include "CpptrajStdio.h"
#include "DataSet_Mesh.h"
#include "ProgressBar.h"

Cluster_DPeaks::Cluster_DPeaks() : epsilon_(-1.0), calc_noise_(false) {}

void Cluster_DPeaks::Help() {
  mprintf("\t[dpeaks epsilon <e> [noise]]\n");
}

int Cluster_DPeaks::SetupCluster(ArgList& analyzeArgs) {
  epsilon_ = analyzeArgs.getKeyDouble("epsilon", -1.0);
  if (epsilon_ <= 0.0) {
    mprinterr("Error: DPeaks requires epsilon to be set and > 0.0\n"
              "Error: Use 'epsilon <e>'\n");
    return 1;
  }
  calc_noise_ = analyzeArgs.hasKey("noise");
  return 0;
}

void Cluster_DPeaks::ClusteringInfo() {
  mprintf("\tDPeaks: Cutoff (epsilon) for determining local density is %g\n", epsilon_);
  if (calc_noise_)
    mprintf("\t\tCalculating noise as all points within epsilon of another cluster.\n");
}

int Cluster_DPeaks::Cluster() {
  mprintf("\tStarting DPeaks clustering.\n");
  Points_.clear();
  // First determine which frames are being clustered.
  for (int frame = 0; frame < (int)FrameDistances_.Nframes(); ++frame)
    if (!FrameDistances_.IgnoringRow( frame ))
      Points_.push_back( Cpoint(frame) );
  // Sanity check.
  if (Points_.size() < 2) {
    mprinterr("Error: Only 1 frame in initial clustering.\n");
    return 1;
  }
  // For each point, determine how many others are within epsilon
  mprintf("\tDetermining local density of each point.\n");
  ProgressBar cluster_progress( Points_.size() );
  for (Carray::iterator point0 = Points_.begin();
                        point0 != Points_.end(); ++point0)
  {
    cluster_progress.Update(point0 - Points_.begin());
    int density = 0;
    for (Carray::const_iterator point1 = Points_.begin();
                                point1 != Points_.end(); ++point1)
    {
      if (point0 != point1) {
        if ( FrameDistances_.GetFdist(point0->Fnum(), point1->Fnum()) < epsilon_ )
          density++;
      }
    }
    point0->SetDensity( density );
  }
  // Sort by density here. Otherwise array indices will be invalid later.
  std::sort( Points_.begin(), Points_.end(), Cpoint::density_sort() );
  // For each point, find the closest point that has higher density.
  mprintf("\tFinding closest neighbor point with higher density for each point.\n");
  cluster_progress.SetupProgress( Points_.size() );
  for (unsigned int idx0 = 0; idx0 != Points_.size(); idx0++)
  {
    cluster_progress.Update( idx0 );
    double min_dist = -1.0;
    double max_dist = -1.0;
    int nearestIdx = -1; // Index of nearest neighbor with higher density
    Cpoint& point0 = Points_[idx0];
    //mprintf("\nDBG:\tSearching for nearest neighbor to idx %u with higher density than %i.\n",
    //        idx0, point0.Density());
    // Since array is sorted by density we can start at the next point.
    for (unsigned int idx1 = idx0+1; idx1 != Points_.size(); idx1++)
    {
        Cpoint const& point1 = Points_[idx1];
        double dist1_2 = FrameDistances_.GetFdist(point0.Fnum(), point1.Fnum());
        max_dist = std::max(max_dist, dist1_2); 
        if (point1.Density() > point0.Density())
        {
          if (min_dist < 0.0) {
            min_dist = dist1_2;
            nearestIdx = (int)idx1;
            //mprintf("DBG:\t\tNeighbor idx %i is first point (density %i), distance %g\n",
            //        nearestIdx, point1.Density(), min_dist);
          } else if (dist1_2 < min_dist) {
            min_dist = dist1_2;
            nearestIdx = (int)idx1;
            //mprintf("DBG:\t\tNeighbor idx %i is closer (density %i, distance %g)\n",
            //        nearestIdx, point1.Density(), min_dist);
          }
        }
    }
    // If min_dist is -1 at this point there is no point with higher density
    // i.e. this point has the highest density. Assign it the maximum observed
    // distance. Check all the points that were skipped.
    if (min_dist < 0.0) {
      for (unsigned int idx1 = 0; idx1 != idx0; idx1++)
        max_dist = std::max(max_dist, FrameDistances_.GetFdist(point0.Fnum(),Points_[idx1].Fnum()));
      point0.SetDist( max_dist );
      //mprintf("DBG:\tPoint %u has no neighbors with higher density."
      //        " Max distance to another point is %g\n", idx0, max_dist);
    } else {
      point0.SetDist( min_dist );
      //mprintf("DBG:\tClosest point to %u with higher density is %i (distance %g)\n",
      //        idx0, nearestIdx, min_dist);
    }
    point0.SetNearestIdx( nearestIdx );
  }
  // DEBUG - Plot density vs distance for each point.
  CpptrajFile output;
  output.OpenWrite("dpeaks.dat");
  output.Printf("%-10s %10s %s %10s %10s\n", "#Density", "Distance", "Frame", "Idx", "Neighbor");
  for (Carray::const_iterator point = Points_.begin();
                              point != Points_.end(); ++point)
    output.Printf("%-10i %10g \"%i\" %10u %10i\n", point->Density(), point->Dist(),
                  point->Fnum()+1, point-Points_.begin(), point->NearestIdx());
  output.CloseFile();
  // Choose points for which the min distance to point with higher density is
  // anomalously high.
  // Currently doing this by calculating the running average of density vs 
  // distance, then choosing points with distance > twice the SD of the 
  // running average.
  // NOTE: Store in a mesh data set for now in case we want to spline etc later.
  unsigned int avg_factor = 10;
  unsigned int window_size = Points_.size() / avg_factor;
  mprintf("DBG:\tRunning avg window size is %u\n", window_size);
  // FIXME: Handle case where window_size < frames
  DataSet_Mesh runavg;
  unsigned int ra_size = Points_.size() - window_size + 1;
  runavg.Allocate1D( ra_size );
  CpptrajFile raOut;
  raOut.OpenWrite("runavg.dpeaks.dat");
  double dwindow = (double)window_size;
  double sumx = 0.0;
  double sumy = 0.0;
  for (unsigned int i = 0; i < window_size; i++) {
    sumx += (double)Points_[i].Density();
    sumy += Points_[i].Dist();
  }
  double avgy = sumy / dwindow;
  runavg.AddXY( sumx / dwindow, avgy );
  raOut.Printf("%g %g\n", sumx / dwindow, avgy );
  for (unsigned int i = 1; i < ra_size; i++) {
    unsigned int nextwin = i + window_size - 1;
    unsigned int prevwin = i - 1;
    sumx = (double)Points_[nextwin].Density() - (double)Points_[prevwin].Density() + sumx;
    sumy =         Points_[nextwin].Dist()    -         Points_[prevwin].Dist()    + sumy;
    avgy = sumy / dwindow;
    runavg.AddXY( sumx / dwindow, avgy );
    raOut.Printf("%g %g\n", sumx / dwindow, avgy );
  }
  raOut.CloseFile();
  double ra_sd;
  double ra_avg = runavg.Avg( ra_sd );
  // Double stdev
  ra_sd *= 2.0;
  mprintf("DBG:\tAvg of running avg set is %g, sd*2.0 is %g\n", ra_avg, ra_sd);
  // For each point, what is the closest running avgd point?
  CpptrajFile raDelta;
  raDelta.OpenWrite("radelta.dat");
  raDelta.Printf("%-10s %10s %10s\n", "#Frame", "RnAvgPos", "Delta");
  unsigned int ra_position = 0;
  unsigned int ra_end = runavg.Size() - 1;
  int cnum = 0;
  for (Carray::iterator point = Points_.begin();
                        point != Points_.end(); ++point)
  {
    if (ra_position != ra_end) {
      // Is the next running avgd point closer to this point?
      while (ra_position != ra_end) {
        double dens  = (double)point->Density();
        double diff0 = fabs( dens - runavg.X(ra_position  ) );
        double diff1 = fabs( dens - runavg.X(ra_position+1) );
        if (diff1 < diff0)
          ++ra_position; // Next running avg position is closer for this point.
        else
          break; // This position is closer.
      }
    }
    double delta = point->Dist() - runavg.Y(ra_position);
    raDelta.Printf("%-10i %10u %10g", point->Fnum()+1, ra_position, delta);
    if (delta > ra_sd) {
      raDelta.Printf(" POTENTIAL CLUSTER %i", cnum);
      point->SetCluster(cnum++);
    }
    raDelta.Printf("\n");
  }
  raDelta.CloseFile();
  int nclusters = cnum;
  mprintf("\tIdentified %i clusters from density vs distance peaks.\n", nclusters);
  // Each remaining point is assigned to the same cluster as its nearest
  // neighbor of higher density. Do this recursively until a cluster
  // center is found.
  cnum = -1;
  for (unsigned int idx = 0; idx != Points_.size(); idx++) {
    if (Points_[idx].Cnum() == -1) {// Point is unassigned.
      AssignClusterNum(idx, cnum);
      //mprintf("Finished recursion for index %i\n\n", idx);
    }
  }
  // Sort by cluster number. NOTE: This invalidates NearestIdx
  std::sort( Points_.begin(), Points_.end(), Cpoint::cnum_sort() );
  // Determine where each cluster starts and stops in Points array
  typedef std::vector<unsigned int> Parray;
  Parray C_start_stop;
  C_start_stop.reserve( nclusters * 2 );
  cnum = -1;
  for (Carray::const_iterator point = Points_.begin(); point != Points_.end(); ++point)
  {
    if (point->Cnum() != cnum) {
      if (!C_start_stop.empty()) C_start_stop.push_back(point - Points_.begin()); // end of cluster
      C_start_stop.push_back(point - Points_.begin()); // beginning of cluster
      cnum = point->Cnum();
    }
  }
  C_start_stop.push_back( Points_.size() ); // end of last cluster
  // Noise calculation.
  if (calc_noise_) {
    // For each cluster find a border region, defined as the set of points
    // assigned to that cluster which are within epsilon of any other
    // cluster.
    // NOTE: Could use a set here to prevent duplicate frames.
    typedef std::vector<Parray> Barray;
    Barray borderIndices( nclusters ); // Hold indices of border points for each cluster.
    for (Parray::const_iterator idx0 = C_start_stop.begin();
                                idx0 != C_start_stop.end(); idx0 += 2)
    {
      int c0 = Points_[*idx0].Cnum();
      mprintf("Cluster %i\n", c0);
      // Check each frame in this cluster.
      for (unsigned int i0 = *idx0; i0 != *(idx0+1); ++i0)
      {
        Cpoint const& point = Points_[i0];
        // Look at each other cluster
        for (Parray::const_iterator idx1 = idx0 + 2;
                                    idx1 != C_start_stop.end(); idx1 += 2)
        {
          int c1 = Points_[*idx1].Cnum();
          // Check each frame in other cluster
          for (unsigned int i1 = *idx1; i1 != *(idx1+1); i1++)
          {
            Cpoint const& other_point = Points_[i1];
            if (FrameDistances_.GetFdist(point.Fnum(), other_point.Fnum()) < epsilon_) {
              mprintf("\tBorder frame: %i (to cluster %i frame %i)\n",
                      point.Fnum() + 1, c1, other_point.Fnum() + 1);
              borderIndices[c0].push_back( i0 );
              borderIndices[c1].push_back( i1 );
            }
          }
        }
      }
    }
    mprintf("Border Frames:\n");
    for (Parray::const_iterator idx = C_start_stop.begin();
                                idx != C_start_stop.end(); idx += 2)
    {
      int c0 = Points_[*idx].Cnum();
      mprintf("\tCluster %u: %u frames: %u border frames:", c0, *(idx+1) - *idx,
              borderIndices[c0].size());
      if (borderIndices[c0].empty())
        mprintf(" No border points.\n");
      else {
        int highestDensity = -1;
        // Find highest density in border region.
        for (Parray::const_iterator bidx = borderIndices[c0].begin();
                                    bidx != borderIndices[c0].end(); ++bidx)
        {
          if (highestDensity == -1)
            highestDensity = Points_[*bidx].Density();
          else
            highestDensity = std::max(highestDensity, Points_[*bidx].Density());
          mprintf(" %i", Points_[*bidx].Fnum()+1);
        }
        mprintf(". Highest density in border= %i\n", highestDensity);
        // Mark any point with density <= highest border density as noise.
        for (unsigned int i = *idx; i != *(idx+1); i++)
        {
          Cpoint& point = Points_[i];
          if (point.Density() <= highestDensity) {
            point.SetCluster( -1 );
            mprintf("\t\tMarking frame %i as noise (density %i)\n", point.Fnum()+1, point.Density());
          }
        }
      }
    }
  }
  // Add the clusters.
  for (Parray::const_iterator idx = C_start_stop.begin();
                              idx != C_start_stop.end(); idx += 2)
  {
    ClusterDist::Cframes frames;
    for (unsigned int i = *idx; i != *(idx+1); i++) {
      if (Points_[i].Cnum() != -1)
        frames.push_back( Points_[i].Fnum() );
    }
    if (!frames.empty())
      AddCluster( frames );
  }
  // Calculate the distances between each cluster based on centroids
  CalcClusterDistances();

  return 0;
}

/** This should never be called for the point with highest density
  * which by definition should be a cluster center.
  */
void Cluster_DPeaks::AssignClusterNum(int idx, int& cnum) {
  // Who is the nearest neighbor with higher density. 
  int neighbor_idx = Points_[idx].NearestIdx();
  //mprintf("Index %i nearest neighbor index %i\n", idx, neighbor_idx);
  // SANITY CHECK
  if (neighbor_idx == -1) {
    mprinterr("Internal Error: In Cluster_DPeaks::AssignClusterNum nearest neighbor is -1.\n");
    return;
  }
  if (Points_[neighbor_idx].Cnum() != -1) {
    // Nearest neighbor has a cluster num assigned.
    cnum = Points_[neighbor_idx].Cnum();
    //mprintf("Neighbor index %i is cluster %i\n", neighbor_idx, cnum);
  } else
    // Ask neighbor to find cluster num.
    AssignClusterNum(neighbor_idx, cnum);
  //mprintf("Index %i cnum %i\n", idx, cnum);
  // At this point cnum should be set. One more sanity check.
  if (cnum == -1) {
    mprinterr("Internal Error: In Cluster_DPeaks::AssignClusterNum could not get"
              " cluster num for index %u.\n", idx);
    return;
  }
  Points_[idx].SetCluster( cnum );
}

void Cluster_DPeaks::ClusterResults(CpptrajFile& outfile) const {
   outfile.Printf("#Algorithm: DPeaks epsilon %g\n", epsilon_);
}

void Cluster_DPeaks::AddSievedFrames() {
  mprintf("FIXME: Adding sieved frames not yet supported.\n");
}
