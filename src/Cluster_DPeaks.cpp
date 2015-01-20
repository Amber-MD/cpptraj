#include <cmath> // fabs
#include <algorithm> // sort
#include "Cluster_DPeaks.h"
#include "CpptrajStdio.h"
#include "DataSet_Mesh.h"

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
  Carray Points; // Frames to cluster.
  // First determine which frames are being clustered.
  for (int frame = 0; frame < (int)FrameDistances_.Nframes(); ++frame)
    if (!FrameDistances_.IgnoringRow( frame ))
      Points.push_back( Cpoint(frame) );
  // Sanity check.
  if (Points.size() < 2) {
    mprinterr("Error: Only 1 frame in initial clustering.\n");
    return 1;
  }
  // For each point, determine how many others are within epsilon
  for (Carray::iterator point0 = Points.begin();
                        point0 != Points.end(); ++point0)
  {
    int density = 0;
    for (Carray::const_iterator point1 = Points.begin();
                                point1 != Points.end(); ++point1)
    {
      if (point0 != point1) {
        if ( FrameDistances_.GetFdist(point0->Fnum(), point1->Fnum()) < epsilon_ )
          density++;
      }
    }
    point0->SetDensity( density );
  }
  // For each point, find the closest point that has higher density.
  for (Carray::iterator point0 = Points.begin();
                        point0 != Points.end(); ++point0)
  {
    double min_dist = -1.0;
    double max_dist = -1.0;
    for (Carray::const_iterator point1 = Points.begin();
                                point1 != Points.end(); ++point1)
    {
      if (point0 != point1) {
        double dist1_2 = FrameDistances_.GetFdist(point0->Fnum(), point1->Fnum());
        max_dist = std::max(max_dist, dist1_2); 
        if (point1->Density() > point0->Density())
        {
          if (min_dist < 0.0)
            min_dist = dist1_2;
          else
            min_dist = std::min(min_dist, dist1_2);
        }
      }
    }
    // If min_dist is -1 at this point there is no point with higher density
    // i.e. this point has the highest density. Assign it the maximum observed
    // distance.
    if (min_dist < 0.0)
      point0->SetDist( max_dist );
    else
      point0->SetDist( min_dist );
  }
  // Sort by density
  std::sort( Points.begin(), Points.end() );
  // DEBUG - Plot density vs distance for each point.
  CpptrajFile output;
  output.OpenWrite("dpeaks.dat");
  for (Carray::const_iterator point = Points.begin();
                              point != Points.end(); ++point)
    output.Printf("%i %g \"%i\"\n", point->Density(), point->Dist(), point->Fnum()+1);
  output.CloseFile();
  // Choose points for which the min distance to point with higher density is
  // anomalously high.
  // Currently doing this by calculating the running average of density vs 
  // distance, then choosing points with distance > twice the SD of the 
  // running average.
  // NOTE: Store in a mesh data set for now in case we want to spline etc later.
  unsigned int avg_factor = 10;
  unsigned int window_size = Points.size() / avg_factor;
  mprintf("DBG:\tRunning avg window size is %u\n", window_size);
  // FIXME: Handle case where window_size < frames
  DataSet_Mesh runavg;
  unsigned int ra_size = Points.size() - window_size + 1;
  runavg.Allocate1D( ra_size );
  mprintf("DBG:\tRunning avg set should be size %u\n", ra_size);
  CpptrajFile raOut;
  raOut.OpenWrite("runavg.dpeaks.dat");
  double dwindow = (double)window_size;
  double sumx = 0.0;
  double sumy = 0.0;
  for (unsigned int i = 0; i < window_size; i++) {
    sumx += (double)Points[i].Density();
    sumy += Points[i].Dist();
  }
  double avgy = sumy / dwindow;
  runavg.AddXY( sumx / dwindow, avgy );
  raOut.Printf("%g %g\n", sumx / dwindow, avgy );
  for (unsigned int i = 1; i < ra_size; i++) {
    unsigned int nextwin = i + window_size - 1;
    unsigned int prevwin = i - 1;
    sumx = (double)Points[nextwin].Density() - (double)Points[prevwin].Density() + sumx;
    sumy =         Points[nextwin].Dist()    -         Points[prevwin].Dist()    + sumy;
    avgy = sumy / dwindow;
    runavg.AddXY( sumx / dwindow, avgy );
    raOut.Printf("%g %g\n", sumx / dwindow, avgy );
  }
  raOut.CloseFile();
  mprintf("DBG:\tRunning avg set is size %zu\n", runavg.Size());
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
  unsigned int ra_end = Points.size() - 1;
  for (Carray::const_iterator point = Points.begin();
                              point != Points.end(); ++point)
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
    if (delta > ra_sd) raDelta.Printf(" POTENTIAL CLUSTER");
    raDelta.Printf("\n");
  }
  raDelta.CloseFile();
  return 0;
}

void Cluster_DPeaks::ClusterResults(CpptrajFile& outfile) const {
   outfile.Printf("#Algorithm: DPeaks epsilon %g\n", epsilon_);
}

void Cluster_DPeaks::AddSievedFrames() {
  mprintf("FIXME: Adding sieved frames not yet supported.\n");
}
