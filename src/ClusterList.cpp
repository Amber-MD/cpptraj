#include <cmath>
#include <cfloat> // DBL_MAX: TODO: Get rid of
#include <vector>
#include "ClusterList.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"

// XMGRACE colors
const char* ClusterList::XMGRACE_COLOR[] = {
  "white", "black", "red", "green", "blue", "yellow", "brown", "grey", "violet",
  "cyan", "magenta", "orange", "indigo", "maroon", "turquoise", "darkgreen"
};

// CONSTRUCTOR
ClusterList::ClusterList() : debug_(0), Cdist_(0) {}

// DESTRUCTOR
ClusterList::~ClusterList() {
  if (Cdist_ != 0) delete Cdist_;
}

// ClusterList::SetDebug()
/** Set the debug level */
void ClusterList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0) mprintf("ClusterList debug set to %i\n",debug_);
}

// ClusterList::Renumber()
/** Sort clusters by size and renumber starting from 0, where cluster 0
  * is the largest. Also determine best representative frame and calculate 
  * anything dependent on ClusterDistances since sorting destroys indexing 
  * into ClusterDistances.
  */
void ClusterList::Renumber() {
  // Before clusters are renumbered, calculate the average distance of 
  // this cluster to every other cluster.
  // Only do this if ClusterDistances has been set.
  if (ClusterDistances_.Nelements() > 0) {
    double numdist = (double) (clusters_.size() - 1);
    for (cluster_it node = clusters_.begin(); node != clusters_.end(); ++node)
    {
      double avgclusterdist = 0.0;
      for (cluster_it node2 = clusters_.begin();
                      node2 != clusters_.end(); node2++)
      {
        if (node == node2) continue;
        //mprintf("DBG:\t\t%i to %i %f\n",(*node).num, (*node2).num, 
        //        ClusterDistances.GetElement( (*node).num, (*node2).num ));
        avgclusterdist += ClusterDistances_.GetElement( (*node).Num(), 
                                                       (*node2).Num() );
      }
      avgclusterdist /= numdist;
      //mprintf("DBG:\tCluster %i avg dist = %f\n",(*node).num,avgclusterdist);
      (*node).SetAvgDist( avgclusterdist ); 
    }
  }
  // Sort clusters by population 
  clusters_.sort( );
  // Renumber clusters and calculate some cluster properties.
  int newNum = 0;
  for (cluster_it node = clusters_.begin(); node != clusters_.end(); ++node) 
  {
    (*node).SetNum( newNum++ );
    // Find the centroid frame. Since FindCentroidFrame uses FrameDistances 
    // and not ClusterDistances its ok to call after sorting/renumbering.
    if ((*node).FindCentroidFrame( FrameDistances_ )) {
      mprinterr("Error: Could not determine centroid frame for cluster %i\n",
                (*node).Num());
    }
  }
  // TODO: Clear ClusterDistances?
}

// ClusterList::Summary()
/** Print a summary of clusters.  */
void ClusterList::Summary(std::string const& summaryfile, int maxframesIn) {
  CpptrajFile outfile;
  std::vector<double> distances;
  float fmax = (float)maxframesIn;
  if (outfile.OpenWrite(summaryfile)) {
    mprinterr("Error: ClusterList::Summary: Could not set up file.\n");
    return;
  }

  outfile.Printf("%-8s %8s %8s %8s %8s %8s %8s\n","#Cluster","Frames","Frac",
                     "AvgDist","Stdev","Centroid","AvgCDist");
  for (cluster_it node = clusters_.begin();
                  node != clusters_.end(); node++)
  {
    // Calculate size and fraction of total size of this cluster
    int numframes = (*node).Nframes();
    float frac = (float)numframes / fmax;
    double internalAvg = 0.0;
    double internalSD = 0.0;
    int ndistances = ((numframes * numframes) - numframes) / 2;
    if (ndistances > 0) {
      // Calculate average distance between all frames in this cluster
      distances.clear();
      distances.reserve( ndistances );
      ClusterNode::frame_iterator frame2_end = (*node).endframe();
      ClusterNode::frame_iterator frame1_end = frame2_end;
      --frame1_end;
      for (ClusterNode::frame_iterator frm1 = (*node).beginframe();
                                       frm1 != frame1_end; ++frm1)
      {
        // Since this can be called after sieved frames are added back in,
        // need to ensure distances were calcd for these frames.
        if (!FrameDistances_.IgnoringRow(*frm1)) {
          ClusterNode::frame_iterator frm2 = frm1;
          ++frm2;
          for (; frm2 != frame2_end; ++frm2) {
            if (!FrameDistances_.IgnoringRow(*frm2)) {
              distances.push_back( FrameDistances_.GetElement(*frm1, *frm2) );
              internalAvg += distances.back();
            }
          }
        }
      }
      if (!distances.empty()) {
        internalAvg /= (double)distances.size();
        // Calculate standard deviation
        for (std::vector<double>::iterator dval = distances.begin();
                                           dval != distances.end(); ++dval)
        {
          double diff = internalAvg - *dval;
          internalSD += (diff * diff);
        }
        internalSD /= (double)distances.size();
        internalSD = sqrt(internalSD);
      }
    }
    // OUTPUT
    outfile.Printf("%8i %8i %8.3f %8.3f %8.3f %8i %8.3f\n",
                   (*node).Num(), numframes, frac, internalAvg, 
                   internalSD, (*node).CentroidFrame()+1, (*node).AvgDist() );
  } // END loop over clusters
  outfile.CloseFile();
}

// ClusterList::Summary_Half
/** Print a summary of the first half of the data to the second half.
  */
void ClusterList::Summary_Half(std::string const& summaryfile, int maxframesIn,
                               int splitFrame) 
{
  CpptrajFile outfile;
  int half;
  float fmax = (float)maxframesIn;
  if (outfile.OpenWrite(summaryfile)) {
    mprinterr("Error: ClusterList::Summary_Half: Could not set up file.\n");
    return;
  }

  // Calculate halfway point
  if (splitFrame < 0)
    half = maxframesIn / 2;
  else
    half = splitFrame;
  // Header
  outfile.Printf("# 1st < %i <= 2nd\n", half + 1); 
  outfile.Printf("#%-7s %8s %6s %2s %10s %8s %8s %6s %6s\n", 
                 "Cluster", "Total", "Frac", "C#", "Color", 
                 "NumIn1st", "NumIn2nd","Frac1","Frac2");
  int color = 1; // xmgrace color, 1-15
  for (cluster_it node = clusters_.begin();
                  node != clusters_.end(); node++)
  {
    // Calculate size and fraction of total size of this cluster
    int numframes = (*node).Nframes();
    float frac = (float)numframes / fmax;
    int numInFirstHalf = 0;
    int numInSecondHalf = 0;
    // DEBUG
    //mprintf("\tCluster %i\n",(*node).num);
    // Count how many frames are in the first half and how many 
    // are in the second half.
    for (ClusterNode::frame_iterator frame1 = (*node).beginframe();
                                     frame1 != (*node).endframe();
                                     frame1++)
    {
      if (*frame1 < half)
        ++numInFirstHalf;
      else
        ++numInSecondHalf;
    }
    float frac1 = (float) numframes;
    frac1 = ((float) numInFirstHalf) / frac1;
    float frac2 = (float) numframes;
    frac2 = ((float) numInSecondHalf) / frac2;
    outfile.Printf("%-8i %8i %6.2f %2i %10s %8i %8i %6.2f %6.2f\n",
                   (*node).Num(), numframes, frac, color, XMGRACE_COLOR[color],
                   numInFirstHalf, numInSecondHalf, frac1, frac2);
    if (color<15) ++color;
  }
  outfile.CloseFile();
}

// ClusterList::PrintClustersToFile()
/** Print list of clusters in a style similar to ptraj; each cluster is
  * given a line maxframes characters long, with X for each frame that is
  * in the clusters and . for all other frames. Also print out the
  * representative frame numbers.
  */
void ClusterList::PrintClustersToFile(std::string const& filename, int maxframesIn) {
  CpptrajFile outfile;
  std::string buffer;
  
  if ( outfile.OpenWrite(filename) ) {
    mprinterr("Error: PrintClustersToFile: Could not set up file %s\n",
              filename.c_str());
    return;
  }
  outfile.Printf("#Clustering: %u clusters %i frames\n",
                 clusters_.size(), maxframesIn);
  ComputeDBI( outfile );
  for (cluster_it C1_it = clusters_.begin(); 
                  C1_it != clusters_.end(); C1_it++)
  {
    buffer.clear();
    buffer.resize(maxframesIn, '.');
    for (ClusterNode::frame_iterator frame1 = (*C1_it).beginframe();
                                     frame1 != (*C1_it).endframe();
                                     frame1++)
    {
      buffer[ *frame1 ] = 'X';
    }
    buffer += '\n';
    outfile.Write((void*)buffer.c_str(), buffer.size());
  }
  // Print representative frames
  outfile.Printf("#Representative frames:");
  for (cluster_it C = clusters_.begin(); C != clusters_.end(); C++)
    outfile.Printf(" %i",(*C).CentroidFrame()+1);
  outfile.Printf("\n");
  
  outfile.CloseFile();
}

// ClusterList::PrintClusters()
/** Print list of clusters and frame numbers belonging to each cluster.
  */
void ClusterList::PrintClusters() {
  mprintf("CLUSTER: %u clusters, %u frames.\n", clusters_.size(), FrameDistances_.Nframes() );
  for (cluster_it C = clusters_.begin(); C != clusters_.end(); C++) {
    mprintf("\t%8i : ",(*C).Num());
    for (ClusterNode::frame_iterator fnum = (*C).beginframe();
                                     fnum != (*C).endframe(); ++fnum)
      mprintf("%i,",(*fnum)+1);
    mprintf("\n");
  }
}

// -----------------------------------------------------------------------------
// ClusterList::AddCluster()
/** Add a cluster made up of frames specified by the given list of frames to 
  * the cluster list. Cluster # is current cluster list size.
  */
int ClusterList::AddCluster( std::list<int> const& framelistIn ) {
  clusters_.push_back( ClusterNode( Cdist_, framelistIn, clusters_.size() ) );
  return 0;
}

// ClusterList::CalcFrameDistances()
int ClusterList::CalcFrameDistances(std::string const& filename, 
                                    ClusterDist::DsArray const& dataSets,
                                    DistModeType mode, bool useDME, bool nofit, 
                                    bool useMass, std::string const& maskexpr,
                                    int sieve) 
{
  if (dataSets.empty()) {
    mprinterr("Internal Error: CalcFrameDistances: No DataSets given.\n");
    return 1;
  }
  // Base everything off of the first DataSet
  // TODO: Check all DataSet sizes?
  DataSet* dsIn = dataSets[0];
  // Set up internal cluster disance calculation
  if (dsIn->Type() == DataSet::COORDS) {
    if (useDME)
      Cdist_ = new ClusterDist_DME(dsIn, maskexpr);
    else
      Cdist_ = new ClusterDist_RMS(dsIn, maskexpr, nofit, useMass);
  } else {
    if (dataSets.size() == 1)
      Cdist_ = new ClusterDist_Num(dsIn);
    else // TODO: More than just euclid
      Cdist_ = new ClusterDist_Euclid(dataSets);
  }
  // Attempt to load pairwise distances from file if specified
  if (mode == USE_FILE && !filename.empty()) {
    mprintf("\tLoading pair-wise distances from %s\n", filename.c_str());
    if (FrameDistances_.LoadFile( filename, dsIn->Size() )) {
      mprintf("\tLoading pair-wise distances failed - regenerating from frames.\n");
      mode = USE_FRAMES;
    }
  }
  // Calculate pairwise distances from input DataSet. The ignore array will
  // be set up to ignore sieved frames.
  if (mode == USE_FRAMES) {
    mprintf("\tCalculating pair-wise distances.\n");
    FrameDistances_ = Cdist_->PairwiseDist(sieve);
  }
  // Save distances if filename specified and not previously loaded.
  if (mode == USE_FRAMES && !filename.empty()) {
    mprintf("\tSaving pair-wise distances to %s\n", filename.c_str());
    FrameDistances_.SaveFile( filename );
  }
  // DEBUG - Print Frame distances
  if (debug_ > 1) {
    mprintf("INTIAL FRAME DISTANCES:\n");
    FrameDistances_.PrintElements();
  }
  
  return 0;
}  

// ClusterList::AddSievedFrames()
void ClusterList::AddSievedFrames() {
  mprintf("\tRestoring non-sieved frames:");
  // Ensure cluster centroids are up-to-date
  for (cluster_it Cnode = clusters_.begin(); Cnode != clusters_.end(); ++Cnode) 
    (*Cnode).CalculateCentroid( Cdist_ );
  for (int frame = 0; frame < (int)FrameDistances_.Nframes(); ++frame) {
    if (FrameDistances_.IgnoringRow(frame)) {
      //mprintf(" %i [", frame + 1); // DEBUG
      // Which clusters centroid is closest to this frame?
      double mindist = DBL_MAX;
      cluster_it  minNode = clusters_.end();
      for (cluster_it Cnode = clusters_.begin(); Cnode != clusters_.end(); ++Cnode) {
        double dist = Cdist_->FrameCentroidDist(frame, (*Cnode).Cent());
        //mprintf(" %i:%-6.2f", (*Cnode).Num(), dist); // DEBUG
        if (dist < mindist) {
          mindist = dist;
          minNode = Cnode;
        }
      }
      //mprintf(" ], to cluster %i\n", (*minNode).Num()); // DEBUG
      // Add sieved frame to the closest cluster.
      (*minNode).AddFrameToCluster( frame );
    }
  }
  mprintf("\n");
}

// -----------------------------------------------------------------------------
/** The Davies-Bouldin Index (DBI) is a measure of clustering merit; the 
  * smaller the DBI, the better. The DBI is defined as the average, for all 
  * clusters X, of fred, where fred(X) = max, across other clusters Y, of 
  * (Cx + Cy)/dXY ...here Cx is the average distance from points in X to the 
  * centroid, similarly Cy, and dXY is the distance between cluster centroids.
  */
double ClusterList::ComputeDBI(CpptrajFile& outfile) {
  std::vector<double> averageDist;
  averageDist.reserve( clusters_.size() );
  for (cluster_it C1 = clusters_.begin(); C1 != clusters_.end(); ++C1) {
    // Make sure centroid for this cluster is up to date
    (*C1).CalculateCentroid( Cdist_ );
    // Calculate average distance to centroid for this cluster
    averageDist.push_back( (*C1).CalcAvgToCentroid( Cdist_ ) );
    if (outfile.IsOpen())
      outfile.Printf("#Cluster %i has average-distance-to-centroid %f\n", (*C1).Num(),
                     averageDist.back());
  }
  double DBITotal = 0.0;
  unsigned int nc1 = 0;
  for (cluster_it c1 = clusters_.begin(); c1 != clusters_.end(); ++c1, ++nc1) {
    double MaxFred = 0;
    unsigned int nc2 = 0;
    for (cluster_it c2 = clusters_.begin(); c2 != clusters_.end(); ++c2, ++nc2) {
      if (c1 == c2) continue;
      double Fred = averageDist[nc1] + averageDist[nc2];
      Fred /= Cdist_->CentroidDist( (*c1).Cent(), (*c2).Cent() );
      if (Fred > MaxFred)
        MaxFred = Fred;
    }
    DBITotal += MaxFred;
  }
  DBITotal /= (double)clusters_.size();
  if (outfile.IsOpen()) outfile.Printf("#DBI: %f\n", DBITotal);
  return DBITotal;
}

/** The pseudo-F statistic is another measure of clustering goodness. HIGH 
  * values are GOOD. Generally, one selects a cluster-count that gives a peak 
  * in the pseudo-f statistic (or pSF, for short).
  * Formula: A/B, where A = (T - P)/(G-1), and B = P / (n-G). Here n is the 
  * number of points, G is the number of clusters, T is the total distance from
  * the all-data centroid, and P is the sum (for all clusters) of the distances
  * from the cluster centroid.
  */
