#include <cmath> // sqrt
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
void ClusterList::Renumber(bool addSievedFrames) {
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
        avgclusterdist += ClusterDistances_.GetCdist( (*node).Num(), (*node2).Num() );
      }
      avgclusterdist /= numdist;
      //mprintf("DBG:\tCluster %i avg dist = %f\n",(*node).num,avgclusterdist);
      (*node).SetAvgDist( avgclusterdist ); 
    }
  }
  // Update cluster centroids.
  bool centroid_error = false;
  for (cluster_it node = clusters_.begin(); node != clusters_.end(); ++node) {
    (*node).SortFrameList();
    // Ensure cluster centroid is up-to-date
    (*node).CalculateCentroid( Cdist_ );
    // Find frame that is closest to the centroid.
    if ((*node).FindCentroidFrame( FrameDistances_ )) {
      mprinterr("Error: Could not determine centroid frame for cluster %i\n",
                (*node).Num());
      centroid_error = true;
    }
  }
  // Add back sieved frames based on distance to cluster centroids.
  if (addSievedFrames) {
    if (centroid_error)
      mprinterr("Error: 1 or more centroids not determined. Cannot add sieved frames.\n");
    else {
      mprintf("\tRestoring sieved frames.\n");
      AddSievedFrames();
    }
    // Re-sort cluster frame lists.
    for (cluster_it node = clusters_.begin(); node != clusters_.end(); ++node)
      (*node).SortFrameList();
  }
  // Sort clusters by population 
  clusters_.sort( );
  // Renumber clusters.
  int newNum = 0;
  for (cluster_it node = clusters_.begin(); node != clusters_.end(); ++node) 
    (*node).SetNum( newNum++ );
  // TODO: Clear ClusterDistances?
}

// ClusterList::Summary()
/** Print a summary of clusters.  */
void ClusterList::Summary(std::string const& summaryfile, int maxframesIn) {
  CpptrajFile outfile;
  double fmax = (double)maxframesIn;
  if (outfile.OpenWrite(summaryfile)) {
    mprinterr("Error: ClusterList::Summary: Could not set up file.\n");
    return;
  }

  outfile.Printf("%-8s %8s %8s %8s %8s %8s %8s\n","#Cluster","Frames","Frac",
                     "AvgDist","Stdev","Centroid","AvgCDist");
  for (cluster_it node = clusters_.begin();
                  node != clusters_.end(); node++)
  {
    // Since there may be a lot of frames do not calculate SD from the
    // mean (which requires either storing distances or two double loops), 
    // instead use SD = sqrt( (SUM[x^2] - ((SUM[x])^2)/N)/N )
    double internalAvg = 0.0;
    double internalSD = 0.0;
    unsigned int Nelements = 0;
    if ((*node).Nframes() > 1) {
      // Calculate average distance between all frames in this cluster
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
              double dist = FrameDistances_.GetFdist(*frm1, *frm2);
              internalAvg += dist;
              internalSD += (dist * dist);
              ++Nelements;
            }
          }
        }
      }
      if (Nelements > 0) {
        double norm = 1.0 / ((double)Nelements);
        internalAvg *= norm;
        internalSD *= norm;
        internalSD -= (internalAvg * internalAvg);
        if (internalSD > 0.0)
          internalSD = sqrt( internalSD );
        else
          internalSD = 0.0;
      }
    }
    // OUTPUT
    outfile.Printf("%8i %8i %8.3f %8.3f %8.3f %8i %8.3f\n",
                   (*node).Num(), (*node).Nframes(), (double)(*node).Nframes()/fmax, internalAvg, 
                   internalSD, (*node).CentroidFrame()+1, (*node).AvgDist() );
  } // END loop over clusters
  outfile.CloseFile();
}

// ClusterList::Summary_Part
/** Print a summary of clustering for specified portions of the overall traj. 
  */
void ClusterList::Summary_Part(std::string const& summaryfile, int maxframesIn,
                               std::vector<int> const& splitFrames)
{
  const char* nExt[] = {"st", "nd", "rd", "th"};
  if (splitFrames.empty()) return; // Sanity check.
  CpptrajFile outfile;
  double fmax = (double)maxframesIn;
  if (outfile.OpenWrite(summaryfile)) {
    mprinterr("Error: Could not open file '%s'.\n", summaryfile.c_str());
    return;
  }

  // Determine number of frames and traj offset for each part.
  outfile.Printf("# 1st");
  std::vector<double> partMax;
  partMax.reserve( splitFrames.size() + 1 );
  std::vector<int> trajOffset;
  trajOffset.reserve( splitFrames.size() + 1);
  trajOffset.push_back( 0 );
  int lastMax = 0;
  unsigned int eidx = 1;
  for (unsigned int sf = 0; sf < splitFrames.size(); sf++)
  {
    partMax.push_back( (double)(splitFrames[sf] - lastMax) );
    lastMax = splitFrames[sf];
    trajOffset.push_back( lastMax );
    outfile.Printf(" <= %i < %u%s", trajOffset.back(), sf+2, nExt[eidx]);
    if (eidx < 3) ++eidx;
  }
  partMax.push_back( (double)(maxframesIn - lastMax) );
  outfile.Printf("\n# ");
  // Print # of frames in each section
  eidx=0;
  for (std::vector<double>::const_iterator pm = partMax.begin(); pm != partMax.end(); ++pm) {
    if (pm != partMax.begin()) outfile.Printf("  ");
    outfile.Printf("%u%s= %.0f", pm - partMax.begin() + 1, nExt[eidx], *pm);
    if (eidx < 3) ++eidx;
  }
  outfile.Printf("\n");
  // DEBUG
  //mprintf("DEBUG: # Frames (offset):");
  //std::vector<int>::const_iterator of = trajOffset.begin();
  //for (std::vector<double>::const_iterator it = partMax.begin();
  //                                         it != partMax.end(); ++it, ++of)
  //  mprintf(" %.0f (%i)", *it, *of);
  //mprintf("\n");
  // Set up bins
  std::vector<int> numInPart(  splitFrames.size() + 1, 0 );
  std::vector<int> firstFrame( splitFrames.size() + 1, -1);

  // Header
  outfile.Printf("#%-7s %8s %8s %2s %10s", "Cluster", "Total", "Frac", "C#", "Color");
  eidx = 0;
  for (unsigned int pm = 1; pm <= partMax.size(); ++pm) {
    outfile.Printf(" %5s%u%2s", "NumIn", pm, nExt[eidx]);
    if (eidx < 3) ++eidx;
  }
  for (unsigned int pm = 1; pm <= partMax.size(); ++pm)
    outfile.Printf(" %7s%u", "Frac", pm);
  for (unsigned int pm = 1; pm <= partMax.size(); ++pm)
    outfile.Printf(" %7s%u", "First", pm);
  outfile.Printf("\n");
  // LOOP OVER CLUSTERS
  int color = 1; // xmgrace color, 1-15
  for (cluster_it node = clusters_.begin();
                  node != clusters_.end(); node++)
  {
    // Calculate size and fraction of total size of this cluster
    int numframes = (*node).Nframes();
    double frac = (double)numframes / fmax;
    std::fill( numInPart.begin(), numInPart.end(), 0 );
    std::fill( firstFrame.begin(), firstFrame.end(), -1 );
    // DEBUG
    //mprintf("\tCluster %i\n",(*node).num);
    // Count how many frames are in each part. 
    for (ClusterNode::frame_iterator frame1 = (*node).beginframe();
                                     frame1 != (*node).endframe();
                                     frame1++)
    {
      unsigned int bin = splitFrames.size();
      for (unsigned int sf = 0; sf < splitFrames.size(); ++sf) {
        if ( *frame1 < splitFrames[sf] ) {
          bin = sf;
          break;
        }
      }
      if (numInPart[ bin ] == 0)
        firstFrame[ bin ] = *frame1 - trajOffset[ bin ] + 1;
      ++numInPart[ bin ];
    }
    outfile.Printf("%-8i %8i %8.4f %2i %10s", (*node).Num(), numframes, frac,
                   color, XMGRACE_COLOR[color]);
    for (std::vector<int>::const_iterator np = numInPart.begin();
                                          np != numInPart.end(); ++np)
      outfile.Printf(" %8i", *np);
    for (unsigned int pm = 0; pm < partMax.size(); ++pm)
      outfile.Printf(" %8.4f", ((double)numInPart[pm]) / partMax[pm]);
    for (std::vector<int>::const_iterator ff = firstFrame.begin();
                                          ff != firstFrame.end(); ++ff)
      outfile.Printf(" %8i", *ff);
    outfile.Printf("\n");
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
  // Call internal info routine.
  ClusterResults( outfile );
  // Do not print trajectory stuff if no filename given (i.e. STDOUT output)
  if (!filename.empty()) {
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
  }
  // Print representative frame numbers
  outfile.Printf("#Representative frames:");
  for (cluster_it C = clusters_.begin(); C != clusters_.end(); C++)
    outfile.Printf(" %i",(*C).CentroidFrame()+1);
  outfile.Printf("\n");
  // Print sieve info if present
  if (FrameDistances_.SieveValue() != 1) {
    if (FrameDistances_.SieveValue() < -1) {
      outfile.Printf("#Sieve value: %i\n#Sieved frames:", -FrameDistances_.SieveValue());
      ClusterSieve::SievedFrames sFrames = FrameDistances_.Sieved();
      for (ClusterSieve::SievedFrames::const_iterator sfrm = sFrames.begin();
                                                      sfrm != sFrames.end(); ++sfrm)
        outfile.Printf(" %i", *sfrm + 1);
      outfile.Printf("\n");
    } else
      outfile.Printf("#Sieve value: %i\n", FrameDistances_.SieveValue());
  } 
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
int ClusterList::AddCluster( ClusterDist::Cframes const& framelistIn ) {
  clusters_.push_back( ClusterNode( Cdist_, framelistIn, clusters_.size() ) );
  return 0;
}

// ClusterList::CalcFrameDistances()
int ClusterList::CalcFrameDistances(std::string const& filename, 
                                    ClusterDist::DsArray const& dataSets,
                                    DistModeType mode, bool useDME, bool nofit, 
                                    bool useMass, std::string const& maskexpr,
                                    int sieve, int sieveSeed) 
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
    // Test that the mask expression is valid
    AtomMask testMask( maskexpr );
    Topology const& dsTop = ((DataSet_Coords*)dsIn)->Top();
    if ( dsTop.SetupIntegerMask( testMask ) ) {
      mprinterr("Error: Could not set up mask [%s] for topology %s\n",
                 maskexpr.c_str(), dsTop.c_str());
      return 1;
    }
    if (testMask.None()) {
      mprinterr("Error: No atoms elected for mask [%s]\n", testMask.MaskString());
      return 1;
    }
    if (useDME)
      Cdist_ = new ClusterDist_DME(dsIn, testMask);
    else
      Cdist_ = new ClusterDist_RMS(dsIn, testMask, nofit, useMass);
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
    // Set up ClusterMatrix with sieve.
    if (FrameDistances_.SetupWithSieve( dsIn->Size(), sieve, sieveSeed )) {
      mprinterr("Error: Could not setup matrix for pair-wise distances.\n");
      return 1; 
    }
    Cdist_->PairwiseDist(FrameDistances_, FrameDistances_.Sieved() );
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
