#include <cfloat> // DBL_MAX
#include <cmath> // sqrt
#include <vector>
#include <algorithm> // sort
#include "ClusterList.h"
#include "CpptrajStdio.h"
#include "Constants.h" // Pseudo-F
#include "ProgressBar.h"
#include "StringRoutines.h"
#ifdef _OPENMP
#  include <omp.h>
#endif
#include "PDBfile.h" // DEBUG

// XMGRACE colors
const char* ClusterList::XMGRACE_COLOR[] = {
  "white", "black", "red", "green", "blue", "yellow", "brown", "grey", "violet",
  "cyan", "magenta", "orange", "indigo", "maroon", "turquoise", "darkgreen"
};

static const char* MetricStringArray[] = {
  "RMSD", "DME", "Symmetry-corrected RMSD", "Data Set(s)"
};

const char* ClusterList::MetricString(DistMetricType dm) {
  return MetricStringArray[dm];
}

// CONSTRUCTOR
ClusterList::ClusterList() : debug_(0), Cdist_(0), frameDistances_(0) {}

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

/** Calculate the distance between the given clusters based on centroids.
  * Centroids MUST be up to date.
  */
double ClusterList::ClusterDistance(ClusterNode const& C1, ClusterNode const& C2) const {
  if (C1.Cent() == 0 || C2.Cent() == 0) {
    mprinterr("Internal Error: One or both centroids are null in ClusterDistance()\n");
    return 0.0;
  }
  return Cdist_->CentroidDist( C1.Cent(), C2.Cent() );
}

// ClusterList::Renumber()
/** Sort clusters by size and renumber starting from 0, where cluster 0
  * is the largest. Also updates cluster centroids and adds back sieved
  * frames if necessary. Ensures cluster node frame lists are sorted.
  */
void ClusterList::Renumber(bool addSievedFrames) {
  // Update cluster centroids in case they need to be used to restore sieved frames.
  for (cluster_it node = clusters_.begin(); node != clusters_.end(); ++node) {
    node->SortFrameList();
    node->CalculateCentroid( Cdist_ );
  }
  // Add back sieved frames
  if (addSievedFrames) {
    mprintf("\tRestoring sieved frames.\n");
    AddSievedFrames();
    // Re-sort cluster frame lists and re-calculate centroids including sieved frames.
    for (cluster_it node = clusters_.begin(); node != clusters_.end(); ++node) {
      node->SortFrameList();
      node->CalculateCentroid( Cdist_ );
    }
  }
  // Sort clusters by population and renumber.
  clusters_.sort( );
  int newNum = 0;
  for (cluster_it node = clusters_.begin(); node != clusters_.end(); ++node)
    node->SetNum( newNum++ );
}

// ClusterList::FindBestRepFrames_CumulativeDist()
/** Find the frame in each cluster that is the best representative by
  * having the lowest cumulative distance to every other point in the cluster.
  */
int ClusterList::FindBestRepFrames_CumulativeDist() {
  int err = 0;
  for (cluster_it node = clusters_.begin(); node != clusters_.end(); ++node) {
    //node->Cent()->Print("centroid." + integerToString(node->Num())); // DEBUG
    //CpptrajFile tmp; // DEBUG
    //tmp.OpenWrite("c"+integerToString(node->Num())+".bestRep.dat"); // DEBUG
    double mindist = DBL_MAX;
    int minframe = -1;
    for (ClusterNode::frame_iterator f1 = node->beginframe();
                                     f1 != node->endframe(); ++f1)
    {
      double cdist = 0.0;
      for (ClusterNode::frame_iterator f2 = node->beginframe(); f2 != node->endframe(); ++f2)
      {
        if (f1 != f2)
          cdist += Frame_Distance(*f1, *f2);
      }
      if (cdist < mindist) {
        mindist = cdist;
        minframe = *f1;
      }
      //tmp.Printf("%i %g %g\n", *f1+1, cdist, Cdist_->FrameCentroidDist(*f1, node->Cent()));
    }
    //tmp.CloseFile();
    if (minframe == -1) {
      mprinterr("Error: Could not determine represenative frame for cluster %i\n",
                node->Num());
      err++;
    }
    node->SetBestRepFrame( minframe );
  }
  return err;
}

/** Find the frame in the cluster that is the best representative by
  * having the lowest distance to the cluster centroid.
  */
int ClusterList::FindBestRepFrames_Centroid() {
  int err = 0;
  for (cluster_it node = clusters_.begin(); node != clusters_.end(); ++node) {
    //mprintf("DEBUG: FindBestRepFrames_Centroid: Cluster %i\n", node->Num());
    double mindist = DBL_MAX;
    int minframe = -1;
    //node->Cent()->Print("centroid." + integerToString(node->Num())); // DEBUG
    for (ClusterNode::frame_iterator f1 = node->beginframe();
                                     f1 != node->endframe(); ++f1)
    {
      double dist = Cdist_->FrameCentroidDist(*f1, node->Cent());
      //mprintf("\t%8i %10.4g %10.4g %i\n", *f1+1, dist, mindist, minframe+1);
      if (dist < mindist) {
        mindist = dist;
        minframe = *f1;
      }
    }
    if (minframe == -1) {
      mprinterr("Error: Could not determine represenative frame for cluster %i\n",
                node->Num());
      err++;
    }
    node->SetBestRepFrame( minframe );
  }
  return err;
}

// ClusterList::DetermineNameWidth()
unsigned int ClusterList::DetermineNameWidth() const {
  // Quick pass through clusters to determine width of cluster names
  unsigned int nWidth = 0;
  for (cluster_iterator node = begincluster(); node != endcluster(); ++node)
    nWidth = std::max(nWidth, (unsigned int)node->Cname().size());
  return nWidth;
}

// ClusterList::Summary()
/** Print a summary of clusters.  */
void ClusterList::Summary(std::string const& summaryfile, bool includeSieveInAvg) const
{
  CpptrajFile outfile;
  double fmax = (double)FrameDistances().OriginalNframes();
  if (outfile.OpenWrite(summaryfile)) {
    mprinterr("Error: ClusterList::Summary: Could not set up file.\n");
    return;
  }
  if (FrameDistances().SieveValue() != 1 && !includeSieveInAvg)
    mprintf("Warning: Within cluster average distance (AvgDist) does not include sieved frames.\n");
  outfile.Printf("%-8s %8s %8s %8s %8s %8s %8s","#Cluster","Frames","Frac",
                     "AvgDist","Stdev","Centroid","AvgCDist");
  unsigned int nWidth = DetermineNameWidth();
  if (nWidth > 0) {
    if (nWidth < 8) nWidth = 8;
    outfile.Printf(" %*s %8s", nWidth, "Name", "RMS");
  }
  outfile.Printf("\n");
  // Calculate distances between clusters.
  Matrix<double> cluster_distances;
  cluster_distances.resize( 0, clusters_.size() );
  for (cluster_iterator c1 = begincluster(); c1 != endcluster(); ++c1)
    for (cluster_iterator c2 = c1; c2 != endcluster(); ++c2)
      if (c2 != c1)
        cluster_distances.addElement( ClusterDistance( *c1, *c2 ) );

  unsigned int idx1 = 0;
  for (cluster_iterator node = begincluster(); node != endcluster(); ++node, ++idx1)
  {
    // Calculate the average distance of this cluster to every other cluster.
    double avgclusterdist = 0.0;
    if (clusters_.size() > 1) {
      unsigned int idx2 = 0;
      for (cluster_iterator node2 = begincluster(); node2 != endcluster(); ++node2, ++idx2)
      {
        if (node != node2)
          avgclusterdist += cluster_distances.element(idx1, idx2);
      }
      avgclusterdist /= (double)(clusters_.size() - 1);
      //mprintf("CLUSTER %i avgclusterdist= %g\n", node->Num(), avgclusterdist);
    }
    // Since there may be a lot of frames do not calculate SD from the
    // mean (which requires either storing distances or two double loops), 
    // instead use SD = sqrt( (SUM[x^2] - ((SUM[x])^2)/N)/N )
    double internalAvg = 0.0;
    double internalSD = 0.0;
    unsigned int Nelements = 0;
    if (node->Nframes() > 1) {
      // Calculate average distance between all frames in this cluster.
      if (includeSieveInAvg) {
        for (ClusterNode::frame_iterator f1 = node->beginframe(); f1 != node->endframe(); ++f1) {
          for (ClusterNode::frame_iterator f2 = f1 + 1; f2 != node->endframe(); ++f2) {
            double dist = Frame_Distance(*f1, *f2);
            internalAvg += dist;
            internalSD += (dist * dist);
            ++Nelements;
          }
        }
      } else {
        for (ClusterNode::frame_iterator f1 = node->beginframe(); f1 != node->endframe(); ++f1) {
          if (!FrameDistances().FrameWasSieved( *f1 )) {
            for (ClusterNode::frame_iterator f2 = f1 + 1; f2 != node->endframe(); ++f2) {
              if (!FrameDistances().FrameWasSieved( *f2 )) {
                double dist = FrameDistances().GetFdist(*f1, *f2);
                internalAvg += dist;
                internalSD += (dist * dist);
                ++Nelements;
              }
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
    outfile.Printf("%8i %8i %8.3f %8.3f %8.3f %8i %8.3f",
                   node->Num(), node->Nframes(), (double)node->Nframes()/fmax, internalAvg, 
                   internalSD, node->BestRepFrame()+1, avgclusterdist );
    if (nWidth > 0)
      outfile.Printf(" %*s %8.3f", nWidth, node->Cname().c_str(), node->RefRms());
    outfile.Printf("\n");
  } // END loop over clusters
  outfile.CloseFile();
}

// ClusterList::Summary_Part
/** Print a summary of clustering for specified portions of the overall traj. 
  */
void ClusterList::Summary_Part(std::string const& summaryfile,
                               std::vector<int> const& splitFrames) const
{
  const char* nExt[] = {"st", "nd", "rd", "th"};
  if (splitFrames.empty()) return; // Sanity check.
  CpptrajFile outfile;
  double fmax = (double)FrameDistances().OriginalNframes();
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
  partMax.push_back( (double)(FrameDistances().OriginalNframes() - lastMax) );
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
  // Determine if cluster names will be output.
  unsigned int nWidth = DetermineNameWidth();
  if (nWidth > 0) {
    if (nWidth < 8) nWidth = 8;
    outfile.Printf(" %*s %6s", nWidth, "Name", "RMS");
  }
  outfile.Printf("\n");
  // LOOP OVER CLUSTERS
  int color = 1; // xmgrace color, 1-15
  for (cluster_iterator node = begincluster(); node != endcluster(); ++node)
  {
    // Calculate size and fraction of total size of this cluster
    int numframes = node->Nframes();
    double frac = (double)numframes / fmax;
    std::fill( numInPart.begin(), numInPart.end(), 0 );
    std::fill( firstFrame.begin(), firstFrame.end(), -1 );
    // DEBUG
    //mprintf("\tCluster %i\n",node->num);
    // Count how many frames are in each part. 
    for (ClusterNode::frame_iterator frame1 = node->beginframe();
                                     frame1 != node->endframe();
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
    outfile.Printf("%-8i %8i %8.4f %2i %10s", node->Num(), numframes, frac,
                   color, XMGRACE_COLOR[color]);
    for (std::vector<int>::const_iterator np = numInPart.begin();
                                          np != numInPart.end(); ++np)
      outfile.Printf(" %8i", *np);
    for (unsigned int pm = 0; pm < partMax.size(); ++pm)
      outfile.Printf(" %8.4f", ((double)numInPart[pm]) / partMax[pm]);
    for (std::vector<int>::const_iterator ff = firstFrame.begin();
                                          ff != firstFrame.end(); ++ff)
      outfile.Printf(" %8i", *ff);
    if (nWidth > 0)
      outfile.Printf(" %*s %6.2f", nWidth, node->Cname().c_str(), node->RefRms());
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
void ClusterList::PrintClustersToFile(std::string const& filename) const {
  CpptrajFile outfile;
  std::string buffer;
  
  if ( outfile.OpenWrite(filename) ) {
    mprinterr("Error: PrintClustersToFile: Could not set up file %s\n",
              filename.c_str());
    return;
  }
  outfile.Printf("#Clustering: %u clusters %i frames\n",
                 clusters_.size(), FrameDistances().OriginalNframes());
  ComputeDBI( outfile );
  ComputePseudoF( outfile );
  // Call internal info routine.
  ClusterResults( outfile );
  // Do not print trajectory stuff if no filename given (i.e. STDOUT output)
  if (!filename.empty()) {
    for (cluster_iterator C1 = begincluster(); C1 != endcluster(); ++C1)
    {
      buffer.clear();
      buffer.resize(FrameDistances().OriginalNframes(), '.');
      for (ClusterNode::frame_iterator f1 = C1->beginframe(); f1 != C1->endframe(); ++f1)
        buffer[ *f1 ] = 'X';
      buffer += '\n';
      outfile.Write((void*)buffer.c_str(), buffer.size());
    }
  }
  // Print representative frame numbers
  outfile.Printf("#Representative frames:");
  for (cluster_iterator C1 = begincluster(); C1 != endcluster(); ++C1)
    outfile.Printf(" %i", C1->BestRepFrame()+1);
  outfile.Printf("\n");
  // Print sieve info if present
  if (FrameDistances().SieveValue() != 1) {
    if (FrameDistances().SieveValue() < -1) {
      outfile.Printf("#Sieve value: %i (random)\n#Sieved frames:", -FrameDistances().SieveValue());
      ClusterSieve::SievedFrames const& sFrames = FrameDistances().FramesToCluster();
      for (ClusterSieve::SievedFrames::const_iterator sfrm = sFrames.begin();
                                                      sfrm != sFrames.end(); ++sfrm)
        outfile.Printf(" %i", *sfrm + 1);
      outfile.Printf("\n");
    } else
      outfile.Printf("#Sieve value: %i\n", FrameDistances().SieveValue());
  }
  outfile.CloseFile();
}

// ClusterList::PrintClusters()
/** Print list of clusters and frame numbers belonging to each cluster.
  */
void ClusterList::PrintClusters() const {
  mprintf("CLUSTER: %u clusters, %u frames.\n", clusters_.size(),
          FrameDistances().OriginalNframes() );
  for (cluster_iterator C = begincluster(); C != endcluster(); C++) {
    mprintf("\t%8i : ",C->Num());
    for (ClusterNode::frame_iterator fnum = C->beginframe();
                                     fnum != C->endframe(); ++fnum)
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

/** Set up the cluster distance calculation for the given input data sets
  * and distance metric type.
  */
int ClusterList::SetupCdist( ClusterDist::DsArray const& dataSets,
                             DistMetricType metric, bool nofit, bool useMass,
                             std::string const& maskexpr )
{
  if (dataSets.empty()) { // SANITY CHECK
    mprinterr("Internal Error: SetupCdist: No DataSets given.\n");
    return 1;
  }
  // Base everything off of the first DataSet
  DataSet* dsIn = dataSets[0];
  // Set up internal cluster disance calculation
  if (metric != DATA) {
    if (dsIn->Group() != DataSet::COORDINATES ) {
      mprinterr("Internal Error: Metric is COORDS base but data set is not.\n");
      return 1;
    }
    // Test that the mask expression is valid
    AtomMask testMask( maskexpr );
    Topology const& dsTop = ((DataSet_Coords*)dsIn)->Top();
    if ( dsTop.SetupIntegerMask( testMask ) ) {
      mprinterr("Error: Could not set up mask '%s' for topology %s\n",
                 maskexpr.c_str(), dsTop.c_str());
      return 1;
    }
    testMask.MaskInfo();
    if (testMask.None()) {
      mprinterr("Error: No atoms elected for mask '%s'\n", testMask.MaskString());
      return 1;
    }
    switch (metric) {
      case DME:   Cdist_ = new ClusterDist_DME(dsIn, testMask); break;
      case RMS:   Cdist_ = new ClusterDist_RMS(dsIn, testMask, nofit, useMass); break;
      case SRMSD: Cdist_ = new ClusterDist_SRMSD(dsIn, testMask, nofit, useMass, debug_); break;
      default: return 1; // Sanity check
    }
  } else { // Metric is DATA
    if (dataSets.size() == 1)
      Cdist_ = new ClusterDist_Num(dsIn);
    else // TODO: More than just euclid
      Cdist_ = new ClusterDist_Euclid(dataSets);
  }
  if (debug_ > 0) mprintf("DEBUG: ClusterDist= %s\n", Cdist_->Description().c_str());
  return 0;
}

// ClusterList::CalcFrameDistances()
/** Set the frame pairwise distance matrix from the given data set and
  * calculate it if need be.
  */
int ClusterList::CalcFrameDistances(DataSet* pwDistMatrixIn, 
                                    ClusterDist::DsArray const& dataSets,
                                    int sieve, int sieveSeed) 
{
  if (dataSets.empty()) { // SANITY CHECK
    mprinterr("Internal Error: CalcFrameDistances: No DataSets given.\n");
    return 1;
  }
  if (Cdist_ == 0) { // SANITY CHECK
    mprinterr("Internal Error: ClusterDist for given metric not yet allocated.\n");
    return 1;
  }
  frameDistances_ = (DataSet_Cmatrix*)pwDistMatrixIn;
  if ( FrameDistances().NeedsSetup() ) {
    // Set up cluster matrix with sieving info. Base total number
    // of frames on first DataSet size.
    if (frameDistances_->SetupWithSieve( Cdist_, dataSets[0]->Size(), sieve, sieveSeed ))
      return 1;
    // If cluster matrix needs calculation (i.e. not NOMEM), perform it.
    if (FrameDistances().NeedsCalc()) {
      mprintf("\tCalculating pair-wise distances.\n");
      ClusterSieve::SievedFrames const& frames = FrameDistances().FramesToCluster();
      int f2end = (int)frames.size();
      int f1end = f2end - 1;
      ParallelProgress progress(f1end);
      int f1, f2;
      // For OMP, every other thread will need its own Cdist.
      ClusterDist* MyCdist = Cdist_;
#     ifdef _OPENMP
#     pragma omp parallel private(MyCdist, f1, f2) firstprivate(progress)
      {
      int mythread = omp_get_thread_num();
      progress.SetThread( mythread );
      if (mythread == 0) {
        mprintf("\tParallelizing pairwise distance calc with %i threads\n", omp_get_num_threads());
        MyCdist = Cdist_;
      } else
        MyCdist = Cdist_->Copy();
#     pragma omp for schedule(dynamic)
#     endif
      for (f1 = 0; f1 < f1end; f1++) {
        progress.Update(f1);
        for (f2 = f1 + 1; f2 < f2end; f2++)
          frameDistances_->SetElement( f1, f2, MyCdist->FrameDist(frames[f1], frames[f2]) );
      }
#     ifdef _OPENMP
      if (mythread > 0)
        delete MyCdist;
      } // END omp parallel
#     endif
      progress.Finish();
    }
    // Currently this is only for DataSet_Cmatrix_DISK
    frameDistances_->Complete();
  } else
    // Pairwise distance matrix already set up
    mprintf("\tUsing existing pairwise distances from '%s'\n", FrameDistances().legend());
  mprintf("\tMemory used by pair-wise matrix and other cluster data: %s\n",
          ByteString(FrameDistances().DataSize(), BYTE_DECIMAL).c_str());
  // DEBUG - Print Frame distances
  if (debug_ > 1) {
    mprintf("INITIAL FRAME DISTANCES:\n");
    FrameDistances().PrintElements();
  }
  return 0;
}

// ClusterList::RemoveEmptyClusters()
void ClusterList::RemoveEmptyClusters() {
  cluster_it cnode = clusters_.begin();
  while (cnode != clusters_.end()) {
    if (cnode->Nframes() == 0)
      cnode = clusters_.erase( cnode );
    else
      ++cnode;
  }
}

// -----------------------------------------------------------------------------
void ClusterList::AddSievedFramesByCentroid() {
    // NOTE: All cluster centroids must be up to date.
  int frame;
  int nframes = (int)FrameDistances().OriginalNframes();
  double mindist, dist;
  cluster_it minNode, Cnode;
  ParallelProgress progress( nframes );
  // For OMP, every other thread will need its own Cdist.
  ClusterDist* MyCdist = Cdist_;
# ifdef _OPENMP
  // For OMP need a temp. array to hold which frame goes to which cluster to avoid clashes
  std::vector<cluster_it> frameToCluster( nframes, clusters_.end() );
# pragma omp parallel private(MyCdist, frame, dist, mindist, minNode, Cnode) firstprivate(progress)
  {
  int mythread = omp_get_thread_num();
  progress.SetThread( mythread );
  if (mythread == 0) {
    mprintf("\tParallelizing sieve restore calc with %i threads\n", omp_get_num_threads());
    MyCdist = Cdist_;
  } else
    MyCdist = Cdist_->Copy();
# pragma omp for schedule(dynamic)
# endif
  for (frame = 0; frame < nframes; ++frame) {
    progress.Update( frame );
    if (FrameDistances().FrameWasSieved(frame)) {
      // Which clusters centroid is closest to this frame?
      mindist = DBL_MAX;
      minNode = clusters_.end();
      for (Cnode = clusters_.begin(); Cnode != clusters_.end(); ++Cnode) {
        dist = MyCdist->FrameCentroidDist(frame, Cnode->Cent());
        if (dist < mindist) {
          mindist = dist;
          minNode = Cnode;
        }
      }
      // Add sieved frame to the closest cluster.
#     ifdef _OPENMP
      frameToCluster[frame] = minNode;
#     else
      minNode->AddFrameToCluster( frame );
#     endif
    }
  } // END loop over frames
# ifdef _OPENMP
  if (mythread > 0)
    delete MyCdist;
  } // END pragma omp parallel
  // Now actually add sieved frames to their appropriate clusters
  for (frame = 0; frame < nframes; frame++)
    if (frameToCluster[frame] != clusters_.end())
      (*frameToCluster[frame]).AddFrameToCluster( frame );
# endif
  progress.Finish();
}

// -----------------------------------------------------------------------------
/** The Davies-Bouldin Index (DBI) measures the average similarity between each
  * cluster and its most similar one; the smaller the DBI, the better. The DBI 
  * is defined as the average, for all clusters X, of fred, where fred(X) = max,
  * across other clusters Y, of (Cx + Cy)/dXY. Here Cx is the average distance
  * from points in X to the centroid, similarly Cy, and dXY is the distance 
  * between cluster centroids.
  * NOTE: To use this, cluster centroids should be fully up-to-date.
  */
double ClusterList::ComputeDBI(CpptrajFile& outfile) const {
  std::vector<double> averageDist;
  averageDist.reserve( clusters_.size() );
  for (cluster_iterator C1 = begincluster(); C1 != endcluster(); ++C1) {
    // Calculate average distance to centroid for this cluster
    averageDist.push_back( C1->CalcAvgToCentroid( Cdist_ ) );
    if (outfile.IsOpen())
      outfile.Printf("#Cluster %i has average-distance-to-centroid %f\n", 
                     C1->Num(), averageDist.back());
  }
  double DBITotal = 0.0;
  unsigned int nc1 = 0;
  for (cluster_iterator c1 = begincluster(); c1 != endcluster(); ++c1, ++nc1) {
    double MaxFred = 0;
    unsigned int nc2 = 0;
    for (cluster_iterator c2 = begincluster(); c2 != endcluster(); ++c2, ++nc2) {
      if (c1 != c2) {
        double Fred = averageDist[nc1] + averageDist[nc2];
        Fred /= Cdist_->CentroidDist( c1->Cent(), c2->Cent() );
        if (Fred > MaxFred)
          MaxFred = Fred;
      }
    }
    DBITotal += MaxFred;
  }
  DBITotal /= (double)clusters_.size();
  if (outfile.IsOpen()) outfile.Printf("#DBI: %f\n", DBITotal);
  return DBITotal;
}

/** The pseudo-F statistic is another measure of clustering goodness. It is 
  * intended to capture the 'tightness' of clusters, and is in essence a ratio
  * of the mean sum of squares between groups to the mean sum of squares within
  * group (Lattin et al., 2003: 291) High values are good. Generally, one 
  * selects a cluster-count that gives a peak in the pseudo-f statistic (or 
  * pSF, for short).
  * Formula: A/B, where A = (T - P)/(G-1), and B = P / (n-G). Here n is the 
  * number of points, G is the number of clusters, T is the total distance from
  * the all-data centroid, and P is the sum (for all clusters) of the distances
  * from the cluster centroid.
  * NOTE: To use this, cluster centroids should be fully up-to-date.
  * NOTE: This calc differs slightly from PTRAJ in that real centroids are used
  *       instead of representative structures.
  */
double ClusterList::ComputePseudoF(CpptrajFile& outfile) const {
  // Calculation makes no sense with fewer than 2 clusters.
  if (Nclusters() < 2) {
    mprintf("Warning: Fewer than 2 clusters. Not calculating pseudo-F.\n");
    return 0.0;
  }

  // Form a cluster with all points to get a centroid. Use only frames that
  // are in clusters, i.e. ignore noise. Assumes all cluster centroids are
  // up to date.
  ClusterNode c_all;
  for (cluster_iterator C1 = begincluster(); C1 != endcluster(); ++C1)
  {
    for (ClusterNode::frame_iterator f1 = C1->beginframe(); f1 != C1->endframe(); ++f1)
      c_all.AddFrameToCluster( *f1 );
  }
  // Pseudo-F makes no sense if # clusters == # frames
  if (Nclusters() == c_all.Nframes()) {
    mprintf("Warning: Each frame is in a separate cluster. Not calculating pseudo-F.\n");
    return 0.0;
  }
  c_all.SortFrameList();
  c_all.CalculateCentroid( Cdist_ );

  // Loop over all clusters
  double gss = 0.0; // between-group sum of squares
  double wss = 0.0; // within-group sum of squares
  for (cluster_iterator C1 = begincluster(); C1 != endcluster(); ++C1)
  {
    for (ClusterNode::frame_iterator f1 = C1->beginframe(); f1 != C1->endframe(); ++f1)
    {
      double dist = Cdist_->FrameCentroidDist(*f1, c_all.Cent());
      gss += (dist * dist);
      dist = Cdist_->FrameCentroidDist(*f1, C1->Cent());
      wss += (dist * dist);
    }
  }
  double d_nclusters = (double)Nclusters();
  double d_ntotal = (double)c_all.Nframes();
  double num = (gss - wss) / (d_nclusters - 1.0);
  double den = wss / (d_ntotal - d_nclusters);
  if (den < Constants::SMALL)
    den = Constants::SMALL;
  double pseudof = num / den;
  if (debug_ > 0)
    mprintf("Pseudo-f: Total distance to centroid is %.4f\n"
            "Pseudo-f: Cluster distance to centroid is %.4f\n"
            "Pseudo-f: Numerator %.4f over denominator %.4f gives %.4f\n", 
            gss, wss, num, den, pseudof);
  if (outfile.IsOpen()) {
    outfile.Printf("#pSF: %f\n", pseudof);
    // This calculation taken directly from ptraj
    double SSRSST = pseudof*(d_nclusters-1)/(d_ntotal-d_nclusters+pseudof*(d_nclusters-1));
    outfile.Printf("#SSR/SST: %f\n", SSRSST);
  }

  return pseudof;
}

/** The cluster silhouette is a measure of how well each point fits within
  * a cluster. Values of 1 indicate the point is very similar to other points
  * in the cluster, i.e. it is well-clustered. Values of -1 indicate the point
  * is dissimilar and may fit better in a neighboring cluster. Values of 0
  * indicate the point is on a border between two clusters. 
  */
void ClusterList::CalcSilhouette(std::string const& prefix) const {
  mprintf("\tCalculating cluster/frame silhouette.\n");
  CpptrajFile Ffile, Cfile;
  if (Ffile.OpenWrite(prefix + ".frame.dat")) return;
  if (Cfile.OpenWrite(prefix + ".cluster.dat")) return;
  Cfile.Printf("%-8s %10s\n", "#Cluster", "<Si>");
  unsigned int idx = 0;
  for (cluster_iterator Ci = begincluster(); Ci != endcluster(); ++Ci)
  {
    Ffile.Printf("#C%-6i %10s\n", Ci->Num(), "Silhouette");
    double avg_si = 0.0;
    int ci_frames = 0;
    std::vector<double> SiVals;
    for (ClusterNode::frame_iterator f1 = Ci->beginframe(); f1 != Ci->endframe(); ++f1)
    {
      // Calculate the average dissimilarity of this frame with all other
      // points in this frames cluster.
      double ai = 0.0;
      int self_frames = 0;
      for (ClusterNode::frame_iterator f2 = Ci->beginframe(); f2 != Ci->endframe(); ++f2)
      {
        if (f1 != f2) {
          ai += Frame_Distance(*f1, *f2);
          ++self_frames;
        }
      }
      if (self_frames > 0)
        ai /= (double)self_frames;
      //mprintf("\t\tFrame %i cluster %i ai = %g\n", *f1+1, Ci->Num(), ai);
      // Determine lowest average dissimilarity of this frame with all
      // other clusters.
      double min_bi = DBL_MAX;
      for (cluster_iterator Cj = begincluster(); Cj != endcluster(); ++Cj)
      {
        if (Ci != Cj)
        {
          double bi = 0.0;
          // NOTE: ASSUMING NO EMPTY CLUSTERS
          for (ClusterNode::frame_iterator f2 = Cj->beginframe(); f2 != Cj->endframe(); ++f2)
            bi += Frame_Distance(*f1, *f2);
          bi /= (double)Cj->Nframes();
          //mprintf("\t\tFrame %i to cluster %i bi = %g\n", *f1 + 1, Cj->Num(), bi);
          if (bi < min_bi)
            min_bi = bi;
        }
      }
      double max_ai_bi = std::max( ai, min_bi );
      if (max_ai_bi == 0.0)
        mprinterr("Error: Divide by zero in silhouette calculation for frame %i\n", *f1 + 1);
      else {
        double si = (min_bi - ai) / max_ai_bi;
        SiVals.push_back( si );
        //Ffile.Printf("%8i %10.4f\n", *f1 + 1, si);
        avg_si += si;
        ++ci_frames;
      }
    }
    std::sort( SiVals.begin(), SiVals.end() );
    for (std::vector<double>::const_iterator it = SiVals.begin(); it != SiVals.end(); ++it, ++idx)
      Ffile.Printf("%8i %g\n", idx, *it);
    Ffile.Printf("\n");
    ++idx;
    if (ci_frames > 0)
      avg_si /= (double)ci_frames;
    Cfile.Printf("%8i %g\n", Ci->Num(), avg_si);
  }
}

// -----------------------------------------------------------------------------
void ClusterList::DrawGraph(bool use_z, DataSet* cnumvtime,
                            double min_tol, int max_iteration) const
{
  if (use_z)
    mprintf("\tCreating PDB of graph points based on pairwise distances. B-factor = cluster #.\n");
  else
    mprintf("\tAttempting to draw graph based on pairwise distances.\n");
  unsigned int nframes = FrameDistances().Nrows();
  std::vector<Vec3> Xarray; // Coords
  std::vector<Vec3> Farray; // Forces
  Xarray.reserve( nframes );
  Farray.assign( nframes, Vec3(0.0) );
  // Initialize coordinates. X and Y only.
  double zcoord = 0.0;
  double theta_deg = 0.0;
  double delta = 360.0 / (double)nframes;
  for (unsigned int n = 0; n != nframes; n++, theta_deg += delta) {
    double theta_rad = Constants::DEGRAD * theta_deg;
    if (use_z)
      zcoord = cos(theta_rad / 2.0);
    Xarray.push_back( Vec3(cos(theta_rad), sin(theta_rad), zcoord) );
  }
  // Write out initial graph
  if (debug_ > 0 && !use_z) {
    CpptrajFile graph0;
    if (graph0.OpenWrite("InitialGraph.dat")) return;
    for (std::vector<Vec3>::const_iterator XV = Xarray.begin();
                                           XV != Xarray.end(); ++XV)
      graph0.Printf("%g %g %u\n", (*XV)[0], (*XV)[1], XV - Xarray.begin() + 1);
    graph0.CloseFile();
  }
  // Degrees of freedom. If Z ever initialized needs to be 3N
  double deg_of_freedom = 2.0 * (double)nframes;
  if (use_z) deg_of_freedom += (double)nframes;
  double fnq = sqrt( deg_of_freedom );
  // Main loop for steepest descent
  const double Rk = 1.0;
  const double dxstm = 1.0E-5;
  const double crits = 1.0E-6;
  double rms = 1.0;
  double dxst = 0.1;
  double last_e = 0.0;
  int iteration = 0;
  mprintf("          \t%8s %12s %12s\n", " ", "ENE", "RMS");
  while (rms > min_tol && iteration < max_iteration) {
    double e_total = 0.0;
    unsigned int idx = 0; // Index into FrameDistances
    for (unsigned int f1 = 0; f1 != nframes; f1++)
    {
      for (unsigned int f2 = f1 + 1; f2 != nframes; f2++)
      {
        //double Req = FrameDistances().GetCdist(f1, f2);
        Vec3 V1_2 = Xarray[f1] - Xarray[f2];
        double r2 = V1_2.Magnitude2();
        double s = sqrt(r2);
        double r = 2.0 / s;
        double db = s - FrameDistances().GetElement(idx++);
        double df = Rk * db;
        double e = df * db;
        e_total += e;
        df *= r;
        // Apply force
        V1_2 *= df;
        Farray[f1] -= V1_2;
        Farray[f2] += V1_2;
      }
    }
    // Calculate the magnitude of the force vector.
    double sum = 0.0;
    for (std::vector<Vec3>::const_iterator FV = Farray.begin(); FV != Farray.end(); ++FV)
      sum += FV->Magnitude2();
    rms = sqrt( sum ) / fnq;
    // Adjust search step size
    if (dxst < crits) dxst = dxstm;
    dxst = dxst / 2.0;
    if (e_total < last_e) dxst = dxst * 2.4;
    double dxsth = dxst / sqrt( sum );
    last_e = e_total;
    // Update positions and reset force array.
    std::vector<Vec3>::iterator FV = Farray.begin();
    for (std::vector<Vec3>::iterator XV = Xarray.begin();
                                     XV != Xarray.end(); ++XV, ++FV)
    {
      *XV += (*FV * dxsth);
      *FV = 0.0;
    }
    // Write out current E.
    mprintf("Iteration:\t%8i %12.4E %12.4E\n", iteration, e_total, rms);
    iteration++;
  }
  // RMS error 
  unsigned int idx = 0; // Index into FrameDistances
  double sumdiff2 = 0.0;
  for (unsigned int f1 = 0; f1 != nframes; f1++)
  {
    for (unsigned int f2 = f1 + 1; f2 != nframes; f2++)
    {
      Vec3 V1_2 = Xarray[f1] - Xarray[f2];
      double r1_2 = sqrt( V1_2.Magnitude2() );
      double Req = FrameDistances().GetElement(idx);
      double diff = r1_2 - Req;
      sumdiff2 += (diff * diff);
      if (debug_ > 0)
        mprintf("\t\t%u to %u: D= %g  Eq= %g  Delta= %g\n",
                f1+1, f2+1, r1_2, Req, fabs(diff));
      ++idx;
    }
  }
  double rms_err = sqrt( sumdiff2 / (double)FrameDistances().Nelements() );
  mprintf("\tRMS error of final graph positions: %g\n", rms_err);
  // Write out final graph with cluster numbers.
  std::vector<int> Nums;
  Nums.reserve( nframes );
  if (cnumvtime != 0) {
    ClusterSieve::SievedFrames const& sievedFrames = FrameDistances().FramesToCluster();
    DataSet_1D const& CVT = static_cast<DataSet_1D const&>( *cnumvtime );
    for (unsigned int n = 0; n != nframes; n++)
      Nums.push_back( (int)CVT.Dval(sievedFrames[n]) );
  } else
    for (int n = 1; n <= (int)nframes; n++)
      Nums.push_back( n );
  if (!use_z) {
    CpptrajFile graph;
    if (graph.OpenWrite("DrawGraph.dat")) return;
    for (std::vector<Vec3>::const_iterator XV = Xarray.begin();
                                           XV != Xarray.end(); ++XV)
      graph.Printf("%g %g %i \"%u\"\n", (*XV)[0], (*XV)[1], 
                   Nums[XV - Xarray.begin()], XV - Xarray.begin() + 1);
    graph.CloseFile();
  } else {
    // Write out PDB with B-factors
    PDBfile pdbout;
    if (pdbout.OpenWrite("DrawGraph.pdb")) return;
    pdbout.WriteTITLE("Cluster points.");
    for (std::vector<Vec3>::const_iterator XV = Xarray.begin();
                                           XV != Xarray.end(); ++XV)
      pdbout.WriteCoord(PDBfile::HETATM, XV - Xarray.begin() + 1, "HE", "HE",
                        XV - Xarray.begin() + 1, (*XV)[0], (*XV)[1], (*XV)[2],
                        1.0, Nums[XV - Xarray.begin()], "HE", 0);
    pdbout.CloseFile();
  }
}
