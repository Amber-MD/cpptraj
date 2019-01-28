#include <cmath> // sqrt
#include <algorithm> // sort, max
#include "Output.h"
#include "../Matrix.h"
// -----------------------------------------------------------------------------
// ClusterList::PrintClustersToFile()
/** Print list of clusters in a style similar to ptraj; each cluster is
  * given a line maxframes characters long, with X for each frame that is
  * in the clusters and . for all other frames. Also print out the
  * representative frame numbers.
  */
void Cpptraj::Cluster::Output::PrintClustersToFile(CpptrajFile& outfile,
                                                   List const& clusters,
//std::string const& filename,
                                                   Algorithm const& algorithmIn,
                                                   Metric* metricIn, int sieve,
                                                   Cframes const& sievedFrames)
{
  //CpptrajFile outfile;
  std::string buffer;
  
  /*if ( outfile.OpenWrite(filename) ) {
    mprinterr("Error: PrintClustersToFile: Could not set up file %s\n",
              filename.c_str());
    return;
  }*/
  outfile.Printf("#Clustering: %i clusters %i frames\n",
                 clusters.Nclusters(), metricIn->Ntotal());
  // DBI
  std::vector<double> averageDist;
  double DBITotal = clusters.ComputeDBI( averageDist, metricIn );
  std::vector<double>::const_iterator avgd = averageDist.begin();
  for (List::cluster_iterator C = clusters.begincluster();
                              C != clusters.endcluster();
                            ++C, ++avgd)
    outfile.Printf("#Cluster %i has average-distance-to-centroid %f\n", C->Num(), *avgd);
  outfile.Printf("#DBI: %f\n", DBITotal);
  // Pseudo-F
  double SSRSST = 0.0;
  double pseudof = clusters.ComputePseudoF( SSRSST, metricIn );
  outfile.Printf("#pSF: %f\n", pseudof);
  outfile.Printf("#SSR/SST: %f\n", SSRSST);

  // Call internal info routine.
  algorithmIn.Results( outfile );
  // Print noise frames.
  if (clusters.Noise().size() > 0) {
    outfile.Printf("#NOISE_FRAMES:");
    unsigned int numNoise = 0;
    for (Cframes::const_iterator frm = clusters.Noise().begin();
                                 frm != clusters.Noise().end(); ++frm)
    {
      outfile.Printf(" %i", *frm+1);
      ++numNoise;
    }
    outfile.Printf("\n");
    outfile.Printf("#Number_of_noise_frames: %u\n", numNoise);
  }
  // Do not print trajectory stuff if no filename given (i.e. STDOUT output)
  if (!outfile.IsStream()) {
    for (List::cluster_iterator C1 = clusters.begincluster(); C1 != clusters.endcluster(); ++C1)
    {
      buffer.clear();
      buffer.resize(metricIn->Ntotal(), '.');
      for (Node::frame_iterator f1 = C1->beginframe(); f1 != C1->endframe(); ++f1)
        buffer[ *f1 ] = 'X';
      buffer += '\n';
      outfile.Write((void*)buffer.c_str(), buffer.size());
    }
  }
  // Print representative frame numbers
  outfile.Printf("#Representative frames:");
  for (List::cluster_iterator C1 = clusters.begincluster(); C1 != clusters.endcluster(); ++C1)
    if (C1->BestReps().size() < 2)
      outfile.Printf(" %i", C1->BestRepFrame()+1);
    else {
      outfile.Printf(" {");
      for (Node::RepPairArray::const_iterator rep = C1->BestReps().begin();
                                              rep != C1->BestReps().end(); ++rep)
        outfile.Printf(" %i %g", rep->first+1, rep->second);
      outfile.Printf(" }");
    }
  outfile.Printf("\n");
  // Print sieve info if present
  if (sieve != 1) {
    if (sieve < -1) {
      outfile.Printf("#Sieve value: %i (random)\n#Sieved frames:", -sieve);
      for (Cframes::const_iterator sfrm = sievedFrames.begin();
                                   sfrm != sievedFrames.end(); ++sfrm)
        outfile.Printf(" %i", *sfrm + 1);
      outfile.Printf("\n");
    } else
      outfile.Printf("#Sieve value: %i\n", sieve);
  }
  outfile.CloseFile();
}

/// For sorting cluster frame silhouettes by silhouette value.
struct sort_by_sil_val {
  typedef Cpptraj::Cluster::Node::SilPair Bpair;
  inline bool operator()(Bpair const& p0, Bpair const& p1)
  {
    if (p0.second == p1.second)
      return (p0.first < p1.first);
    else
      return (p0.second < p1.second);
  }
};

/** Print cluster silhouette frame values, sorted by silhouette. */
int Cpptraj::Cluster::Output::PrintSilhouetteFrames(CpptrajFile& Ffile, List const& clusters)
{
  // TODO different ways of writing out cluster frame silhouettes
  unsigned int idx = 0;
  for (List::cluster_iterator Ci = clusters.begincluster();
                              Ci != clusters.endcluster(); ++Ci, ++idx)
  {
    Ffile.Printf("#C%-6i %10s\n", Ci->Num(), "Silhouette");
    Node::SilPairArray spaTemp = Ci->FrameSilhouettes();
    std::sort( spaTemp.begin(), spaTemp.end(), sort_by_sil_val() );
    for (Node::SilPairArray::const_iterator it = spaTemp.begin();
                                            it != spaTemp.end(); ++it, ++idx)
      Ffile.Printf("%8u %g\n", idx, it->second);
    Ffile.Printf("\n");
  }
  return 0;
}

/** Print average cluster silhouette values. */
int Cpptraj::Cluster::Output::PrintSilhouettes(CpptrajFile& Cfile, List const& clusters)
{
  // TODO is it ok to assume clusters are in order?
  Cfile.Printf("%-8s %10s\n", "#Cluster", "<Si>");
  for (List::cluster_iterator Ci = clusters.begincluster();
                              Ci != clusters.endcluster(); ++Ci)
    Cfile.Printf("%8i %g\n", Ci->Num(), Ci->Silhouette());
  return 0;
}

/** Quick pass through clusters to determine max width of cluster names. */
unsigned int Cpptraj::Cluster::Output::DetermineNameWidth(List const& clusters)
{
  unsigned int nWidth = 0;
  for (List::cluster_iterator node = clusters.begincluster();
                              node != clusters.endcluster(); ++node)
    nWidth = std::max(nWidth, (unsigned int)node->Cname().size());
  return nWidth;
}

/** Print a summary of clusters. */
int Cpptraj::Cluster::Output::Summary(CpptrajFile& outfile, List const& clusters,
                                      Algorithm const& algorithm,
                                      PairwiseMatrix const& pmatrix, bool includeSieved,
                                      Cframes const& sievedOut)
{
  double fmax = (double)pmatrix.DistMetric().Ntotal();
  //if (FrameDistances().SieveValue() != 1 && !includeSieveInAvg)
  //  mprintf("Warning: Within cluster average distance (AvgDist) does not include sieved frames.\n");
  outfile.Printf("%-8s %8s %8s %8s %8s","#Cluster","Frames","Frac", "AvgDist","Stdev");
  if (!clusters.empty() && clusters.front().BestReps().size() > 1) {
    int nBestReps = clusters.front().BestReps().size();
    for (int i = 0; i != nBestReps; i++)
      outfile.Printf(" %8s %8s", "Rep", "RepScore");
  } else
    outfile.Printf(" %8s", "Centroid");
  outfile.Printf(" %8s", "AvgCDist");
  unsigned int nWidth = DetermineNameWidth(clusters);
  if (nWidth > 0) {
    if (nWidth < 8) nWidth = 8;
    outfile.Printf(" %*s %8s", nWidth, "Name", "RMS");
  }
  outfile.Printf("\n");
  //Timer t_fdist; // DEBUG
  //Timer t_cdist; // DEBUG
  //t_cdist.Start();
  // Calculate distances between clusters.
  Matrix<double> cluster_distances;
  cluster_distances.resize( 0, clusters.Nclusters() );
  for (List::cluster_iterator c1 = clusters.begincluster(); c1 != clusters.endcluster(); ++c1)
    for (List::cluster_iterator c2 = c1; c2 != clusters.endcluster(); ++c2)
      if (c2 != c1)
        cluster_distances.addElement( algorithm.ClusterDistance( *c1, *c2, pmatrix,
                                                                includeSieved, sievedOut ) );
  //t_cdist.Stop();

  unsigned int idx1 = 0;
  for (List::cluster_iterator node = clusters.begincluster();
                              node != clusters.endcluster(); ++node, ++idx1)
  {
    // Calculate the average distance of this cluster to every other cluster.
    //t_cdist.Start();
    double avgclusterdist = 0.0;
    if (clusters.Nclusters() > 1) {
      unsigned int idx2 = 0;
      for (List::cluster_iterator node2 = clusters.begincluster();
                                  node2 != clusters.endcluster(); ++node2, ++idx2)
      {
        if (node != node2)
          avgclusterdist += cluster_distances.element(idx1, idx2);
      }
      avgclusterdist /= (double)(clusters.Nclusters() - 1);
      //mprintf("CLUSTER %i avgclusterdist= %g\n", node->Num(), avgclusterdist);
    }
    //t_cdist.Stop();
    // Since there may be a lot of frames do not calculate SD from the
    // mean (which requires either storing distances or two double loops), 
    // instead use SD = sqrt( (SUM[x^2] - ((SUM[x])^2)/N)/N )
    //t_fdist.Start();
    double internalAvg = 0.0;
    double internalSD = 0.0;
    unsigned int Nelements = 0;
    if (node->Nframes() > 1) {
      // Calculate average distance between all frames in this cluster.
      if (includeSieved) {
        for (Node::frame_iterator f1 = node->beginframe(); f1 != node->endframe(); ++f1) {
          for (Node::frame_iterator f2 = f1 + 1; f2 != node->endframe(); ++f2) {
            double dist = pmatrix.Frame_Distance(*f1, *f2);
            internalAvg += dist;
            internalSD += (dist * dist);
            ++Nelements;
          }
        }
      } else {
        for (Node::frame_iterator f1 = node->beginframe(); f1 != node->endframe(); ++f1) {
          if (!sievedOut.HasFrame( *f1 )) {
            for (Node::frame_iterator f2 = f1 + 1; f2 != node->endframe(); ++f2) {
              if (!sievedOut.HasFrame( *f2 )) {
                double dist = pmatrix.GetFdist(*f1, *f2);
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
      //t_fdist.Stop();
    }
    // OUTPUT - TODO handle case when clusters dont have same number best reps
    outfile.Printf("%8i %8i %8.3f %8.3f %8.3f",
                   node->Num(), node->Nframes(), (double)node->Nframes()/fmax,
                   internalAvg, internalSD);
    if (node->BestReps().size() < 2)
      outfile.Printf(" %8i", node->BestRepFrame()+1);
    else {
      for (Node::RepPairArray::const_iterator rep = node->BestReps().begin();
                                              rep != node->BestReps().end(); ++rep)
        outfile.Printf(" %8i %8.3f", rep->first+1, rep->second);
    }
    outfile.Printf(" %8.3f", avgclusterdist);
    if (nWidth > 0)
      outfile.Printf(" %*s %8.3f", nWidth, node->Cname().c_str(), node->RefRms());
    outfile.Printf("\n");
  } // END loop over clusters
  //t_cdist.WriteTiming(1, "Between-cluster distance calc.");
  //t_fdist.WriteTiming(1, "Within-cluster distance calc.");
  return 0;
}

