#include <algorithm> // sort
#include "Output.h"
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
