#include <limits> // double max
#include "Silhouette.h"
#include "List.h"
#include "MetricArray.h"
#include "Node.h"
#include "../CpptrajStdio.h"

using namespace Cpptraj::Cluster;

/** CONSTRUCTOR - set debug level */
Silhouette::Silhouette(int dbgIn) :
  debug_(dbgIn)
{}

/** The cluster silhouette is a measure of how well each point fits within
  * a cluster. Values of 1 indicate the point is very similar to other points
  * in the cluster, i.e. it is well-clustered. Values of -1 indicate the point
  * is dissimilar and may fit better in a neighboring cluster. Values of 0
  * indicate the point is on a border between two clusters. 
  */
int Silhouette::CalcSilhouette(List const& clusters,
                               MetricArray& metrics,
                               std::vector<bool> const& frameIsPresent,
                               bool includeSieved)
{
  clusterFrameSil_.clear();
  clusterAvgSil_.clear();
  // Loop over all clusters
  for (List::cluster_iterator Ci = clusters.begincluster(); Ci != clusters.endcluster(); ++Ci)
  {
    clusterFrameSil_.push_back( SilPairArray() );
    SilPairArray& SiVals = clusterFrameSil_.back();
    SiVals.reserve( Ci->Nframes() );
    double avg_si = 0.0;
    int ci_frames = 0;
    for (Node::frame_iterator f1 = Ci->beginframe(); f1 != Ci->endframe(); ++f1)
    {
      if (includeSieved || frameIsPresent[ *f1 ]) {
        // Calculate the average dissimilarity of this frame with all other
        // points in this frames cluster.
        double ai = 0.0;
        int self_frames = 0;
        if (includeSieved) {
          for (Node::frame_iterator f2 = Ci->beginframe(); f2 != Ci->endframe(); ++f2)
          {
            if (f1 != f2) {
              ai += metrics.Frame_Distance(*f1, *f2);
              ++self_frames;
            }
          }
        } else {
          for (Node::frame_iterator f2 = Ci->beginframe(); f2 != Ci->endframe(); ++f2)
          {
            if (f1 != f2 && frameIsPresent[ *f2 ]) {
              ai += metrics.Frame_Distance(*f1, *f2);
              ++self_frames;
            }
          }
        }
        if (self_frames > 0)
          ai /= (double)self_frames;
        //mprintf("\t\tFrame %i cluster %i ai = %g\n", *f1+1, Ci->Num(), ai);
        // Determine lowest average dissimilarity of this frame with all
        // other clusters.
        double min_bi = std::numeric_limits<double>::max();
        for (List::cluster_iterator Cj = clusters.begincluster(); Cj != clusters.endcluster(); ++Cj)
        {
          if (Ci != Cj)
          {
            double bi = 0.0;
            // NOTE: ASSUMING NO EMPTY CLUSTERS
            if (includeSieved) {
              for (Node::frame_iterator f2 = Cj->beginframe(); f2 != Cj->endframe(); ++f2)
                bi += metrics.Frame_Distance(*f1, *f2);
              bi /= (double)Cj->Nframes();
            } else {
              int cj_frames = 0;
              for (Node::frame_iterator f2 = Cj->beginframe(); f2 != Cj->endframe(); ++f2)
              {
                if (frameIsPresent[ *f2 ]) {
                  bi += metrics.Frame_Distance(*f1, *f2);
                  ++cj_frames;
                }
              }
              bi /= (double)cj_frames;
            }
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
          SiVals.push_back( SilPair(*f1, si) );
          avg_si += si;
          ++ci_frames;
        }
      } // END if frame should be calcd
    } // END loop over cluster frames
    //std::sort( SiVals.begin(), SiVals.end() );
    // DEBUG
    if (debug_ > 1) {
      mprintf("DEBUG: Cluster frame silhouette values for cluster %i\n", Ci->Num());
      for (SilPairArray::const_iterator it = SiVals.begin(); it != SiVals.end(); ++it)
        mprintf("\t%8i %g\n", it->first+1, it->second);
    }
    if (ci_frames > 0)
      avg_si /= (double)ci_frames;
    //mprintf("DEBUG: Cluster silhouette: %8i %g\n", Ci->Num(), avg_si);
    clusterAvgSil_.push_back( avg_si );
  } // END outer loop over clusters
  return 0;
}
