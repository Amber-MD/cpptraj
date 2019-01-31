#include <limits> // double max 
#include <algorithm> // sort
#include "Algorithm_DBscan.h"
#include "../CpptrajStdio.h"
#include "../ProgressBar.h"
#include "../StringRoutines.h" // integerToString
#ifdef _OPENMP
#  include <omp.h>
#endif
#include "../DataSet_MatrixDbl.h" // DEBUG
#include "../DataFile.h" // DEBUG

// CONSTRUCTOR
Cpptraj::Cluster::Algorithm_DBscan::Algorithm_DBscan() :
  Algorithm(DBSCAN),
  minPoints_(-1),
  epsilon_(-1.0),
  sieveToCentroid_(true)
{}

// Algorithm_DBscan::Help()
void Cpptraj::Cluster::Algorithm_DBscan::Help() {
  mprintf("\t[dbscan minpoints <n> epsilon <e> [sievetoframe] [kdist <k> [kfile <prefix>]]]\n");
}

// Cluster_DBSCAN::SetupCluster()
int Cpptraj::Cluster::Algorithm_DBscan::Setup(ArgList& analyzeArgs) {
  kdist_.SetRange(analyzeArgs.GetStringKey("kdist"));
  if (kdist_.Empty()) {
    minPoints_ = analyzeArgs.getKeyInt("minpoints", -1);
    if (minPoints_ < 1) {
      mprinterr("Error: DBSCAN requires minimum # of points to be set and >= 1\n"
                "Error: Use 'minpoints <N>'\n");
      return 1;
    }
    epsilon_ = analyzeArgs.getKeyDouble("epsilon", -1.0);
    if (epsilon_ <= 0.0) {
      mprinterr("Error: DBSCAN requires epsilon to be set and > 0.0\n"
                "Error: Use 'epsilon <e>'\n");
      return 1;
    }
    sieveToCentroid_ = !analyzeArgs.hasKey("sievetoframe");
  } else {
    k_prefix_ = analyzeArgs.GetStringKey("kfile");
    if (!k_prefix_.empty() && k_prefix_.at(k_prefix_.size()-1) != '/')
      k_prefix_ += '/';
  }
  return 0;
}

// Cluster_DBSCAN::ClusteringInfo()
void Cpptraj::Cluster::Algorithm_DBscan::Info() const {
  mprintf("\tDBSCAN:\n");
  if (!kdist_.Empty()) {
    mprintf("\t\tOnly calculating Kdist graph for K=%s\n", kdist_.RangeArg());
    if (!k_prefix_.empty()) mprintf("\t\tKdist file prefix: %s\n", k_prefix_.c_str());
  } else {
    mprintf("\t\tMinimum pts to form cluster= %i\n", minPoints_);
    mprintf("\t\tCluster distance criterion= %.3f\n", epsilon_);
    if (sieveToCentroid_)
      mprintf("\t\tSieved frames will be added back solely based on their\n"
              "\t\t  closeness to cluster centroids.\n"
              "\t\t  (This option is less accurate but faster.)\n");
    else
      mprintf("\t\tSieved frames will only be added back if they are within\n"
              "\t\t  %.3f of a frame in an existing cluster.\n"
              "\t\t  (This option is more accurate and will identify sieved\n"
              "\t\t  frames as noise but is slower.)\n", epsilon_);
  }
}

// Cluster_DBSCAN::Cluster()
int Cpptraj::Cluster::Algorithm_DBscan::DoClustering(List& clusters,
                                                     Cframes const& framesToCluster,
                                                     PairwiseMatrix const& pmatrix)
{
  // Check if only need to calculate Kdist function(s)
  if (!kdist_.Empty()) {
    if (kdist_.Size() == 1)
      ComputeKdist( kdist_.Front(), framesToCluster, pmatrix );
    else
      ComputeKdistMap( kdist_, framesToCluster, pmatrix );
    return 0;
  }
  // Actual clustering
  unsigned int nPtsToCluster = framesToCluster.size();
  // SetOfPoints is UNCLASSIFIED
  Status_.assign( nPtsToCluster, UNCLASSIFIED );
  int ClusterId = 0;
  // If incoming clusters are not empty, set status of frames as appropriate.
  if (!clusters.empty()) {
    // Find highest cluster num
    for (List::cluster_iterator node = clusters.begincluster();
                                node != clusters.endcluster(); ++node)
      if (node->Num() > ClusterId) ClusterId = node->Num();
    ClusterId++;
    // See if any frames to be clustered are already in clusters.
    for (Cframes::const_iterator it = framesToCluster.begin(); it != framesToCluster.end(); ++it)
    {
      for (List::cluster_iterator node = clusters.begincluster();
                                  node != clusters.endcluster(); ++node)
      {
        if (node->HasFrame( *it ))
          Status_[it-framesToCluster.begin()] = node->Num();
      }
    }
  }
  ProgressBar progress( nPtsToCluster );
  for (unsigned int idx = 0; idx != nPtsToCluster; idx++)
  {
    progress.Update(idx);
    //Point := SetOfPoints.get(i);
    //IF Point.ClId = UNCLASSIFIED THEN
    if ( Status_[idx] == UNCLASSIFIED )
    {
      //IF ExpandCluster(SetOfPoints, Point, ClusterId, Eps, MinPts) THEN
      if (ExpandCluster(framesToCluster, pmatrix, idx, ClusterId))
      //ClusterId := nextId(ClusterId)
        ClusterId++;
    }
  }
  // Create clusters based on point statuses
  if (ClusterId > 0) {
    std::vector<Cframes> C0( ClusterId );
    // Make sure to store actual frame #s
    for (unsigned int idx = 0; idx != Status_.size(); idx++)
    {
      int cnum = Status_[idx];
      if (cnum == UNCLASSIFIED)
        mprintf("Warning: point %u was unclassified.\n", idx);
      else if (cnum == NOISE)
        clusters.AddNoise( framesToCluster[idx] );
      else
        C0[ cnum ].push_back( framesToCluster[idx] );
    }
    // Add clusters. 
    for (unsigned int cnum = 0; cnum != C0.size(); cnum++)
      clusters.AddCluster( Node(pmatrix.MetricPtr(), C0[cnum], cnum) );
  }
  return 0;
}

// Cluster_DBSCAN::ExpandCluster()
bool Cpptraj::Cluster::Algorithm_DBscan::ExpandCluster(Cframes const& framesToCluster,
                                                       PairwiseMatrix const& pmatrix,
                                                       unsigned int point, int ClusterId)
{
  //seeds:=SetOfPoints.regionQuery(Point,Eps);
  RegionQuery(seeds_, framesToCluster, pmatrix, point);

  //IF seeds.size<MinPts THEN // no core point
  if ((int)seeds_.size() < minPoints_)
  {
    //SetOfPoint.changeClId(Point,NOISE);
    Status_[point] = NOISE;
    //RETURN False;
    return false;
  }
  else
  {
    // all points in seeds are density-reachable from Point
    //SetOfPoints.changeClIds(seeds,ClId);
    Status_[point] = ClusterId;
    for (Iarray::const_iterator pt = seeds_.begin(); pt != seeds_.end(); ++pt)
      Status_[*pt] = ClusterId;
    //seeds.delete(Point);
    //WHILE seeds <> Empty DO
    unsigned int endIdx = seeds_.size();
    for (unsigned int idx = 0; idx < endIdx; idx++)
    {
      //currentP := seeds.first();
      int otherpoint = seeds_[idx];
      //result := SetOfPoints.regionQuery(currentP, Eps);
      RegionQuery(result_, framesToCluster, pmatrix, otherpoint);
      //IF result.size >= MinPts THEN
      if ( (int)result_.size() >= minPoints_ )
      {
        //FOR i FROM 1 TO result.size DO
        //  resultP := result.get(i);
        //  IF resultP.ClId IN {UNCLASSIFIED, NOISE} THEN
        //    IF resultP.ClId = UNCLASSIFIED THEN
        //      seeds.append(resultP);
        //    END IF;
        //    SetOfPoints.changeClId(resultP,ClId);
        //  END IF; // UNCLASSIFIED or NOISE
        //END FOR;
        for (Iarray::const_iterator rt = result_.begin(); rt != result_.end(); ++rt)
        {
          if (Status_[*rt] == UNCLASSIFIED || Status_[*rt] == NOISE)
          {
            if (Status_[*rt] == UNCLASSIFIED)
            {
              seeds_.push_back( *rt );
              endIdx = seeds_.size();
            }
            Status_[*rt] = ClusterId;
          }
        }
      }
      //END IF; // result.size >= MinPts
      //seeds.delete(currentP);
    }
    //END WHILE; // seeds <> Empty
    return true;
  }
}

// Cluster_DBSCAN::RegionQuery()
void Cpptraj::Cluster::Algorithm_DBscan::RegionQuery(Iarray& NeighborPts,
                                                     Cframes const& framesToCluster,
                                                     PairwiseMatrix const& pmatrix,
                                                     int point) const
{
  NeighborPts.clear();
  // point and otherpoint are indices, not frame #s
  int f1 = framesToCluster[ point ];
  for (int otherpoint = 0; otherpoint < (int)Status_.size(); ++otherpoint)
  {
    if (point != otherpoint) {
      int f2 = framesToCluster[ otherpoint ];
      if ( pmatrix.Frame_Distance(f1, f2) < epsilon_ )
        NeighborPts.push_back( otherpoint );
    }
  }
}

// Cluster_DBSCAN::ClusterResults()
void Cpptraj::Cluster::Algorithm_DBscan::Results(CpptrajFile& outfile) const {
  outfile.Printf("#Algorithm: DBSCAN minpoints %i epsilon %g sieveToCentroid %i\n",
                 minPoints_, epsilon_, (int)sieveToCentroid_);
/* TODO implement elsewhere
  // List the number of noise points.
  outfile.Printf("#NOISE_FRAMES:");
  unsigned int numNoise = 0;
  for (unsigned int idx = 0; idx != Status_.size(); ++idx)
  {
    if ( Status_[idx] == NOISE ) {
      outfile.Printf(" %u", FrameDistances().FramesToCluster()[idx]+1);
      ++numNoise;
    }
  }
  outfile.Printf("\n");
  outfile.Printf("#Number_of_noise_frames: %u\n", numNoise); 
*/
}

/*
// Cluster_DBSCAN::AddSievedFrames()
void Cpptraj::Cluster::Algorithm_DBscan::AddSievedFrames() {
  // NOTE: All cluster centroids must be up to date!
  if (sieveToCentroid_)
    mprintf("\tRestoring sieved frames if within %.3f of cluster centroid.\n", epsilon_);
  else
    mprintf("\tRestoring sieved frames if within %.3f of frame in nearest cluster.\n",
            epsilon_);
  // Vars allocated here in case of OpenMP
  int n_sieved_noise = 0;
  int Nsieved = 0;
  int frame;
  int nframes = (int)FrameDistances().OriginalNframes();
  ParallelProgress progress( nframes );
  // Need a temporary array to hold which frame belongs to which cluster. 
  // Otherwise we could be comparoing sieved frames to other sieved frames.
  std::vector<cluster_it> frameToCluster( nframes, clusters_.end() );
  // For OMP, every other thread will need its own Cdist.
  ClusterDist* MyCdist = Cdist_;
# ifdef _OPENMP
# pragma omp parallel private(MyCdist, frame) firstprivate(progress) reduction(+ : Nsieved, n_sieved_noise)
  {
  int mythread = omp_get_thread_num();
  progress.SetThread( mythread );
  if (mythread == 0) {
    mprintf("\tParallelizing calculation with %i threads\n", omp_get_num_threads());
    MyCdist = Cdist_;
  } else
    MyCdist = Cdist_->Copy();
# pragma omp for schedule(dynamic)
# endif
  for (frame = 0; frame < nframes; ++frame) {
    progress.Update( frame );
    if (FrameDistances().FrameWasSieved(frame)) {
      // Which clusters centroid is closest to this frame?
      double mindist = std::numeric_limits<double>::max();
      cluster_it minNode = clusters_.end();
      for (cluster_it Cnode = clusters_.begin(); Cnode != clusters_.end(); ++Cnode) {
        double dist = MyCdist->FrameCentroidDist(frame, Cnode->Cent());
        if (dist < mindist) {
          mindist = dist;
          minNode = Cnode;
        }
      }
      bool goodFrame = false;
      if ( mindist < epsilon_ )
        // Frame is already within epsilon, accept.
        goodFrame = true;
      else if ( !sieveToCentroid_ ) {
        // Check if any frames in the cluster are closer than epsilon to sieved frame.
        for (int cidx=0; cidx < minNode->Nframes(); cidx++)
        {
          if ( MyCdist->FrameDist(frame, minNode->ClusterFrame(cidx)) < epsilon_ )
          {
            goodFrame = true;
            break;
          }
        }
      }
      // Add sieved frame to the closest cluster if closest distance is
      // less than epsilon.
      ++Nsieved;
      if ( goodFrame )
        frameToCluster[frame] = minNode;
      else
        n_sieved_noise++;
    }
  } // END loop over frames
# ifdef _OPENMP
  if (mythread > 0)
    delete MyCdist;
  } // END pragma omp parallel
# endif
  progress.Finish();
  // Now actually add sieved frames to their appropriate clusters
  for (frame = 0; frame < nframes; frame++)
    if (frameToCluster[frame] != clusters_.end())
      frameToCluster[frame]->AddFrameToCluster( frame );
  mprintf("\t%i of %i sieved frames were discarded as noise.\n", 
          n_sieved_noise, Nsieved);
}
*/

/** For each point p, calculate function Kdist(p) which is the distance of
  * the Kth nearest point to p.
  */
void Cpptraj::Cluster::Algorithm_DBscan::ComputeKdist( int Kval,
                                                       Cframes const& FramesToCluster,
                                                       PairwiseMatrix const& pmatrix )
const
{
  std::vector<double> dists;
  std::vector<double> Kdist;
  dists.reserve( FramesToCluster.size() ); 
  Kdist.reserve( FramesToCluster.size() );
  std::string outfilename = k_prefix_ + "Kdist." + integerToString(Kval) + ".dat";
  mprintf("\tDBSCAN: Calculating Kdist(%i), output to %s\n", Kval, outfilename.c_str());
  for (std::vector<int>::const_iterator point = FramesToCluster.begin();
                                        point != FramesToCluster.end();
                                        ++point)
  {
    // Store distances from this point
    dists.clear();
    for (std::vector<int>::const_iterator otherpoint = FramesToCluster.begin();
                                          otherpoint != FramesToCluster.end();
                                          ++otherpoint)
      dists.push_back( pmatrix.Frame_Distance(*point, *otherpoint) );
    // Sort distances - first dist should always be 0
    std::sort(dists.begin(), dists.end());
    Kdist.push_back( dists[Kval] );
  }
  std::sort( Kdist.begin(), Kdist.end() );
  CpptrajFile Outfile;
  Outfile.OpenWrite(outfilename);
  Outfile.Printf("%-8s %1i%-11s\n", "#Point", Kval,"-dist");
  // Write out largest to smallest
  unsigned int ik = 0;
  for (std::vector<double>::reverse_iterator k = Kdist.rbegin(); 
                                             k != Kdist.rend(); ++k, ++ik)
    Outfile.Printf("%8u %12.4f\n", ik, *k);
  Outfile.CloseFile();
}

// Cluster_DBSCAN::ComputeKdistMap()
void Cpptraj::Cluster::Algorithm_DBscan::ComputeKdistMap( Range const& Kvals, 
                                                          Cframes const& FramesToCluster,
                                                          PairwiseMatrix const& pmatrix )
const
{
  int pt1_idx, pt2_idx, d_idx, point;
  mprintf("\tCalculating Kdist map for %s\n", Kvals.RangeArg());
  double* kdist_array; // Store distance of pt1 to every other point.
  int nframes = (int)FramesToCluster.size();
  // Ensure all Kdist points are within proper range
  Range::const_iterator kval;
  for (kval = Kvals.begin(); kval != Kvals.end(); ++kval)
    if (*kval < 1 || *kval >= nframes) {
      mprinterr("Error: Kdist value %i is out of range (1 <= Kdist < %i)\n",
                 *kval, nframes);
      return;
    }
  int nvals = (int)Kvals.Size();
  double** KMAP; // KMAP[i] has the ith nearest point for each point.
  KMAP = new double*[ nvals ];
  for (int i = 0; i != nvals; i++)
    KMAP[i] = new double[ nframes ];
  ParallelProgress progress( nframes );
# ifdef _OPENMP
# pragma omp parallel private(pt1_idx, pt2_idx, d_idx, kval, point, kdist_array) firstprivate(progress)
  {
  progress.SetThread( omp_get_thread_num() );
#endif
  kdist_array = new double[ nframes ];
# ifdef _OPENMP
# pragma omp for
# endif
  for (pt1_idx = 0; pt1_idx < nframes; pt1_idx++) // X
  {
    progress.Update( pt1_idx );
    point = FramesToCluster[pt1_idx];
    d_idx = 0;
    // Store distances from pt1 to pt2
    for (pt2_idx = 0; pt2_idx != nframes; pt2_idx++)
      kdist_array[d_idx++] = pmatrix.Frame_Distance(point, FramesToCluster[pt2_idx]);
    // Sort distances; will be smallest to largest
    std::sort( kdist_array, kdist_array + nframes );
    // Save the distance of specified nearest neighbors to this point.
    d_idx = 0;
    for (kval = Kvals.begin(); kval != Kvals.end(); ++kval) // Y
      KMAP[d_idx++][pt1_idx] = kdist_array[ *kval ];
  }
  delete[] kdist_array;
# ifdef _OPENMP
  } // END omp parallel
# endif
  progress.Finish();
  // Sort all of the individual kdist plots, smallest to largest.
  for (int i = 0; i != nvals; i++)
    std::sort(KMAP[i], KMAP[i] + nframes);
  // Save in matrix, largest to smallest.
  DataSet_MatrixDbl kmatrix;
  kmatrix.Allocate2D( FramesToCluster.size(), Kvals.Size() );
  for (int y = 0; y != nvals; y++) {
    for (int x = nframes - 1; x != -1; x--)
      kmatrix.AddElement( KMAP[y][x] );
    delete[] KMAP[y];
  }
  delete[] KMAP;
  // Write matrix to file
  DataFile outfile;
  ArgList outargs("usemap");
  outfile.SetupDatafile(k_prefix_ + "Kmatrix.gnu", outargs, debug_);
  outfile.AddDataSet( (DataSet*)&kmatrix );
  outfile.WriteDataOut();
  // Write out the largest and smallest values for each K.
  // This means for each value of K the point with the furthest Kth-nearest
  // neighbor etc.
  CpptrajFile maxfile;
  if (maxfile.OpenWrite(k_prefix_ + "Kmatrix.max.dat")) return;
  maxfile.Printf("%-12s %12s %12s\n", "#Kval", "MaxD", "MinD");
  d_idx = 0;
  for (kval = Kvals.begin(); kval != Kvals.end(); ++kval, d_idx++)
    maxfile.Printf("%12i %12g %12g\n", *kval, kmatrix.GetElement(0, d_idx),
                   kmatrix.GetElement(nframes-1, d_idx));
  maxfile.CloseFile();
}
