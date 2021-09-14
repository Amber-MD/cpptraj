#include "DrawGraph.h"
#include "Cframes.h"
#include "MetricArray.h"
#include "../Constants.h"
#include "../CpptrajStdio.h"
#include "../DataSet_1D.h"
#include "../PDBfile.h"
#include "../Vec3.h"
#include <vector>
#include <cmath>

void Cpptraj::Cluster::DrawGraph(Cframes const& framesToCluster, MetricArray& pmatrix,
                                 GraphType graphTypeIn, DataSet* cnumvtime,
                                 double min_tol, int max_iteration, int debug)
{
  bool use_z;
  if (graphTypeIn == NO_DRAWGRAPH)
    return;
  else if (graphTypeIn == THREED)
    use_z = true;
  else
    use_z = false;
  if (use_z)
    mprintf("\tCreating PDB of graph points based on pairwise distances. B-factor = cluster #.\n");
  else
    mprintf("\tAttempting to draw graph based on pairwise distances.\n");
  unsigned int nframes = framesToCluster.size();
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
  if (debug > 0 && !use_z) {
    CpptrajFile graph0;
    if (graph0.OpenWrite("InitialGraph.dat")) return;
    for (std::vector<Vec3>::const_iterator XV = Xarray.begin();
                                           XV != Xarray.end(); ++XV)
      graph0.Printf("%g %g %li\n", (*XV)[0], (*XV)[1], XV - Xarray.begin() + 1);
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
    //unsigned int idx = 0; // Index into FrameDistances
    for (unsigned int f1 = 0; f1 != nframes; f1++)
    {
      for (unsigned int f2 = f1 + 1; f2 != nframes; f2++)
      {
        //double Req = FrameDistances().GetCdist(f1, f2);
        Vec3 V1_2 = Xarray[f1] - Xarray[f2];
        double r2 = V1_2.Magnitude2();
        double s = sqrt(r2);
        double r = 2.0 / s;
        double db = s - pmatrix.Frame_Distance(framesToCluster[f1], framesToCluster[f2]);//FrameDistances().GetElement(idx++);
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
  //unsigned int idx = 0; // Index into FrameDistances
  double sumdiff2 = 0.0;
  for (unsigned int f1 = 0; f1 != nframes; f1++)
  {
    for (unsigned int f2 = f1 + 1; f2 != nframes; f2++)
    {
      Vec3 V1_2 = Xarray[f1] - Xarray[f2];
      double r1_2 = sqrt( V1_2.Magnitude2() );
      double Req = pmatrix.Frame_Distance(framesToCluster[f1], framesToCluster[f2]);//FrameDistances().GetElement(idx);
      double diff = r1_2 - Req;
      sumdiff2 += (diff * diff);
      if (debug > 0)
        mprintf("\t\t%u to %u: D= %g  Eq= %g  Delta= %g\n",
                f1+1, f2+1, r1_2, Req, fabs(diff));
      //++idx;
    }
  }
  unsigned int Nelements = nframes * (nframes-1)/2;
  double rms_err = sqrt( sumdiff2 / (double)Nelements );
  mprintf("\tRMS error of final graph positions: %g\n", rms_err);
  // Write out final graph with cluster numbers.
  std::vector<int> Nums;
  Nums.reserve( nframes );
  if (cnumvtime != 0) {
    //ClusterSieve::SievedFrames const& sievedFrames = FrameDistances().FramesToCluster();
    DataSet_1D const& CVT = static_cast<DataSet_1D const&>( *cnumvtime );
    for (unsigned int n = 0; n != nframes; n++)
      Nums.push_back( (int)CVT.Dval(framesToCluster[n]) );
  } else
    for (int n = 1; n <= (int)nframes; n++)
      Nums.push_back( n );
  if (!use_z) {
    CpptrajFile graph;
    if (graph.OpenWrite("DrawGraph.dat")) return;
    for (std::vector<Vec3>::const_iterator XV = Xarray.begin();
                                           XV != Xarray.end(); ++XV)
      graph.Printf("%g %g %i \"%li\"\n", (*XV)[0], (*XV)[1], 
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
