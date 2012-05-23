#include <cstdio> // sscanf
#include "Action_ClusterDihedral.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_ClusterDihedral::Action_ClusterDihedral() :
  phibins_(0),
  psibins_(0),
  CUT_(0)
{}

int Action_ClusterDihedral::ReadDihedrals(std::string const& fname) {
  CpptrajFile infile;
  char buffer[256];
  int a1, a2, a3, a4, bins;

  if ( infile.OpenRead( fname ) ) return 1;
  mprintf("\tReading dihedral information from %s\n", fname.c_str());
  while (infile.IO->Gets(buffer, 256)==0) {
    // Expected line format: At#1 At#2 At#3 At#4 Bins
    // ATOM NUMBERS SHOULD START FROM 0!
    int nvals = sscanf(buffer, "%i %i %i %i %i", &a1, &a2, &a3, &a4, &bins);
    if (nvals < 5) {
      mprinterr("Error: Dihedral file %s: Expected 5 values, got %i\n", fname.c_str(), nvals);
      mprinterr("Error: Problem line: [%s]\n",buffer);
      return 1; // This should automatically close infile through destructor.
    }
    DCmasks_.push_back( DCmask(a1, a2, a3, a4, bins) );
    mprintf("\t\t(%i)-(%i)-(%i)-(%i) Bins=%i\n",a1,a2,a3,a4,bins);
  }
  mprintf("\tRead %zu dihedrals.\n", DCmasks_.size());
  infile.CloseFile();
  return 0;
}

int Action_ClusterDihedral::init() {
  // # of phi and psi bins
  phibins_ = actionArgs.getKeyInt("phibins", 10);
  psibins_ = actionArgs.getKeyInt("psibins", 10);
  if ( phibins_>360 || phibins_<=1 || psibins_>360 || psibins_<=1 ) {
    mprinterr("Error: clusterdihedral: phi or psi bins out of range 1 <= x < 360 (%i, %i)\n",
              phibins_, psibins_);
    return 1;
  }
  // Cluster pop cutoff
  CUT_ = actionArgs.getKeyInt("cut",0);
  // Output files
  outfile_ = actionArgs.GetStringKey("out");
  framefile_ = actionArgs.GetStringKey("framefile");
  infofile_ = actionArgs.GetStringKey("clusterinfo");
  cvtfile_ = actionArgs.GetStringKey("clustervtime");
  // Input dihedral file or scan mask
  std::string dihedralIn = actionArgs.GetStringKey("dihedralfile");
  if (!dihedralIn.empty()) {
    if ( ReadDihedrals( dihedralIn ) != 0) return 1;
  } else {
    mask_.SetMaskString( actionArgs.getNextMask() );
  }

  // INFO
  mprintf("    DIHEDRAL CLUSTERING:");
  if (DCmasks_.empty()) {
    mprintf(" PHI and PSI dihedrals will be scanned for using mask [%s]\n", mask_.MaskString());
    mprintf("\t\t# phi bins = %i   # psi bins = %i\n",phibins_,psibins_);
  } else {
    mprintf(" Clustering on %zu dihedral angles.\n", DCmasks_.size());
  }
  if (CUT_>0)
    mprintf("\tOnly clusters with population > %i will be printed.\n", CUT_);
  if (!framefile_.empty())
    mprintf("\tFrame-Cluster data will be output to %s\n", framefile_.c_str());
  if (!infofile_.empty())
    mprintf("\tCluster information (pop. & ID) will be output to %s\n", infofile_.c_str());
  if (!cvtfile_.empty())
    mprintf("\tNumber of clusters v time will be output to %s\n", cvtfile_.c_str());
  
  return 0;
}
