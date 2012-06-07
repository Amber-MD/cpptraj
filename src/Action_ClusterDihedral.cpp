#include <cstdio> // sscanf
#include <algorithm> // sort
#include "Action_ClusterDihedral.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG

// CONSTRUCTOR
Action_ClusterDihedral::Action_ClusterDihedral() :
  phibins_(0),
  psibins_(0),
  CUT_(0),
  lastframe_(0),
  dcparm_(0),
  CVT_(0)
{}

// Action_ClusterDihedral::ReadDihedrals()
int Action_ClusterDihedral::ReadDihedrals(std::string const& fname) {
  CpptrajFile infile;
  char buffer[256];
  int a1, a2, a3, a4, bins;

  if ( infile.OpenRead( fname ) ) return 1;
  mprintf("\tReading dihedral information from %s\n", fname.c_str());
  while (infile.Gets(buffer, 256)==0) {
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

// Action_ClusterDihedral::init()
/** Usage: clusterdihedral [phibins <N>] [psibins <M>] [out <outfile>]
  *                        [framefile <framefile>] [clusterinfo <infofile>]
  *                        [clustervtime <cvtfile>] [cut <CUT>] 
  *                        [dihedralfile <dfile> | <mask>]
  */
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
  std::string cvtfile = actionArgs.GetStringKey("clustervtime");
  // Input dihedral file or scan mask
  std::string dihedralIn = actionArgs.GetStringKey("dihedralfile");
  if (!dihedralIn.empty()) {
    if ( ReadDihedrals( dihedralIn ) != 0) return 1;
  } else {
    mask_.SetMaskString( actionArgs.getNextMask() );
  }

  // CVT dataset
  if (!cvtfile.empty()) {
    CVT_ = DSL->Add(DataSet::INT, actionArgs.getNextString(), "DCVT");
    if (CVT_ == NULL) return 1;
    DFL->Add(cvtfile.c_str(), CVT_);
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
  if (!cvtfile.empty())
    mprintf("\tNumber of clusters v time will be output to %s\n", cvtfile.c_str());
  
  return 0;
}

// Action_ClusterDihedral::setup()
int Action_ClusterDihedral::setup() {
  // Currently setup can only be performed based on first prmtop
  if (dcparm_!=0) {
    mprintf("Warning: clusterdihedral is only setup based on the first prmtop\n");
    mprintf("Warning: read in. Skipping setup for this prmtop.\n");
    return 0;
  }
  // Set up backbone dihedral angles if none were read
  if (DCmasks_.empty()) {
    // Setup mask
    if ( currentParm->SetupIntegerMask( mask_ ) ) return 1;
    if ( mask_.None()) {
      mprinterr("Error clusterdihedral: No atoms selected by mask [%s]\n", mask_.MaskString());
      return 1;
    }
    // NOTE: This code relies on selection having contiguous residues!
    int C1 = -1;
    int N2 = -1;
    int CA = -1;
    int C2 = -1;
    for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
    {
      if ( C2 > -1 ) {
        // If we have already found the last C in phi dihedral, this N is 
        // the last atom in psi dihedral - store both.
        if ( (*currentParm)[*atom].Name() == "N   " ) {
          DCmasks_.push_back( DCmask(C1, N2, CA, C2,    phibins_) ); // PHI
          DCmasks_.push_back( DCmask(N2, CA, C2, *atom, psibins_) ); // PSI
          if (debug > 0)
            mprintf("DIHEDRAL PAIR FOUND: C1= %i, N2= %i, CA= %i, C2= %i, N3= %li\n",
                    C1, N2, CA, C2, *atom);
          // Since the carbonyl C/amide N probably starts a new dihedral,
          // reset to those.
          C1 = C2;
          N2 = *atom;
          C2 = -1;
          CA = -1;
        }
      } else if ( C1 > -1 ) {
        // If we've already found the first carbonyl, look for other atoms
        // in the dihedral pair.
        if ( (*currentParm)[*atom].Name() == "N   " ) N2 = *atom;
        if ( (*currentParm)[*atom].Name() == "CA  " ) CA = *atom;
        if ( (*currentParm)[*atom].Name() == "C   " ) C2 = *atom;
      } else if ( (*currentParm)[*atom].Name() == "C   " ) C1 = *atom; // 1st carbon
    } // End loop over selected atoms
    mprintf("\tFound %zu dihedral angles.\n", DCmasks_.size());
  }

  if (DCmasks_.empty()) {
    mprinterr("Error: clusterdihedral: No dihedral angles defined.\n");
    return 1;
  }

  // Allocate space to hold Bin IDs each frame
  Bins_.resize( DCmasks_.size() );
  dcparm_ = currentParm;

  return 0;
}

// Action_ClusterDihedral::action()
int Action_ClusterDihedral::action() {
  // For each dihedral, calculate which bin it should go into and store bin#
  int bidx = 0;
  for (std::vector<DCmask>::iterator dih = DCmasks_.begin(); 
                                     dih != DCmasks_.end(); ++dih)
  {
    double PHI = currentFrame->DIHEDRAL( (*dih).A1(), (*dih).A2(), (*dih).A3(), (*dih).A4() );
    // NOTE: Torsion is in radians; should bins be converted to rads as well?
    PHI *= RADDEG;
    PHI += 180;
    double phistep = 360 / (double)(*dih).Bins();
    PHI /= phistep;
    int phibin = (int)PHI;
    Bins_[bidx++] = phibin;
  }
  // DEBUG - print bins
  //mprintf("[");
  //for (std::vector<int>::iterator bin = Bins_.begin(); bin != Bins_.end(); ++bin)
  //  mprintf("%3i,",*bin);
  //mprintf("]\n");
  // Now search for this bin combo in the DCarray
  std::vector<DCnode>::iterator DC = dcarray_.begin();
  for (; DC != dcarray_.end(); ++DC)
  {
    if ( (*DC).BinMatch( Bins_ ) ) break;
  }
  if (DC == dcarray_.end()) {
    // No match; create new bin combo and store frame num
    //mprintf("NEW DCNODE.\n");
    dcarray_.push_back( DCnode( Bins_, frameNum ) );
  } else {
    // Match; increment bin count and store frame num
    //mprintf("DCNODE ALREADY PRESENT.\n");
    (*DC).Increment();
    (*DC).AddFrame( frameNum );
  }
  // Store frame number
  lastframe_ = frameNum;
  // Update cvt dataset if specified
  if (CVT_ != 0) {
    int cvtdata = (int)dcarray_.size();
    CVT_->Add(frameNum, &cvtdata);
  }
  return 0;
}

// Action_ClusterDihedral::print()
void Action_ClusterDihedral::print() {
  // Setup output file
  CpptrajFile output;
  if (output.OpenWrite( outfile_ )) return;
  mprintf("\tPrinting Dihedral Clustering Results.\n");
  // Print bin information
  output.Printf("DIHEDRAL CLUSTER RESULTS");
  if (mask_.MaskStringSet())
    output.Printf(" for %s", mask_.MaskString());
  output.Printf("\n");
  long int num = 0;
  for (std::vector<DCmask>::iterator dih = DCmasks_.begin();
                                     dih != DCmasks_.end(); ++dih)
  {
    output.Printf("    %6li ", num++);
    output.Printf("%-s(%i)", (*dcparm_)[ (*dih).A1() ].c_str(), (*dih).A1() + 1);
    output.Printf("%-s(%i)", (*dcparm_)[ (*dih).A2() ].c_str(), (*dih).A2() + 1);
    output.Printf("%-s(%i)", (*dcparm_)[ (*dih).A3() ].c_str(), (*dih).A3() + 1);
    output.Printf("%-s(%i)", (*dcparm_)[ (*dih).A4() ].c_str(), (*dih).A4() + 1);
    output.Printf(" [Bins=%i]\n", (*dih).Bins());
  }
  output.Printf("%zu clusters.\n", dcarray_.size());

  // Sort array by count
  std::sort( dcarray_.begin(), dcarray_.end() );

  // Allocate space for storing cluster #s for each frame
  std::vector<long int> framecluster( lastframe_ + 1 );

  // Print sorted cluster array
  if (CUT_ > 0)
    output.Printf("Only printing clusters with pop > %i\n",CUT_);
  num = 0;
  for (std::vector<DCnode>::iterator DC = dcarray_.begin();
                                     DC != dcarray_.end(); ++DC)
  {
    if ( (*DC).Count() > CUT_ ) {
      //mprintf("DEBUG: Cluster %li has %i frames.\n", num, (*DC).NumFrames());
      output.Printf("Cluster %10li %10li [ ", num, (*DC).Count());
      for (DCnode::bin_it binid = (*DC).binbegin(); binid != (*DC).binend(); ++binid)
        output.Printf("%3i ", *binid);
      output.Printf(" ]\n");
      for (DCnode::frame_it frame = (*DC).framebegin(); frame != (*DC).frameend(); ++frame)
      {
        output.Printf("%i ", *frame + OUTPUTFRAMESHIFT);
        // store which cluster each frame belongs to. Not neccesary if user
        // didn't specify this option, but avoids a second loop if they did.
        framecluster[ *frame ] = num;
      }
      output.Printf("\n");
    }
    ++num;
  }
  output.CloseFile();

  // Print cluster for each frame
  if (!framefile_.empty()) {
    if (output.OpenWrite( framefile_ )) return;
    mprintf("\tPrinting cluster number for each frame.\n");
    num = 1;
    for (std::vector<long int>::iterator cnum = framecluster.begin(); 
                                         cnum != framecluster.end(); ++cnum)
    {
      // Frame, cluster num, cluster count
      output.Printf("%10li %10i %10li ", num++, *cnum, dcarray_[*cnum].Count());
      // Print binID
      for (DCnode::bin_it binid = dcarray_[*cnum].binbegin();
                          binid != dcarray_[*cnum].binend(); ++binid)
        output.Printf("%03i", *binid);
      output.Printf("\n");
    }
    output.CloseFile();
  }

  // Print cluster information file
  if (!infofile_.empty()) {
    if (output.OpenWrite( infofile_ )) return;
    mprintf("\tPrinting cluster information.\n");
    output.Printf("%zu\n", DCmasks_.size());
    for (std::vector<DCmask>::iterator dih = DCmasks_.begin();
                                       dih != DCmasks_.end(); ++dih)
    {
      output.Printf("%10i %10i %10i %10i %10i\n", (*dih).A1()+1, (*dih).A2()+1,
                    (*dih).A3()+1, (*dih).A4()+1, (*dih).Bins());
    }
    output.Printf("%zu\n", dcarray_.size());
    num = 0;
    for (std::vector<DCnode>::iterator DC = dcarray_.begin();
                                       DC != dcarray_.end(); ++DC)
    {
      output.Printf("%10li %10li ", num++, (*DC).Count());
      for (DCnode::bin_it binid = (*DC).binbegin(); binid != (*DC).binend(); ++binid)
        output.Printf("%03i", *binid);
      output.Printf("\n");
    }
    output.CloseFile();
  }

}
