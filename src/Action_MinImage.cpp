#include <cmath>
#include <cfloat>
#include "Action_MinImage.h"
#include "CpptrajStdio.h"
#ifdef _OPENMP
#  include "omp.h"
#endif
#include "StringRoutines.h" // DEBUG
#include "Constants.h" // DEBUG

// CONSTRUCTOR
Action_MinImage::Action_MinImage() : 
  dist_(0),
  useMass_(true),
  calcUsingMask_(false)
{} 

void Action_MinImage::Help() {
  mprintf("\t[<name>] <mask1> <mask2> [out <filename>] [geom] [mask]\n"
          "  Calculate minimum non-self imaged distance between atoms in <mask1> and <mask2>\n");
}

// Action_MinImage::Init()
Action::RetType Action_MinImage::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get Keywords
  // Require imaging.
  image_.InitImaging( true );
  useMass_ = !(actionArgs.hasKey("geom"));
  calcUsingMask_ = actionArgs.hasKey("mask");
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  // Get Masks
  std::string mask1 = actionArgs.GetMaskNext();
  std::string mask2 = actionArgs.GetMaskNext();
  if (mask1.empty() || mask2.empty()) {
    mprinterr("Error: Requires 2 masks\n");
    return Action::ERR;
  }
  Mask1_.SetMaskString(mask1);
  Mask2_.SetMaskString(mask2);

  // Dataset to store distances
  dist_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "MID");
  if (dist_==0) return Action::ERR;
  dist_->SetScalar( DataSet::M_DISTANCE );
  // Add dataset to data file
  if (outfile != 0) outfile->AddSet( dist_ );
  int numthreads = 1;
# ifdef _OPENMP
#pragma omp parallel
{
  if (omp_get_thread_num()==0)
    numthreads = omp_get_num_threads();
}
# endif
  minDist_.resize( numthreads );

  mprintf("    MINIMAGE: Looking for closest approach of");
  if (calcUsingMask_) {
    mprintf(" center of mask %s\n\tto images of center of mask %s\n",
            Mask1_.MaskString(), Mask2_.MaskString());
    if (useMass_) 
      mprintf("\tUsing center of mass of masks.\n");
    else
      mprintf("\tUsing geometric center of masks.\n");
  } else
    mprintf(" atoms in %s\n\tto images of atoms in %s\n",
            Mask1_.MaskString(), Mask2_.MaskString());
  return Action::OK;
}

// Action_MinImage::Setup()
/** Determine what atoms each mask pertains to for the current parm file.
  */
Action::RetType Action_MinImage::Setup(Topology* currentParm, Topology** parmAddress) {
  if (currentParm->SetupIntegerMask( Mask1_ )) return Action::ERR;
  if (currentParm->SetupIntegerMask( Mask2_ )) return Action::ERR;
  mprintf("\t%s (%i atoms) to %s (%i atoms)\n",Mask1_.MaskString(), Mask1_.Nselected(),
          Mask2_.MaskString(),Mask2_.Nselected());
  if (Mask1_.None() || Mask2_.None()) {
    mprintf("Warning: One or both masks have no atoms.\n");
    return Action::ERR;
  }
  // Set up imaging info for this parm
  image_.SetupImaging( currentParm->BoxType() );
  if (!image_.ImagingEnabled()) {
    mprintf("Warning: Imaging cannot be performed for topology %s\n", currentParm->c_str());
    return Action::ERR;
  }

  return Action::OK;  
}
/*
static void WriteMatrix(Matrix_3x3 const& ucell, PDBfile& pdbout, const char* name, int res)
{
  Vec3 xvec = ucell.Row1();
  Vec3 yvec = ucell.Row2();
  Vec3 zvec = ucell.Row3();
  pdbout.WriteATOM("O", res, 0.0, 0.0, 0.0, name, 1.0);
  pdbout.WriteATOM("X", res, xvec[0], xvec[1], xvec[2], name, 1.0);
  mprintf("DEBUG: X={ %g %g %g }\n", xvec[0], xvec[1], xvec[2]);
  pdbout.WriteATOM("Y", res, yvec[0], yvec[1], yvec[2], name, 1.0);
  pdbout.WriteATOM("Z", res, zvec[0], zvec[1], zvec[2], name, 1.0);
}
*/

double Action_MinImage::MinNonSelfDist2(Vec3 const& a1, Vec3 const& a2) {
#ifndef NONSENSE
//  a1.Print("A1");
//  a2.Print("A2");
  Vec3 frac1 = recip_ * a1; // a1 in fractional coords
  // NOTE: Do not use 'floor' here since we want to keep a1/a2 in same ref frame
//  frac1.Print("A1 fractional coords");
  Vec3 T1 = ucell_.TransposeMult(frac1); // a1 back in Cartesian space
//  pdbout_.WriteATOM("O1", 6, T1[0], T1[1], T1[2], "O1", 1.0); // DEBUG
  Vec3 frac2 = recip_ * a2; // a2 in fractional coords 
//  frac2.Print("A2 fractional coords");
  // Floor
//  Vec3 floor1(frac1[0] - floor(frac1[0]), frac1[1] - floor(frac1[1]), frac1[2] - floor(frac1[2]));
//  Vec3 floor2(frac2[0] - floor(frac2[0]), frac2[1] - floor(frac2[1]), frac2[2] - floor(frac2[2]));
//  floor1.Print("Floor1 fractional coords");
//  floor2.Print("Floor2 fractional coords");
//  int ndist = 0; // DEBUG
  double minDist2 = DBL_MAX;
  for (int ix = -1; ix < 2; ix++) {
    for (int iy = -1; iy < 2; iy++) {
      for (int iz = -1; iz < 2; iz++) {
        if (ix != 0 || iy != 0 || iz != 0) { // Ignore a2 self
          Vec3 ixyz(ix, iy, iz);
//          Vec3 t1 = ucell_.TransposeMult(frac1 + ixyz); // DEBUG
          Vec3 t2 = ucell_.TransposeMult(frac2 + ixyz); // a2 image back in Cartesian space
          // Write out imaged coordinates
//          std::string name1 = "T1" + integerToString(ndist); // DEBUG
//          std::string name2 = "T2" + integerToString(ndist); // DEBUG
//          pdbout_.WriteATOM(name1.c_str(), ndist*2+7, t1[0], t1[1], t1[2], "T1", 1.0);
//          pdbout_.WriteATOM(name2.c_str(), ndist*2+8, t2[0], t2[1], t2[2], "T2", 1.0);
          // Distance from imaged a2 to original a1
          Vec3 dxyz = t2 - T1;
          double dist2 = dxyz.Magnitude2();
          minDist2 = std::min(minDist2, dist2);
        }
//        ndist++;
      }
    }
  }
#else
  // Project a1 in fractional coords axes.
  Vec3 f1 = recip_ * a1;
  // Ensure fractional coords.
  Vec3 frac1(f1[0] - floor(f1[0]), f1[1] - floor(f1[1]), f1[2] - floor(f1[2]));
  // Project back in original unit cell axes
  Vec3 T1 = ucell_.TransposeMult(frac1);
  pdbout_.WriteATOM("O1", 6, T1[0], T1[1], T1[2], "O1", 1.0); // DEBUG
  // Determine which projected a1 is closest to self.
/*  double minDist2 = DBL_MAX;
  Vec3 minT1(0.0);
  for (int ix = -1; ix < 2; ix++) {
    for (int iy = -1; iy < 2; iy++) {
      for (int iz = -1; iz < 2; iz++) {
        Vec3 ixyz(ix, iy, iz);
        // Project back in original unit cell axes
        Vec3 t1 = ucell_.TransposeMult(frac1 + ixyz);
        Vec3 dxyz = t1 - a1;
        double dist2 = dxyz.Magnitude2();
        if (dist2 < dxyz.Magnitude2()) {
          minT1 = t1;
          minDist2 = dist2;
        }
      }
    }
  }
  pdbout_.WriteATOM("T1", 6, minT1[0], minT1[1], minT1[2], "T1", 1.0); // DEBUG*/

  // Project a2 in fractional coords axes.
  Vec3 f2 = recip_ * a2;
  // Ensure fractional coords.
  Vec3 frac2(f2[0] - floor(f2[0]), f2[1] - floor(f2[1]), f2[2] - floor(f2[2]));
  int ndist = 0; // DEBUG
  double minDist2 = DBL_MAX;
  // Loop over all images of a2, ignoring self. Calc dist to T1.
  for (int ix = -1; ix < 2; ix++) {
    for (int iy = -1; iy < 2; iy++) {
      for (int iz = -1; iz < 2; iz++) {
        if (ix != 0 || iy != 0 || iz != 0) { // Ignore self
          Vec3 ixyz(ix, iy, iz);
          Vec3 t1 = ucell_.TransposeMult(frac1 + ixyz);
          Vec3 t2 = ucell_.TransposeMult(frac2 + ixyz);
          std::string name1 = "T1" + integerToString(ndist);
          std::string name2 = "T2" + integerToString(ndist); // DEBUG
//          pdbout_.WriteATOM(name2.c_str(), ndist+7, t2[0], t2[1], t2[2], "T2", 1.0); // DEBUG
          pdbout_.WriteATOM(name1.c_str(), ndist*2+7, t1[0], t1[1], t1[2], "T1", 1.0);
          pdbout_.WriteATOM(name2.c_str(), ndist*2+8, t2[0], t2[1], t2[2], "T2", 1.0);
          //Vec3 dxyz = ucell_.TransposeMult(frac2 + ixyz) - a1;
          Vec3 dxyz = t2 - T1;
//        std::string name1 = "O" + integerToString(ndist);
//        std::string name2 = "D" + integerToString(ndist);
//        std::string rname = "V" + integerToString(ndist);
//        pdbout.WriteATOM(name1.c_str(), ndist*2+7, a1[0], a1[1], a1[2], rname.c_str(), 1.0);
//        pdbout.WriteATOM(name2.c_str(), ndist*2+8, a1[0] + dxyz[0],
//                                                   a1[1] + dxyz[1],
//                                                   a2[2] + dxyz[2], rname.c_str(), 1.0);
          double dist2 = dxyz.Magnitude2();
          //minDist2 = std::min(minDist2, dist2);
          if (dist2 < minDist2) {
            minDist2 = dist2;
            minxyz_ = t2;
          }
//          mprintf("DEBUG: Map: %2i = {%2g %2g %2g} %g\n", ndist, ixyz[0], ixyz[1], ixyz[2],
//                  sqrt(dist2));
        }
        ndist++;
      }
    }
  }
#endif
//  mprintf("DEBUG: Minimum distance= %g\n", sqrt(minDist2));
  return minDist2;
}


// Action_MinImage::DoAction()
Action::RetType Action_MinImage::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress)
{
//  PDBfile pdbout;
//  pdbout_.OpenWrite("minimage.pdb");
  double Dist2;

  currentFrame->BoxCrd().ToRecip(ucell_, recip_);
/*
  // Test alternate recip calc. Cartesian->Frac
  recip_.Print("Recip");
  Matrix_3x3 rtest;
  double a = currentFrame->BoxCrd().BoxX();
  double b = currentFrame->BoxCrd().BoxY();
  double c = currentFrame->BoxCrd().BoxZ();
  double alpha = currentFrame->BoxCrd().Alpha() * Constants::DEGRAD;
  double beta  = currentFrame->BoxCrd().Beta()  * Constants::DEGRAD;
  double gamma = currentFrame->BoxCrd().Gamma() * Constants::DEGRAD;
  double volume = sqrt( 1 - cos(alpha)*cos(alpha) - cos(beta)*cos(beta)
                                 - cos(gamma)*cos(gamma)
                                 + 2*cos(alpha)*cos(beta)*cos(gamma) );
  volume *= a * b * c;
  rtest[0] = 1 / a;
  rtest[1] = -cos( gamma ) / (a * sin( gamma ));
  rtest[2] = cos( alpha ) * cos( gamma ) - cos( beta );
  rtest[3] = 0.0;
  rtest[4] = 1 / (b * sin( gamma ));
  rtest[5] = (cos(beta)*cos(gamma) - cos(alpha)) / (b * volume * sin(gamma));
  rtest[6] = 0.0;
  rtest[7] = 0.0;
  rtest[8] = (a * b * sin( gamma )) / volume;
  rtest.Print("Rtest");
  // Test alternate ucell calc. Frac->Cartesian 
  ucell_.Print("Ucell");
  Matrix_3x3 utest;
  utest[0] = a;
  utest[1] = b * cos( gamma );
  utest[2] = c * cos( beta );
  utest[3] = 0.0;
  utest[4] = b * sin( gamma );
  utest[5] = c * ((cos(alpha)-cos(beta)*cos(gamma)) / sin(gamma));
  utest[6] = 0.0;
  utest[7] = 0.0;
  utest[8] = volume / (a * b * sin(gamma));
  utest.Print("Utest");
*/
  if (calcUsingMask_) {
    // Use center of mask1 and mask2
    Vec3 a1, a2;
    if (useMass_) {
      a1 = currentFrame->VCenterOfMass( Mask1_ );
      a2 = currentFrame->VCenterOfMass( Mask2_ );
    } else {
      a1 = currentFrame->VGeometricCenter( Mask1_ );
      a2 = currentFrame->VGeometricCenter( Mask2_ );
    }

/*
  Dist = DIST2_NoSelfImageNonOrtho(a1, a2, ucell, recip);
  Dist = sqrt(Dist);
*/
  // Project A1 & A2
//  Vec3 f1 = recip * a1;
//  Vec3 f2 = recip * a2;
//  pdbout.WriteATOM("F1", 5, f1[0], f1[1], f1[2], "F1", 1.0);
//  pdbout.WriteATOM("F2", 6, f2[0], f2[1], f2[2], "F2", 1.0);
  // Fractional coords?
//  Vec3 frac1(f1[0] - floor(f1[0]), f1[1] - floor(f1[1]), f1[2] - floor(f1[2]));
//  Vec3 frac2(f2[0] - floor(f2[0]), f2[1] - floor(f2[1]), f2[2] - floor(f2[2]));
//  frac1.Print("frac1");
//  frac2.Print("frac2");

//  pdbout.WriteATOM("A1", 1, a1[0], a1[1], a1[2], "A1", 1.0);
//  pdbout.WriteATOM("A2", 2, a2[0], a2[1], a2[2], "A2", 1.0);
  // Unit cell parameters
//  WriteMatrix( ucell, pdbout, "UNT", 3);
//  WriteMatrix( recip, pdbout, "RCP", 4);

    Dist2 = MinNonSelfDist2(a1, a2);
    Dist2 = sqrt( Dist2 );
  } else {
    // Look at all atoms in mask1/mask2
    minDist_.assign(minDist_.size(), DBL_MAX);
    int m1end = Mask1_.Nselected();
    int m2end = Mask2_.Nselected();
    int m1, m2;
    int mythread = 0;
    Vec3 a1;
#   ifdef _OPENMP
#   pragma omp parallel private(m1, m2, mythread, Dist2) firstprivate(a1)
    {
    mythread = omp_get_thread_num();
#   pragma omp for
#   endif 
    for (m1 = 0; m1 < m1end; m1++)
    {
      a1 = Vec3(currentFrame->XYZ( Mask1_[m1] ));
      for (m2 = 0; m2 < m2end; m2++)
      {
        Dist2 = MinNonSelfDist2( a1, Vec3(currentFrame->XYZ(Mask2_[m2])) );
        if (Dist2 < minDist_[mythread]) {
          minDist_[mythread] = Dist2;
//          rprintf("DEBUG: New Min Dist: Atom %i to %i (%g)\n",
//                  Mask1_[m1]+1,Mask2_[m2]+1,sqrt(Dist2));
          //pdbout_.WriteHET(1, minxyz_[0], minxyz_[1], minxyz_[2]);
        }
      }
    }
#   ifdef _OPENMP
    } // END openmp pragma
#   endif
    // Find lowest minDist
    double globalMin = minDist_[0];
    for (unsigned int n = 1; n != minDist_.size(); n++)
      if (minDist_[n] < globalMin)
        globalMin = minDist_[n];
    Dist2 = sqrt( globalMin );
  }

  dist_->Add(frameNum, &Dist2);

  return Action::OK;
}
