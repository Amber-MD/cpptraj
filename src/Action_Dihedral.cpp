// Action_Dihedral
#include "Action_Dihedral.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG
#include "TorsionRoutines.h"

// CONSTRUCTOR
Action_Dihedral::Action_Dihedral() :
  dih_(NULL)
{
  //fprintf(stderr,"Action_Dihedral Con\n");
  useMass_=false;
} 

// Action_Dihedral::init()
/** Expected call: dihedral <name> <mask1> <mask2> <mask3> <mask4> [out filename]
  *                         [mass]
  *                         [type {alpha|beta|gamma|delta|epsilon|zeta|chi|c2p
  *                                h1p  |phi |psi  |pchi}]
  */
int Action_Dihedral::init() {
  // Get keywords
  ArgList::ConstArg dihedralFile = actionArgs.getKeyString("out");
  useMass_ = actionArgs.hasKey("mass");
  DataSet::scalarType stype = DataSet::UNDEFINED;
  std::string stypename = actionArgs.GetStringKey("type");
  if      ( stypename == "alpha"   ) stype = DataSet::ALPHA;
  else if ( stypename == "beta"    ) stype = DataSet::BETA;
  else if ( stypename == "gamma"   ) stype = DataSet::GAMMA;
  else if ( stypename == "delta"   ) stype = DataSet::DELTA;
  else if ( stypename == "epsilon" ) stype = DataSet::EPSILON;
  else if ( stypename == "zeta"    ) stype = DataSet::ZETA;
  else if ( stypename == "chi"     ) stype = DataSet::CHI;
  else if ( stypename == "c2p"     ) stype = DataSet::C2P;
  else if ( stypename == "h1p"     ) stype = DataSet::H1P;
  else if ( stypename == "phi"     ) stype = DataSet::PHI;
  else if ( stypename == "psi"     ) stype = DataSet::PSI;
  else if ( stypename == "pchi"    ) stype = DataSet::PCHI;

  // Get Masks
  ArgList::ConstArg mask1 = actionArgs.getNextMask();
  ArgList::ConstArg mask2 = actionArgs.getNextMask();
  ArgList::ConstArg mask3 = actionArgs.getNextMask();
  ArgList::ConstArg mask4 = actionArgs.getNextMask();
  if (mask1==NULL || mask2==NULL || mask3==NULL || mask4==NULL) {
    mprinterr("Error: dihedral: Requires 4 masks\n");
    return 1;
  }
  M1_.SetMaskString(mask1);
  M2_.SetMaskString(mask2);
  M3_.SetMaskString(mask3);
  M4_.SetMaskString(mask4);

  // Setup dataset
  dih_ = DSL->Add(DataSet::DOUBLE, actionArgs.getNextString(),"Dih");
  if (dih_==NULL) return 1;
  dih_->SetScalar( DataSet::M_TORSION, stype );
  // Add dataset to datafile list
  DFL->Add(dihedralFile, dih_);

  mprintf("    DIHEDRAL: [%s]-[%s]-[%s]-[%s]\n", M1_.MaskString(), 
          M2_.MaskString(), M3_.MaskString(), M4_.MaskString());
  if (useMass_)
    mprintf("              Using center of mass of atoms in masks.\n");

  return 0;
}

// Action_Dihedral::setup
int Action_Dihedral::setup() {
  if (currentParm->SetupIntegerMask(M1_)) return 1;
  if (currentParm->SetupIntegerMask(M2_)) return 1;
  if (currentParm->SetupIntegerMask(M3_)) return 1;
  if (currentParm->SetupIntegerMask(M4_)) return 1;
  M1_.MaskInfo();
  M2_.MaskInfo();
  M3_.MaskInfo();
  M4_.MaskInfo();
  if ( M1_.None() || M2_.None() || M3_.None() || M4_.None() ) {
    mprintf("Warning: dihedral: One or more masks have no atoms.\n");
    return 1;
  }

  return 0;  
}

// Action_Dihedral::action()
int Action_Dihedral::action() {
  Vec3 a1, a2, a3, a4;

  if (useMass_) {
    a1 = currentFrame->VCenterOfMass( M1_ );
    a2 = currentFrame->VCenterOfMass( M2_ );
    a3 = currentFrame->VCenterOfMass( M3_ );
    a4 = currentFrame->VCenterOfMass( M4_ );
  } else {
    a1 = currentFrame->VGeometricCenter( M1_ );
    a2 = currentFrame->VGeometricCenter( M2_ );
    a3 = currentFrame->VGeometricCenter( M3_ );
    a4 = currentFrame->VGeometricCenter( M4_ );
  }
  double torsion = Torsion(a1.Dptr(), a2.Dptr(), a3.Dptr(), a4.Dptr());

  torsion *= RADDEG;

  dih_->Add(frameNum, &torsion);

  //fprintf(outfile,"%10i %10.4lf\n",frameNum,D);
  
  return 0;
} 

