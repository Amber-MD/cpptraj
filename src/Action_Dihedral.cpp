// Dihedral
#include "Action_Dihedral.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG

// CONSTRUCTOR
Dihedral::Dihedral() {
  //fprintf(stderr,"Dihedral Con\n");
  dih=NULL;
  useMass=false;
} 

// Dihedral::init()
/** Expected call: dihedral <name> <mask1> <mask2> <mask3> <mask4> [out filename]
  *                         [mass]
  */
int Dihedral::init() {
  char *mask1, *mask2, *mask3, *mask4;
  char *dihedralFile;

  // Get keywords
  dihedralFile = actionArgs.getKeyString("out",NULL);
  useMass = actionArgs.hasKey("mass");

  // Get Masks
  mask1 = actionArgs.getNextMask();
  mask2 = actionArgs.getNextMask();
  mask3 = actionArgs.getNextMask();
  mask4 = actionArgs.getNextMask();
  if (mask1==NULL || mask2==NULL || mask3==NULL || mask4==NULL) {
    mprinterr("Error: Dihedral::init: Requires 4 masks\n");
    return 1;
  }
  M1.SetMaskString(mask1);
  M2.SetMaskString(mask2);
  M3.SetMaskString(mask3);
  M4.SetMaskString(mask4);

  // Setup dataset
  dih = DSL->Add(DOUBLE, actionArgs.getNextString(),"Dih");
  if (dih==NULL) return 1;
  // Add dataset to datafile list
  DFL->Add(dihedralFile,dih);

  mprintf("    DIHEDRAL: [%s]-[%s]-[%s]-[%s]\n", M1.MaskString(),M2.MaskString(),
          M3.MaskString(), M4.MaskString());
  if (useMass)
    mprintf("              Using center of mass of atoms in masks.\n");

  return 0;
}

// Dihedral::setup
int Dihedral::setup() {

  if (currentParm->SetupIntegerMask(M1, activeReference)) return 1;
  if (currentParm->SetupIntegerMask(M2, activeReference)) return 1;
  if (currentParm->SetupIntegerMask(M3, activeReference)) return 1;
  if (currentParm->SetupIntegerMask(M4, activeReference)) return 1;
  mprintf("\t%s (%i atoms)\n",M1.MaskString(),M1.Nselected);
  mprintf("\t%s (%i atoms)\n",M2.MaskString(),M2.Nselected);
  mprintf("\t%s (%i atoms)\n",M3.MaskString(),M3.Nselected);
  mprintf("\t%s (%i atoms)\n",M4.MaskString(),M4.Nselected);
  if ( M1.None() || M2.None() || M3.None() || M4.None() ) {
    mprintf("    Error: Dihedral::setup: One or more masks have no atoms.\n");
    return 1;
  }

  return 0;  
}

// Dihedral::action()
int Dihedral::action() {
  double D;

  D=currentFrame->DIHEDRAL(&M1,&M2,&M3,&M4,useMass);

  D *= RADDEG;

  dih->Add(frameNum, &D);

  //fprintf(outfile,"%10i %10.4lf\n",frameNum,D);
  
  return 0;
} 

