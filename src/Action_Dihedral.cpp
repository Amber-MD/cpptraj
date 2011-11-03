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

// DESTRUCTOR
Dihedral::~Dihedral() { }

/* Dihedral::init()
 * Expected call: dihedral <name> <mask1> <mask2> <mask3> <mask4> [out filename]
 *                         [mass]
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
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
    mprintf("    Error: Dihedral::init: Requires 4 masks\n");
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

  mprintf("    DIHEDRAL: %s-%s-%s-%s\n", M1.maskString,M2.maskString,
          M3.maskString, M4.maskString);
  if (useMass)
    mprintf("              Using center of mass of atoms in masks.\n");

  return 0;
}

/* Dihedral::setup
 */
int Dihedral::setup() {

  if ( M1.SetupMask(currentParm,activeReference,debug) ) return 1;
  if ( M2.SetupMask(currentParm,activeReference,debug) ) return 1;
  if ( M3.SetupMask(currentParm,activeReference,debug) ) return 1;
  if ( M4.SetupMask(currentParm,activeReference,debug) ) return 1;
  if ( M1.None() || M2.None() || M3.None() || M4.None() ) {
    mprintf("    Error: Dihedral::setup: One or more masks have no atoms.\n");
    return 1;
  }

  return 0;  
}

/* Dihedral::action()
 */
int Dihedral::action() {
  double D;

  D=currentFrame->DIHEDRAL(&M1,&M2,&M3,&M4,useMass);

  D *= RADDEG;

  dih->Add(frameNum, &D);

  //fprintf(outfile,"%10i %10.4lf\n",frameNum,D);
  
  return 0;
} 


