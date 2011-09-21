// Pucker
#include "Action_Pucker.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG

// CONSTRUCTOR
Pucker::Pucker() {
  //fprintf(stderr,"Pucker Con\n");
  puck=NULL;
  puckerMethod=0;
  amplitude=false;
  useMass=true;
  offset=0;
  puckermin = -180.0;
  puckermax = 180.0;
} 

// DESTRUCTOR
Pucker::~Pucker() { }

/* Pucker::init()
 * Expected call: pucker <name> <mask1> <mask2> <mask3> <mask4> <mask5> out <filename>
 *                [range360] [amplitude] [altona | cremer] [offset <offset>]
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int Pucker::init() {
  char *mask1, *mask2, *mask3, *mask4, *mask5;
  char *puckerFile;

  // Get keywords
  puckerFile = A->getKeyString("out",NULL);
  if (A->hasKey("altona")) puckerMethod=0;
  else if (A->hasKey("cremer")) puckerMethod=1;
  if (A->hasKey("amplitude")) amplitude=true;
  offset = A->getKeyDouble("offset",0.0);
  if (A->hasKey("range360")) {
    puckermax=360.0;
    puckermin=0.0;
  }

  // Get Masks
  mask1 = A->getNextMask();
  mask2 = A->getNextMask();
  mask3 = A->getNextMask();
  mask4 = A->getNextMask();
  mask5 = A->getNextMask();
  if (mask1==NULL || mask2==NULL || mask3==NULL || mask4==NULL || mask5==NULL) {
    mprintf("    Error: Pucker::init: Requires 5 masks\n");
    return 1;
  }
  M1.SetMaskString(mask1);
  M2.SetMaskString(mask2);
  M3.SetMaskString(mask3);
  M4.SetMaskString(mask4);
  M5.SetMaskString(mask5);

  // Setup dataset
  puck = DSL->Add(DOUBLE, A->getNextString(),"Pucker");
  if (puck==NULL) return 1;
  // Add dataset to datafile list
  DFL->Add(puckerFile,puck);

  //dih->Info();
  mprintf("    PUCKER: %s-%s-%s-%s-%s\n", M1.maskString,M2.maskString,
          M3.maskString, M4.maskString, M5.maskString);
  if (puckerMethod==0) 
    mprintf("            Using Altona & Sundaralingam method.\n");
  else if (puckerMethod==1)
    mprintf("            Using Cremer & Pople method.\n");
  if (puckerFile!=NULL) 
    mprintf("            Data will be written to %s\n",puckerFile);
  if (amplitude)
    mprintf("            Amplitudes will be stored instead of psuedorotation.\n");
  if (offset!=0)
    mprintf("            Offset: %lf will be added to values.\n");
  mprintf  ("            Values will range from %.1lf to %.1lf\n",puckermin,puckermax);

  return 0;
}

/* Pucker::setup
 */
int Pucker::setup() {

  if ( M1.SetupMask(P,debug) ) return 1;
  if ( M2.SetupMask(P,debug) ) return 1;
  if ( M3.SetupMask(P,debug) ) return 1;
  if ( M4.SetupMask(P,debug) ) return 1;
  if ( M5.SetupMask(P,debug) ) return 1;
  if ( M1.None() || M2.None() || M3.None() || M4.None() || M5.None() ) {
    mprintf("    Error: Pucker::setup: One or more masks have no atoms.\n");
    return 1;
  }

  return 0;  
}

/* Pucker::action()
 */
int Pucker::action() {
  double D;

  D=F->PUCKER(&M1,&M2,&M3,&M4,&M5,puckerMethod,amplitude,useMass);
  D *= RADDEG;

  // Deal with offset
  D += offset;

  // Wrap values > puckermax or < puckermin
  if      (D > puckermax) D -= 360.0;
  else if (D < puckermin) D += 360.0;

  puck->Add(currentFrame, &D);

  //fprintf(outfile,"%10i %10.4lf\n",currentFrame,D);
  
  return 0;
} 

