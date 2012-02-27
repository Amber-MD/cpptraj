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

// Pucker::init()
/** Expected call: pucker <name> <mask1> <mask2> <mask3> <mask4> <mask5> out <filename>
  *                [range360] [amplitude] [altona | cremer] [offset <offset>]
  */
int Pucker::init() {
  char *mask1, *mask2, *mask3, *mask4, *mask5;
  char *puckerFile;

  // Get keywords
  puckerFile = actionArgs.getKeyString("out",NULL);
  if (actionArgs.hasKey("altona")) puckerMethod=0;
  else if (actionArgs.hasKey("cremer")) puckerMethod=1;
  if (actionArgs.hasKey("amplitude")) amplitude=true;
  offset = actionArgs.getKeyDouble("offset",0.0);
  if (actionArgs.hasKey("range360")) {
    puckermax=360.0;
    puckermin=0.0;
  }

  // Get Masks
  mask1 = actionArgs.getNextMask();
  mask2 = actionArgs.getNextMask();
  mask3 = actionArgs.getNextMask();
  mask4 = actionArgs.getNextMask();
  mask5 = actionArgs.getNextMask();
  if (mask1==NULL || mask2==NULL || mask3==NULL || mask4==NULL || mask5==NULL) {
    mprinterr("Error: Pucker::init: Requires 5 masks\n");
    return 1;
  }
  M1.SetMaskString(mask1);
  M2.SetMaskString(mask2);
  M3.SetMaskString(mask3);
  M4.SetMaskString(mask4);
  M5.SetMaskString(mask5);

  // Setup dataset
  puck = DSL->Add(DOUBLE, actionArgs.getNextString(),"Pucker");
  if (puck==NULL) return 1;
  // Add dataset to datafile list
  DFL->Add(puckerFile,puck);

  //dih->Info();
  mprintf("    PUCKER: [%s]-[%s]-[%s]-[%s]-[%s]\n", M1.MaskString(),M2.MaskString(),
          M3.MaskString(), M4.MaskString(), M5.MaskString());
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

// Pucker::setup
int Pucker::setup() {

  if ( currentParm->SetupIntegerMask( M1, activeReference) ) return 1;
  if ( currentParm->SetupIntegerMask( M2, activeReference) ) return 1;
  if ( currentParm->SetupIntegerMask( M3, activeReference) ) return 1;
  if ( currentParm->SetupIntegerMask( M4, activeReference) ) return 1;
  if ( currentParm->SetupIntegerMask( M5, activeReference) ) return 1;
  mprintf("\t%s (%i atoms)\n",M1.MaskString(),M1.Nselected);
  mprintf("\t%s (%i atoms)\n",M2.MaskString(),M2.Nselected);
  mprintf("\t%s (%i atoms)\n",M3.MaskString(),M3.Nselected);
  mprintf("\t%s (%i atoms)\n",M4.MaskString(),M4.Nselected);
  mprintf("\t%s (%i atoms)\n",M5.MaskString(),M5.Nselected);

  if ( M1.None() || M2.None() || M3.None() || M4.None() || M5.None() ) {
    mprintf("Warning: Pucker::setup: One or more masks have no atoms.\n");
    return 1;
  }

  return 0;  
}

// Pucker::action()
int Pucker::action() {
  double D;

  D=currentFrame->PUCKER(&M1,&M2,&M3,&M4,&M5,puckerMethod,amplitude,useMass);
  D *= RADDEG;

  // Deal with offset
  D += offset;

  // Wrap values > puckermax or < puckermin
  if      (D > puckermax) D -= 360.0;
  else if (D < puckermin) D += 360.0;

  puck->Add(frameNum, &D);

  //fprintf(outfile,"%10i %10.4lf\n",frameNum,D);
  
  return 0;
} 

