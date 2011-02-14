// Distance
#include <cmath>
#include "Action_Distance.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Distance::Distance() {
  //fprintf(stderr,"Distance Con\n");
  dist=NULL;
  noimage=false;
  useMass=true;
  imageType=0;
} 

// DESTRUCTOR
Distance::~Distance() {
  //fprintf(stderr,"Distance Destructor.\n");
}

/*
 * Distance::init()
 * Expected call: distance <name> <mask1> <mask2> [out filename] [geom] [noimage]
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int Distance::init( ) {
  char *mask1, *mask2;
  char *distanceFile;

  // Get Keywords
  noimage = A->hasKey("noimage");
  useMass = !(A->hasKey("geom"));
  distanceFile = A->getKeyString("out",NULL);

  // Get Masks
  mask1 = A->getNextMask();
  mask2 = A->getNextMask();
  //fprintf(stdout,"    Mask 1: %s\n",mask1);
  //fprintf(stdout,"    Mask 2: %s\n",mask2);
  if (mask1==NULL || mask2==NULL) {
    mprintf("    Error: Distance::init: Requires 2 masks\n");
    return 1;
  }
  Mask1.SetMaskString(mask1);
  Mask2.SetMaskString(mask2);

  // Dataset to store distances
  dist = DSL->Add(DOUBLE, A->getNextString(),"Dis");
  if (dist==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(distanceFile,dist);


  //dist->Info();
  mprintf("    DISTANCE: %s to %s",Mask1.maskString, Mask2.maskString);
  if (noimage) 
    mprintf(", non-imaged");
  else
    mprintf(", imaged");
  if (useMass) 
    mprintf(", center of mass");
  else
    mprintf(", geometric center");
  mprintf(".\n");

  return 0;
}

int Distance::setup() {

  if ( Mask1.SetupMask(P,debug) ) return 1;
  if ( Mask2.SetupMask(P,debug) ) return 1;

  if (Mask1.None() || Mask2.None()) {
    mprintf("    Error: Distance::setup: One or both masks have no atoms.\n");
    return 1;
  }
  // Set mass from parm file
  //if (!geom) Mass=P->mass;

  if (P->mass==NULL && useMass) {
    mprintf("    Warning: Distance::setup: Mass for this parm is NULL.\n");
    mprintf("             Geometric center of mass will be used.\n");
    useMass=false;
  }

  // Check imaging - check box based on prmtop box
  imageType = 0;
  if (!noimage) {
    imageType = P->ifbox;
    if (P->ifbox==0 && debug>0) {
      mprintf("    Warning: No box info in %s, disabling imaging.\n",P->parmName);
    }
  }
        
  return 0;  
}

/*
 * Distance::action()
 */
int Distance::action() {
  double D, ucell[9], recip[9];

  if (imageType>0) F->BoxToRecip(ucell,recip);
  D = F->DIST2(&Mask1, &Mask2, useMass, imageType, ucell, recip);
  D = sqrt(D);

  dist->Add(currentFrame, &D);

  //fprintf(outfile,"%10i %10.4lf\n",currentFrame,D);
  
  return 0;
} 


