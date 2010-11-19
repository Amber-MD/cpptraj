// Distance
#include <cstdio>
#include <cstdlib>
#include "Action_Distance.h"

// CONSTRUCTOR
Distance::Distance() {
  //fprintf(stderr,"Distance Con\n");
  dist=NULL;
  noimage=false;
  useMass=true;
  imageType=NONE;
  //currentType=DISTANCE;
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
    fprintf(stdout,"    Error: Distance::init: Requires 2 masks\n");
    return 1;
  }
  Mask1.SetMaskString(mask1);
  Mask2.SetMaskString(mask2);

  // Dataset to store distances
  dist = DSL->Add(DOUBLE, A->getNextString());
  if (dist==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(distanceFile,dist);


  //dist->Info();
  fprintf(stdout,"    DISTANCE: %s to %s",Mask1.maskString, Mask2.maskString);
  if (noimage) 
    fprintf(stdout,", non-imaged");
  else
    fprintf(stdout,", imaged");
  if (useMass) 
    fprintf(stdout,", center of mass");
  else
    fprintf(stdout,", geometric center");
  fprintf(stdout,".\n");

  return 0;
}

int Distance::setup() {
//  int m1atoms, m2atoms;

  // DEBUG - Test of AtomMask
  if ( Mask1.SetupMask(P,debug) ) return 1;
  if ( Mask2.SetupMask(P,debug) ) return 1;

  if (Mask1.None() || Mask2.None()) {
    fprintf(stdout,"    Error: Distance::setup: One or both masks have no atoms.\n");
    return 1;
  }
  // Set mass from parm file
  //if (!geom) Mass=P->mass;

  if (P->mass==NULL && useMass) {
    fprintf(stdout,"    Warning: Distance::setup: Mass for this parm is NULL.\n");
    fprintf(stdout,"             Geometric center of mass will be used.\n");
    useMass=false;
  }

  // Figure out imaging - check box based on prmtop box
  // NOTE: Should box be figured out from read-in coords?
  imageType = NONE;
  if (!noimage) {
    switch (P->ifbox) {
      case 0: 
        fprintf(stdout,"    Warning: Distance::setup: ");
        fprintf(stdout," Imaging specified but no box information in prmtop %s\n",P->parmName);
        fprintf(stdout,"             Disabling imaging.\n");
        noimage = true;
        break;
      case 1: imageType = ORTHO;    break;
      case 2: imageType = NONORTHO; break;
      default:
        fprintf(stdout,"    Error: Distance::setup: Unrecognized box type (%i) in %s\n.",
                P->ifbox, P->parmName);
        return 1;
    }
  }
        
  return 0;  
}

/*
 * Distance::action()
 */
int Distance::action() {
  double D;
  double ucell[9], recip[9];

  switch (imageType) {
    case NONE:     D = F->DIST(&Mask1, &Mask2, useMass); break;
    case ORTHO:    D = F->DIST_ImageOrtho(&Mask1, &Mask2, useMass); break;
    case NONORTHO: 
      F->BoxToRecip(ucell, recip);
      D = F->DIST_ImageNonOrtho(&Mask1, &Mask2, useMass, ucell, recip);
      break;
  }

  dist->Add(currentFrame, &D);

  //fprintf(outfile,"%10i %10.4lf\n",currentFrame,D);
  
  return 0;
} 


