// Radius of Gyration 
#include <cstdio>
#include <cstdlib>
#include "Action_Radgyr.h"

// CONSTRUCTOR
Radgyr::Radgyr() {
  //fprintf(stderr,"Radgyr Con\n");
  rog=NULL;
  useMass = false;
} 

// DESTRUCTOR
Radgyr::~Radgyr() { }

/*
 * Radgyr::init()
 * Expected call: radgyr <name> <mask1> [out filename] [mass]
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int Radgyr::init() {
  char *mask1;
  char *rogFile;

  // Get keywords
  rogFile = A->getKeyString("out",NULL);
  useMass = A->hasKey("mass");

  // Get Masks
  mask1 = A->getNextMask();
  Mask1.SetMaskString(mask1);

  // Dataset to store radius of gyration 
  rog = DSL->Add(DOUBLE, A->getNextString());
  if (rog==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(rogFile,rog);

  fprintf(stdout,"    RADGYR: Calculating for atoms in mask %s",Mask1.maskString);
  if (useMass)
    fprintf(stdout," using mass weighting");
  fprintf(stdout,".\n");

  return 0;
}

/*
 * Radgyr::setup()
 * Set radius of gyration up for this parmtop. Get masks etc.
 * P is set in Action::Setup
 */
int Radgyr::setup() {

  if ( Mask1.SetupMask(P,debug) ) return 1;
  if (Mask1.None()) {
    fprintf(stdout,"    Error: Radgyr::setup: Mask contains 0 atoms.\n");
    return 1;
  }

  return 0;  
}

/*
 * Radgyr::action()
 */
int Radgyr::action() {
  double Rog;

  //Ang=F->ANGLE(&Mask1,&Mask2,&Mask3);

  rog->Add(currentFrame, &Rog);

  //fprintf(outfile,"%10i %10.4lf\n",currentFrame,Rog);
  
  return 0;
} 


