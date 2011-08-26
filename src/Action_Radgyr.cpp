// Radius of Gyration 
#include "Action_Radgyr.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Radgyr::Radgyr() {
  //fprintf(stderr,"Radgyr Con\n");
  rog=NULL;
  rogmax=NULL;
  useMass = false;
  calcRogmax=true;
} 

// DESTRUCTOR
Radgyr::~Radgyr() { }

/* Radgyr::init()
 * Expected call: radgyr <name> <mask1> [out filename] [mass] [nomax]
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
  if (A->hasKey("nomax")) calcRogmax=false;

  // Get Masks
  mask1 = A->getNextMask();
  Mask1.SetMaskString(mask1);

  // Datasets to store radius of gyration and max
  // Also add datasets to data file list 
  rog = DSL->Add(DOUBLE, A->getNextString(),"RoG");
  if (rog==NULL) return 1;
  DFL->Add(rogFile,rog);
  if (calcRogmax) {
    rogmax = DSL->Add(DOUBLE, NULL, "RoGMax");
    if (rogmax == NULL) return 1; 
    DFL->Add(rogFile,rogmax);
  }

  mprintf("    RADGYR: Calculating for atoms in mask %s",Mask1.maskString);
  if (useMass)
    mprintf(" using mass weighting");
  mprintf(".\n");
  if (!calcRogmax)
    mprintf("            RoG max will not be stored.\n");

  return 0;
}

/* Radgyr::setup()
 * Set radius of gyration up for this parmtop. Get masks etc.
 * P is set in Action::Setup
 */
int Radgyr::setup() {

  if ( Mask1.SetupMask(P,debug) ) return 1;
  if (Mask1.None()) {
    mprintf("    Error: Radgyr::setup: Mask contains 0 atoms.\n");
    return 1;
  }

  return 0;  
}

/* Radgyr::action()
 */
int Radgyr::action() {
  double Rog, max;

  Rog = F->RADGYR(&Mask1, useMass, &max);

  rog->Add(currentFrame, &Rog);
  if (calcRogmax)
    rogmax->Add(currentFrame, &max);

  //fprintf(outfile,"%10i %10.4lf %10.4lf\n",currentFrame,Rog,max);
  
  return 0;
} 

