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

// Radgyr::init()
/** Expected call: radgyr <name> <mask1> [out filename] [mass] [nomax]
  */
int Radgyr::init() {
  char *mask1, *rogname;
  char *rogFile;

  // Get keywords
  rogFile = actionArgs.getKeyString("out",NULL);
  useMass = actionArgs.hasKey("mass");
  if (actionArgs.hasKey("nomax")) calcRogmax=false;

  // Get Masks
  mask1 = actionArgs.getNextMask();
  Mask1.SetMaskString(mask1);

  // Datasets to store radius of gyration and max
  // Also add datasets to data file list
  rogname = actionArgs.getNextString();
  rog = DSL->Add(DOUBLE, rogname, "RoG");
  if (rog==NULL) return 1;
  DFL->Add(rogFile,rog);
  if (calcRogmax) {
    rogmax = DSL->AddMulti(DOUBLE, rogname, "Max");
    //rogmax = DSL->Add(DOUBLE, NULL, "RoGMax");
    if (rogmax == NULL) return 1; 
    DFL->Add(rogFile,rogmax);
  }

  mprintf("    RADGYR: Calculating for atoms in mask %s",Mask1.MaskString());
  if (useMass)
    mprintf(" using mass weighting");
  mprintf(".\n");
  if (!calcRogmax)
    mprintf("            RoG max will not be stored.\n");

  return 0;
}

// Radgyr::setup()
/** Set radius of gyration up for this parmtop. Get masks etc. */
// currentParm is set in Action::Setup
int Radgyr::setup() {

  if ( currentParm->SetupIntegerMask(Mask1, activeReference)) return 1;
  mprintf("\t%s (%i atoms).\n",Mask1.MaskString(),Mask1.Nselected);
  if (Mask1.None()) {
    mprintf("Warning: Radgyr::setup: Mask contains 0 atoms.\n");
    return 1;
  }

  return 0;  
}

// Radgyr::action()
int Radgyr::action() {
  double Rog, max;

  Rog = currentFrame->RADGYR(&Mask1, useMass, &max);

  rog->Add(frameNum, &Rog);
  if (calcRogmax)
    rogmax->Add(frameNum, &max);

  //fprintf(outfile,"%10i %10.4lf %10.4lf\n",frameNum,Rog,max);
  
  return 0;
} 

