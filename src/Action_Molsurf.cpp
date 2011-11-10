// Molsurf 
#include "Action_Molsurf.h"
#include "CpptrajStdio.h"
#include "molsurf.h"

// CONSTRUCTOR
Molsurf::Molsurf() {
  //fprintf(stderr,"Angle Con\n");
  sasa=NULL;
} 

// DESTRUCTOR
Molsurf::~Molsurf() { }

// Molsurf::init()
/** Expected call: molsurf [<name>] [<mask1>] [out filename] 
  * Dataset name will be the last arg checked for. Check order is:
  *    1) Keywords
  *    2) Masks
  *    3) Dataset name
  */
int Molsurf::init() {
  char *mask1;
  char *molsurfFile;

  // Get keywords
  molsurfFile = actionArgs.getKeyString("out",NULL);

  // Get Masks
  mask1 = actionArgs.getNextMask();
  Mask1.SetMaskString(mask1);

  // Dataset to store angles
  sasa = DSL->Add(DOUBLE, actionArgs.getNextString(),"MSURF");
  if (sasa==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(molsurfFile,sasa);

  mprintf("    MOLSURF: [%s]\n",Mask1.maskString);

  return 0;
}

// Molsurf::setup()
/** Set mask up for this parmtop.
  * currentParm is set in Action::Setup
  */
int Molsurf::setup() {

  if ( Mask1.SetupMask(currentParm,activeReference,debug) ) return 1;
  if (Mask1.None()) {
    mprintf("    Error: Molsurf::setup: Mask contains 0 atoms.\n");
    return 1;
  }

  return 0;  
}

// Molsurf::action()
int Molsurf::action() {
  double molsurf_sasa;

  molsurf_sasa = molsurf( 1.4, currentFrame->X, currentFrame->natom, 
                          currentParm->names, currentParm->resnames,
                          currentParm->resnums, currentParm->charge,
                          currentParm->GB_radii() );

  sasa->Add(frameNum, &molsurf_sasa);

  //fprintf(outfile,"%10i %10.4lf\n",frameNum,molsurf_sasa);
  
  return 0;
} 


