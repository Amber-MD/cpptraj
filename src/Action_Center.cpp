// Center 
#include "Action_Center.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Center::Center() {
  //fprintf(stderr,"Center Con\n");
//  mass=NULL;
  box[0]=0.0;
  box[1]=0.0;
  box[2]=0.0;
  origin=false;
  useMass=false;
} 

// DESTRUCTOR
Center::~Center() { }

/*
 * Center::init()
 * Expected call: center <mask> [origin] [mass] 
 * Check order is:
 *    1) Keywords
 *    2) Masks
 */
int Center::init() {
  char *mask1;

  // Get keywords
  origin = A->hasKey("origin");
  useMass = A->hasKey("mass");

  // Get Masks
  mask1 = A->getNextMask();
  Mask1.SetMaskString(mask1);

  mprintf("    CENTER: To");
  if (origin)
    mprintf(" origin");
  else
    mprintf(" box center");
  mprintf(" via center of");
  if (useMass)
    mprintf(" mass");
  else
    mprintf(" geometry");
  mprintf(" using atoms in mask %s\n",Mask1.maskString);

  return 0;
}

/*
 * Center::setup()
 * Set angle up for this parmtop. Get masks etc.
 * P is set in Action::Setup
 */
int Center::setup() {

  if ( Mask1.SetupMask(P,debug) ) return 1;
  if (Mask1.None()) {
    mprintf("    Error: Center::setup: Mask contains 0 atoms.\n");
    return 1;
  }

  if (useMass && P->mass==NULL) {
    mprintf("    Warning: Center::setup: usemass: Parm %s contains no mass info.\n",
            P->parmName);
    mprintf("             Geometric center will be used instead.\n");
    useMass=false;
  }
//    mass = P->mass;
// else
//    mass = NULL;

  if (!origin && P->BoxType==0) {
    mprintf("    Error: Center::setup: Box center specified but no box information.\n");
    //fprintf(stdout,"                            Centering on origin.\n");
    return 1;
  }

  return 0;  
}

/*
 * Center::action()
 */
int Center::action() {

  // Set up box
  if (!origin) {
    box[0] = F->box[0] / 2.0;
    box[1] = F->box[1] / 2.0;
    box[2] = F->box[2] / 2.0;
  }
  
  F->Center(&Mask1, box, useMass);

  return 0;
} 


