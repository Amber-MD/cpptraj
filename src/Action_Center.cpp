// Center 
#include "Action_Center.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Center::Center() {
  //fprintf(stderr,"Center Con\n");
  box[0]=0.0;
  box[1]=0.0;
  box[2]=0.0;
  origin=false;
  useMass=false;
} 

// Center::init()
/** Expected call: center <mask> [origin] [mass] 
  */
int Center::init() {
  char *mask1;

  // Get keywords
  origin = actionArgs.hasKey("origin");
  useMass = actionArgs.hasKey("mass");

  // Get Masks
  mask1 = actionArgs.getNextMask();
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
  mprintf(" using atoms in mask %s\n",Mask1.MaskString());

  return 0;
}

// Center::setup()
/** Set angle up for this parmtop. Get masks etc.
  */
// currentParm is set in Action::Setup
int Center::setup() {

  if ( currentParm->SetupIntegerMask(Mask1, activeReference) ) return 1;
  if (Mask1.None()) {
    mprintf("Warning: Center::setup: Mask contains 0 atoms.\n");
    return 1;
  }
  mprintf("\t%s (%i atoms)\n",Mask1.MaskString(),Mask1.Nselected);

  if (!origin && currentParm->boxType==NOBOX) {
    mprintf("Warning: Center::setup: Box center specified but no box information.\n");
    //fprintf(stdout,"                            Centering on origin.\n");
    return 1;
  }

  return 0;  
}

// Center::action()
/** Center coordinates in frame to coord origin or box origin (corner).
  */
int Center::action() {

  // Set up box
  if (!origin) {
    box[0] = currentFrame->box[0] / 2.0;
    box[1] = currentFrame->box[1] / 2.0;
    box[2] = currentFrame->box[2] / 2.0;
  }
  
  currentFrame->Center(&Mask1, box, useMass);

  return 0;
} 

