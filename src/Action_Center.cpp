// Center 
#include "Action_Center.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Center::Center() :
  origin_(false)
{
  //fprintf(stderr,"Center Con\n");
  useMass_ = false;
} 

// Center::init()
/** Expected call: center <mask> [origin] [mass] 
  */
int Center::init() {
  // Get keywords
  origin_ = actionArgs.hasKey("origin");
  useMass_ = actionArgs.hasKey("mass");

  // Get Masks
  char *mask1 = actionArgs.getNextMask();
  Mask_.SetMaskString(mask1);

  mprintf("    CENTER: To");
  if (origin_)
    mprintf(" origin");
  else
    mprintf(" box center");
  mprintf(" via center of");
  if (useMass_)
    mprintf(" mass");
  else
    mprintf(" geometry");
  mprintf(" using atoms in mask %s\n",Mask_.MaskString());

  return 0;
}

// Center::setup()
/** Set angle up for this parmtop. Get masks etc.
  */
// currentParm is set in Action::Setup
int Center::setup() {

  if ( currentParm->SetupIntegerMask(Mask_) ) return 1;
  if (Mask_.None()) {
    mprintf("Warning: Center::setup: Mask contains 0 atoms.\n");
    return 1;
  }
  mprintf("\t%s (%i atoms)\n",Mask_.MaskString(), Mask_.Nselected());

  if (!origin_ && currentParm->BoxType()==Box::NOBOX) {
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

  currentFrame->Center(Mask_, origin_, useMass_);

  return 0;
} 

