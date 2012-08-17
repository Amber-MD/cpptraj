// Action_Center 
#include "Action_Center.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Center::Action_Center() :
  origin_(false)
{
  //fprintf(stderr,"Center Con\n");
  useMass_ = false;
} 

// Action_Center::init()
/** Expected call: center <mask> [origin] [mass] 
  */
int Action_Center::init() {
  // Get keywords
  origin_ = actionArgs.hasKey("origin");
  useMass_ = actionArgs.hasKey("mass");

  // Get Masks
  Mask_.SetMaskString( actionArgs.getNextMask() );

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

// Action_Center::setup()
/** Set angle up for this parmtop. Get masks etc. */
// currentParm is set in Action::Setup
int Action_Center::setup() {

  if ( currentParm->SetupIntegerMask(Mask_) ) return 1;
  if (Mask_.None()) {
    mprintf("Warning: center:: Mask contains 0 atoms.\n");
    return 1;
  }
  mprintf("\t%s (%i atoms)\n",Mask_.MaskString(), Mask_.Nselected());

  if (!origin_ && currentParm->BoxType()==Box::NOBOX) {
    mprintf("Warning: center: Box center specified but no box information.\n");
    //fprintf(stdout,"                            Centering on origin.\n");
    return 1;
  }

  return 0;  
}

// Action_Center::action()
/** Center coordinates in frame to coord origin or box origin (corner).
  */
int Action_Center::action() {

  currentFrame->Center(Mask_, origin_, useMass_);

  return 0;
} 
