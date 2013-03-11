// Action_Center 
#include "Action_Center.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Center::Action_Center() :
  origin_(false),
  useMass_(false)
{ } 

void Action_Center::Help() {
  mprintf("\t<mask> [origin] [mass]\n\tCenter coordinates in <mask>.\n");
}

// Action_Center::init()
Action::RetType Action_Center::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  origin_ = actionArgs.hasKey("origin");
  useMass_ = actionArgs.hasKey("mass");

  // Get Masks
  Mask_.SetMaskString( actionArgs.GetMaskNext() );

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

  return Action::OK;
}

// Action_Center::setup()
/** Set angle up for this parmtop. Get masks etc. */
// currentParm is set in Action::Setup
Action::RetType Action_Center::Setup(Topology* currentParm, Topology** parmAddress) {

  if ( currentParm->SetupIntegerMask(Mask_) ) return Action::ERR;
  Mask_.MaskInfo();
  if (Mask_.None()) {
    mprintf("Warning: center:: Mask contains 0 atoms.\n");
    return Action::ERR;
  }

  if (!origin_ && currentParm->BoxType()==Box::NOBOX) {
    mprintf("Warning: center: Box center specified but no box information.\n");
    //fprintf(stdout,"                            Centering on origin.\n");
    return Action::ERR;
  }

  return Action::OK;  
}

// Action_Center::action()
/** Center coordinates in frame to coord origin or box origin (corner).
  */
Action::RetType Action_Center::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {

  currentFrame->Center(Mask_, origin_, useMass_);

  return Action::OK;
} 
