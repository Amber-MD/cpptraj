// Radius of Gyration 
#include "Action_Radgyr.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Radgyr::Action_Radgyr() :
  rog_(NULL),
  rogmax_(NULL),
  calcRogmax_(true)
{
  //fprintf(stderr,"Radgyr Con\n");
  useMass_ = false;
} 

// Action_Radgyr::init()
/** Expected call: radgyr <name> <mask1> [out filename] [mass] [nomax]
  */
int Action_Radgyr::init() {
  // Get keywords
  char* rogFile = actionArgs.getKeyString("out",NULL);
  useMass_ = actionArgs.hasKey("mass");
  calcRogmax_ = !actionArgs.hasKey("nomax");

  // Get Masks
  Mask1_.SetMaskString( actionArgs.getNextMask() );

  // Datasets to store radius of gyration and max
  // Also add datasets to data file list
  rog_ = DSL->Add(DataSet::DOUBLE, actionArgs.getNextString(), "RoG");
  if (rog_==NULL) return 1;
  DFL->Add(rogFile, rog_);
  if (calcRogmax_) {
    rogmax_ = DSL->AddSetAspect(DataSet::DOUBLE, rog_->Name(), "Max");
    if (rogmax_ == NULL) return 1; 
    DFL->Add(rogFile, rogmax_);
  }

  mprintf("    RADGYR: Calculating for atoms in mask %s",Mask1_.MaskString());
  if (useMass_)
    mprintf(" using mass weighting");
  mprintf(".\n");
  if (!calcRogmax_)
    mprintf("            RoG max will not be stored.\n");

  return 0;
}

// Action_Radgyr::setup()
/** Set radius of gyration up for this parmtop. Get masks etc. */
// currentParm is set in Action::Setup
int Action_Radgyr::setup() {
  if ( currentParm->SetupIntegerMask(Mask1_)) return 1;
  mprintf("\t%s (%i atoms).\n",Mask1_.MaskString(),Mask1_.Nselected());
  if (Mask1_.None()) {
    mprintf("Warning: Radgyr::setup: Mask contains 0 atoms.\n");
    return 1;
  }

  return 0;  
}

// Action_Radgyr::action()
int Action_Radgyr::action() {
  double max;

  double Rog = currentFrame->RADGYR(&Mask1_, useMass_, &max);

  rog_->Add(frameNum, &Rog);
  if (calcRogmax_)
    rogmax_->Add(frameNum, &max);

  //fprintf(outfile,"%10i %10.4lf %10.4lf\n",frameNum,Rog,max);
  
  return 0;
} 

