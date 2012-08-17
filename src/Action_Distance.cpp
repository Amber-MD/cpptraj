// Action_Distance
#include <cmath>
#include "Action_Distance.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Distance::Action_Distance() :
  dist_(NULL)
{
  //fprintf(stderr,"Action_Distance Con\n");
  useImage_=true;
  useMass_=true;
} 

// Action_Distance::init()
/** Expected call: distance <name> <mask1> <mask2> [out filename] [geom] [noimage]
  */
int Action_Distance::init( ) {
  // Get Keywords
  useImage_ = !(actionArgs.hasKey("noimage"));
  useMass_ = !(actionArgs.hasKey("geom"));
  ArgList::ConstArg distanceFile = actionArgs.getKeyString("out");
  DataSet::scalarType stype = DataSet::UNDEFINED;
  std::string stypename = actionArgs.GetStringKey("type");
  if      ( stypename == "hbond" ) stype = DataSet::HBOND;
  else if (stypename == "noe"    ) stype = DataSet::NOE; // TODO: Grab bound and boundh

  // Get Masks
  ArgList::ConstArg mask1 = actionArgs.getNextMask();
  ArgList::ConstArg mask2 = actionArgs.getNextMask();
  //fprintf(stdout,"    Mask 1: %s\n",mask1);
  //fprintf(stdout,"    Mask 2: %s\n",mask2);
  if (mask1==NULL || mask2==NULL) {
    mprinterr("Error: distance: Requires 2 masks\n");
    return 1;
  }
  Mask1_.SetMaskString(mask1);
  Mask2_.SetMaskString(mask2);

  // Dataset to store distances
  dist_ = DSL->Add(DataSet::DOUBLE, actionArgs.getNextString(),"Dis");
  if (dist_==NULL) return 1;
  dist_->SetScalar( DataSet::M_DISTANCE, stype );
  // Add dataset to data file list
  DFL->Add(distanceFile, dist_);

  mprintf("    DISTANCE: %s to %s",Mask1_.MaskString(), Mask2_.MaskString());
  if (!useImage_) 
    mprintf(", non-imaged");
  if (useMass_) 
    mprintf(", center of mass");
  else
    mprintf(", geometric center");
  mprintf(".\n");

  return 0;
}

// Action_Distance::setup()
/** Determine what atoms each mask pertains to for the current parm file.
  * Imaging is checked for in Action::Setup. 
  */
int Action_Distance::setup() {
  if (currentParm->SetupIntegerMask( Mask1_ )) return 1;
  if (currentParm->SetupIntegerMask( Mask2_ )) return 1;

  // Print mask and imaging info for this parm
  mprintf("\t%s (%i atoms) to %s (%i atoms)",Mask1_.MaskString(), Mask1_.Nselected(),
          Mask2_.MaskString(),Mask2_.Nselected());
  if (imageType_ != Frame::NOIMAGE)
    mprintf(", imaged");
  else
    mprintf(", imaging off");
  mprintf(".\n");

  if (Mask1_.None() || Mask2_.None()) {
    mprintf("Warning: distance: One or both masks have no atoms.\n");
    return 1;
  }
       
  return 0;  
}

// Action_Distance::action()
int Action_Distance::action() {
  double ucell[9], recip[9];

  if (imageType_==Frame::NONORTHO) currentFrame->BoxToRecip(ucell,recip);
  double Dist = currentFrame->DIST2(Mask1_, Mask2_, useMass_, imageType_, ucell, recip);
  Dist = sqrt(Dist);

  dist_->Add(frameNum, &Dist);

  //fprintf(outfile,"%10i %10.4lf\n",frameNum,D);
  
  return 0;
} 

