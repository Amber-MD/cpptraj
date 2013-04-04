// Action_Distance
#include <cmath>
#include "Action_Distance.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Distance::Action_Distance() :
  dist_(0),
  useMass_(true)
{ } 

void Action_Distance::Help() {
  mprintf("\t[<name>] <mask1> <mask2> [out <filename>] [geom] [noimage]\n");
}

// Action_Distance::init()
Action::RetType Action_Distance::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get Keywords
  InitImaging( !(actionArgs.hasKey("noimage")) );
  useMass_ = !(actionArgs.hasKey("geom"));
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  DataSet::scalarType stype = DataSet::UNDEFINED;
  std::string stypename = actionArgs.GetStringKey("type");
  if      ( stypename == "hbond" ) stype = DataSet::HBOND;
  else if (stypename == "noe"    ) stype = DataSet::NOE; // TODO: Grab bound and boundh

  // Get Masks
  std::string mask1 = actionArgs.GetMaskNext();
  std::string mask2 = actionArgs.GetMaskNext();
  if (mask1.empty() || mask2.empty()) {
    mprinterr("Error: distance: Requires 2 masks\n");
    return Action::ERR;
  }
  Mask1_.SetMaskString(mask1);
  Mask2_.SetMaskString(mask2);

  // Dataset to store distances
  dist_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "Dis");
  if (dist_==0) return Action::ERR;
  dist_->SetScalar( DataSet::M_DISTANCE, stype );
  // Add dataset to data file
  if (outfile != 0) outfile->AddSet( dist_ );

  mprintf("    DISTANCE: %s to %s",Mask1_.MaskString(), Mask2_.MaskString());
  if (!UseImage()) 
    mprintf(", non-imaged");
  if (useMass_) 
    mprintf(", center of mass");
  else
    mprintf(", geometric center");
  mprintf(".\n");

  return Action::OK;
}

// Action_Distance::setup()
/** Determine what atoms each mask pertains to for the current parm file.
  * Imaging is checked for in Action::Setup. 
  */
Action::RetType Action_Distance::Setup(Topology* currentParm, Topology** parmAddress) {
  if (currentParm->SetupIntegerMask( Mask1_ )) return Action::ERR;
  if (currentParm->SetupIntegerMask( Mask2_ )) return Action::ERR;
  mprintf("\t%s (%i atoms) to %s (%i atoms)",Mask1_.MaskString(), Mask1_.Nselected(),
          Mask2_.MaskString(),Mask2_.Nselected());
  if (Mask1_.None() || Mask2_.None()) {
    mprintf("\nWarning: distance: One or both masks have no atoms.\n");
    return Action::ERR;
  }
  // Set up imaging info for this parm
  SetupImaging( currentParm->BoxType() );
  if (ImagingEnabled())
    mprintf(", imaged");
  else
    mprintf(", imaging off");
  mprintf(".\n");

  return Action::OK;  
}

// Action_Distance::action()
Action::RetType Action_Distance::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  double Dist;
  Matrix_3x3 ucell, recip;
  Vec3 a1, a2;

  if (useMass_) {
    a1 = currentFrame->VCenterOfMass( Mask1_ );
    a2 = currentFrame->VCenterOfMass( Mask2_ );
  } else {
    a1 = currentFrame->VGeometricCenter( Mask1_ );
    a2 = currentFrame->VGeometricCenter( Mask2_ );
  }

  switch ( ImageType() ) {
    case NONORTHO:
      currentFrame->BoxCrd().ToRecip(ucell, recip);
      Dist = DIST2_ImageNonOrtho(a1, a2, ucell, recip);
      break;
    case ORTHO:
      Dist = DIST2_ImageOrtho(a1, a2, currentFrame->BoxCrd());
      break;
    case NOIMAGE:
      Dist = DIST2_NoImage(a1, a2);
      break;
  }
  Dist = sqrt(Dist);

  dist_->Add(frameNum, &Dist);

  //fprintf(outfile,"%10i %10.4lf\n",frameNum,D);
  
  return Action::OK;
} 
