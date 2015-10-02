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
  mprintf("\t[<name>] <mask1> <mask2> [out <filename>] [geom] [noimage] [type noe]\n"
          "\tOptions for 'type noe':\n"
          "\t  %s\n"
          "  Calculate distance between atoms in <mask1> and <mask2>\n", 
          AssociatedData_NOE::HelpText);
}

// Action_Distance::Init()
Action::RetType Action_Distance::Init(ArgList& actionArgs, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  AssociatedData_NOE noe;
  // Get Keywords
  image_.InitImaging( !(actionArgs.hasKey("noimage")) );
  useMass_ = !(actionArgs.hasKey("geom"));
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  MetaData::scalarType stype = MetaData::UNDEFINED;
  std::string stypename = actionArgs.GetStringKey("type");
  if ( stypename == "noe" ) {
    stype = MetaData::NOE;
    if (noe.NOE_Args( actionArgs )) return Action::ERR;
  }
  // Get Masks
  std::string mask1 = actionArgs.GetMaskNext();
  std::string mask2 = actionArgs.GetMaskNext();
  if (mask1.empty() || mask2.empty()) {
    mprinterr("Error: distance: Requires 2 masks\n");
    return Action::ERR;
  }
  Mask1_.SetMaskString(mask1);
  Mask2_.SetMaskString(mask2);

  // Dataset to store distances TODO store masks in data set?
  dist_ = DSL->AddSet(DataSet::DOUBLE, MetaData(actionArgs.GetStringNext(),
                                                MetaData::M_DISTANCE, stype), "Dis");
  if (dist_==0) return Action::ERR;
  if ( stype == MetaData::NOE ) {
    dist_->AssociateData( &noe );
    dist_->SetLegend(Mask1_.MaskExpression() + " and " + Mask2_.MaskExpression());
  }
  // Add dataset to data file
  if (outfile != 0) outfile->AddDataSet( dist_ );

  mprintf("    DISTANCE: %s to %s",Mask1_.MaskString(), Mask2_.MaskString());
  if (!image_.UseImage()) 
    mprintf(", non-imaged");
  if (useMass_) 
    mprintf(", center of mass");
  else
    mprintf(", geometric center");
  mprintf(".\n");

  return Action::OK;
}

// Action_Distance::Setup()
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
  image_.SetupImaging( currentParm->BoxType() );
  if (image_.ImagingEnabled())
    mprintf(", imaged");
  else
    mprintf(", imaging off");
  mprintf(".\n");

  return Action::OK;  
}

// Action_Distance::DoAction()
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

  switch ( image_.ImageType() ) {
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
