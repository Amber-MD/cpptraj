#include <cmath>
#include "Action_Distance.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Distance::Action_Distance() :
  dist_(0),
  mode_(NORMAL),
  plane_(XY),
  useMass_(true)
{}

// Action_Distance::Help()
void Action_Distance::Help() const {
  mprintf("\t[<name>] <mask1> <mask2> [out <filename>] [geom] [noimage] [type noe]\n"
          "\tOptions for 'type noe':\n"
          "\t  %s\n"
          "  Calculate distance between atoms in <mask1> and <mask2>\n", 
          AssociatedData_NOE::HelpText);
}

// Action_Distance::Init()
Action::RetType Action_Distance::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  AssociatedData_NOE noe;
  // Get Keywords
  image_.InitImaging( !(actionArgs.hasKey("noimage")) );
  useMass_ = !(actionArgs.hasKey("geom"));
  std::string modestring = actionArgs.GetStringKey("mode");
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
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
    mprinterr("Error: distance requires 2 masks\n");
    return Action::ERR;
  }
  Mask1_.SetMaskString(mask1);
  Mask2_.SetMaskString(mask2);

  // Dataset to store distances TODO store masks in data set?
  dist_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(actionArgs.GetStringNext(),
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
Action::RetType Action_Distance::Setup(ActionSetup& setup) {
  if (setup.Top().SetupIntegerMask( Mask1_ )) return Action::ERR;
  if (setup.Top().SetupIntegerMask( Mask2_ )) return Action::ERR;
  mprintf("\t%s (%i atoms) to %s (%i atoms)",Mask1_.MaskString(), Mask1_.Nselected(),
          Mask2_.MaskString(),Mask2_.Nselected());
  if (Mask1_.None() || Mask2_.None()) {
    mprintf("\nWarning: One or both masks have no atoms.\n");
    return Action::SKIP;
  }
  // Set up imaging info for this parm
  image_.SetupImaging( setup.CoordInfo().TrajBox().Type() );
  if (image_.ImagingEnabled())
    mprintf(", imaged");
  else
    mprintf(", imaging off");
  mprintf(".\n");

  return Action::OK;  
}

// Action_Distance::DoAction()
Action::RetType Action_Distance::DoAction(int frameNum, ActionFrame& frm) {
  double Dist;
  Matrix_3x3 ucell, recip;
  Vec3 a1, a2;

  if (useMass_) {
    a1 = frm.Frm().VCenterOfMass( Mask1_ );
    a2 = frm.Frm().VCenterOfMass( Mask2_ );
  } else {
    a1 = frm.Frm().VGeometricCenter( Mask1_ );
    a2 = frm.Frm().VGeometricCenter( Mask2_ );
  }

  switch ( image_.ImageType() ) {
    case NONORTHO:
      frm.Frm().BoxCrd().ToRecip(ucell, recip);
      Dist = DIST2_ImageNonOrtho(a1, a2, ucell, recip);
      break;
    case ORTHO:
      Dist = DIST2_ImageOrtho(a1, a2, frm.Frm().BoxCrd());
      break;
    case NOIMAGE:
      Dist = DIST2_NoImage(a1, a2);
      break;
  }
  Dist = sqrt(Dist);

  dist_->Add(frameNum, &Dist);

  return Action::OK;
}
