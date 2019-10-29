#include <cmath>
#include "Action_Distance.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Distance::Action_Distance() :
  a2_(0.0),
  dist_(0),
  mode_(NORMAL),
  useMass_(true)
{}

// Action_Distance::Help()
void Action_Distance::Help() const {
  mprintf("\t[<name>] <mask1> [<mask2>] [point <X> <Y> <Z>]\n"
          "\t[ %s ]\n", DataSetList::RefArgs);
  mprintf("\t[out <filename>] [geom] [noimage] [type noe]\n"
          "\tOptions for 'type noe':\n"
          "\t  %s\n"
          "  Calculate distance between atoms in <mask1> and <mask2>, between\n"
          "  atoms in <mask1> and atoms in <mask2> in specified reference, or\n"
          "  atoms in <mask1> and the specified point.\n",
          AssociatedData_NOE::HelpText);
}

// Action_Distance::Init()
Action::RetType Action_Distance::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  AssociatedData_NOE noe;
  // Get Keywords
  image_.InitImaging( !(actionArgs.hasKey("noimage")) );
  useMass_ = !(actionArgs.hasKey("geom"));
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  MetaData::scalarType stype = MetaData::UNDEFINED;
  std::string stypename = actionArgs.GetStringKey("type");
  if ( stypename == "noe" ) {
    stype = MetaData::NOE;
    if (noe.NOE_Args( actionArgs )) return Action::ERR;
  }
  // Determine mode.
  ReferenceFrame refFrm = init.DSL().GetReferenceFrame( actionArgs );
  if (refFrm.error()) return Action::ERR;
  mode_ = NORMAL;
  if (!refFrm.empty())
    mode_ = REF;
  else if (actionArgs.hasKey("point")) {
    mode_ = POINT;
    a2_[0] = actionArgs.getNextDouble(0.0);
    a2_[1] = actionArgs.getNextDouble(0.0);
    a2_[2] = actionArgs.getNextDouble(0.0);
  }

  // Get Masks
  std::string maskexp = actionArgs.GetMaskNext();
  if (maskexp.empty()) {
    mprinterr("Error: Need at least 1 atom mask.\n");
    return Action::ERR;
  }
  if (Mask1_.SetMaskString(maskexp)) return Action::ERR;
  if (mode_ != POINT) {
    maskexp = actionArgs.GetMaskNext();
    if (maskexp.empty()) {
      mprinterr("Error: Need 2 atom masks.\n");
      return Action::ERR;
    }
    if (Mask2_.SetMaskString(maskexp)) return Action::ERR;
  }

  // Set up reference and get reference point
  if (mode_ == REF) {
    if (refFrm.Parm().SetupIntegerMask( Mask2_, refFrm.Coord() ))
      return Action::ERR;
    if (useMass_)
      a2_ = refFrm.Coord().VCenterOfMass( Mask2_ );
    else
      a2_ = refFrm.Coord().VGeometricCenter( Mask2_ );
  }

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

  mprintf("    DISTANCE:");
  if (mode_ == NORMAL)
    mprintf(" %s to %s", Mask1_.MaskString(), Mask2_.MaskString());
  else if (mode_ == REF)
    mprintf(" %s to %s (%i atoms) in %s", Mask1_.MaskString(),
            Mask2_.MaskString(), Mask2_.Nselected(), refFrm.refName());
  else if (mode_ == POINT)
    mprintf(" %s to point {%g %g %g}", Mask1_.MaskString(), a2_[0], a2_[1], a2_[2]);
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
  if (mode_ == NORMAL) {
    if (setup.Top().SetupIntegerMask( Mask2_ )) return Action::ERR;
    mprintf("\t%s (%i atoms) to %s (%i atoms)",
            Mask1_.MaskString(), Mask1_.Nselected(),
            Mask2_.MaskString(),Mask2_.Nselected());
    if (Mask1_.None() || Mask2_.None()) {
      mprintf("\nWarning: One or both masks have no atoms.\n");
      return Action::SKIP;
    }
  } else {
    mprintf("\t%s (%i atoms)", Mask1_.MaskString(), Mask1_.Nselected());
    if (Mask1_.None()) {
      mprintf("\nWarning: Mask has no atoms.\n");
      return Action::SKIP;
    }
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
  Vec3 a1;

  if ( mode_ == NORMAL ) {
    if (useMass_) {
      a1  = frm.Frm().VCenterOfMass( Mask1_ );
      a2_ = frm.Frm().VCenterOfMass( Mask2_ );
    } else {
      a1  = frm.Frm().VGeometricCenter( Mask1_ );
      a2_ = frm.Frm().VGeometricCenter( Mask2_ );
    }
  } else { // REF, POINT
    if (useMass_)
      a1 = frm.Frm().VCenterOfMass( Mask1_ );
    else
      a1 = frm.Frm().VGeometricCenter( Mask1_ );
  }

  switch ( image_.ImageType() ) {
    case NONORTHO:
      frm.Frm().BoxCrd().ToRecip(ucell, recip);
      Dist = DIST2_ImageNonOrtho(a1, a2_, ucell, recip);
      break;
    case ORTHO:
      Dist = DIST2_ImageOrtho(a1, a2_, frm.Frm().BoxCrd());
      break;
    case NOIMAGE:
      Dist = DIST2_NoImage(a1, a2_);
      break;
  }
  Dist = sqrt(Dist);

  dist_->Add(frameNum, &Dist);

  return Action::OK;
}
