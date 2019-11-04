// Action_Center 
#include "Action_Center.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Center::Action_Center() :
  centerMode_(BOXCTR),
  useMass_(false)
{ } 

void Action_Center::Help() const {
  mprintf("\t<mask> [mass] {origin | point <X> <Y> <Z>|\n"
          "\t %s [<refmask>]}\n"
          "  Center coordinates in <mask> to specified point. If 'origin' is specified\n"
          "  center on coordinate origin (0,0,0). If 'point' is specified center on\n"
          "  point specified by XYZ. If a reference is specified center on selected\n"
          "  atoms in the given reference. Otherwise center on box center.\n",
           DataSetList::RefArgs);
}

// Action_Center::Init()
Action::RetType Action_Center::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get keywords
  useMass_ = actionArgs.hasKey("mass");
  ReferenceFrame refFrm = init.DSL().GetReferenceFrame( actionArgs );
  if (refFrm.error()) return Action::ERR;
  // Determine center mode.
  if (!refFrm.empty())
    centerMode_ = REF;
  else if (actionArgs.hasKey("origin"))
    centerMode_ = ORIGIN;
  else if (actionArgs.hasKey("point")) {
    centerMode_ = POINT;
    refCenter_[0] = actionArgs.getNextDouble(0.0);
    refCenter_[1] = actionArgs.getNextDouble(0.0);
    refCenter_[2] = actionArgs.getNextDouble(0.0);
  } else
    centerMode_ = BOXCTR;

  // Get Masks
  if (Mask_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;
  // Get reference mask if reference specified.
  AtomMask refMask;
  if (centerMode_ == REF) {
    std::string rMaskExpr = actionArgs.GetMaskNext();
    if (rMaskExpr.empty())
      rMaskExpr = Mask_.MaskExpression();
    if (refMask.SetMaskString( rMaskExpr )) return Action::ERR;
    if (refFrm.Parm().SetupIntegerMask( refMask, refFrm.Coord() ))
      return Action::ERR;
    // Get center of mask in reference
    if (useMass_)
      refCenter_ = refFrm.Coord().VCenterOfMass( refMask );
    else
      refCenter_ = refFrm.Coord().VGeometricCenter( refMask );
  }

  mprintf("    CENTER: Centering coordinates using");
  if (useMass_)
    mprintf(" center of mass");
  else
    mprintf(" geometric center");
  mprintf(" of atoms in mask (%s) to\n", Mask_.MaskString());
  switch (centerMode_) {
    case ORIGIN: mprintf("\tcoordinate origin.\n"); break;
    case BOXCTR: mprintf("\tbox center.\n"); break;
    case REF:
      mprintf("\tcenter of mask (%s) in reference '%s'.\n", refMask.MaskString(),
              refFrm.refName());
      break;
    case POINT:  mprintf("\tpoint (%g, %g, %g).\n",
                         refCenter_[0], refCenter_[1], refCenter_[2]);
      break;
  }
  return Action::OK;
}

// Action_Center::Setup()
Action::RetType Action_Center::Setup(ActionSetup& setup) {

  if ( setup.Top().SetupIntegerMask(Mask_) ) return Action::ERR;
  Mask_.MaskInfo();
  if (Mask_.None()) {
    mprintf("Warning: Mask contains 0 atoms.\n");
    return Action::SKIP;
  }

  if (centerMode_ == BOXCTR && setup.CoordInfo().TrajBox().Type()==Box::NOBOX) {
    mprintf("Warning: Box center specified but no box information.\n");
    return Action::SKIP;
  }

  return Action::OK;  
}

// Action_Center::DoAction()
/** Center coordinates in frame according to specified mode. */
Action::RetType Action_Center::DoAction(int frameNum, ActionFrame& frm) {
  Vec3 center;
  if (useMass_)
    center = frm.Frm().VCenterOfMass(Mask_);
  else
    center = frm.Frm().VGeometricCenter(Mask_);
  //mprinterr("  FRAME CENTER: %lf %lf %lf\n",center[0],center[1],center[2]); //DEBUG
  switch (centerMode_) {
    case ORIGIN: // Shift to coordinate origin (0,0,0)
      center.Neg(); break;
    case BOXCTR: // Shift to box center
      center = frm.Frm().BoxCrd().Center() - center; break;
    case POINT:  // Shift to specified point
    case REF:    // Shift to reference point
      center = refCenter_ - center; break;
  }
  frm.ModifyFrm().Translate(center);

  return Action::MODIFY_COORDS;
}
