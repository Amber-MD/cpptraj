#include "Action_Unwrap.h"
#include "CpptrajStdio.h"
#include "ImageRoutines.h"

// CONSTRUCTOR
Action_Unwrap::Action_Unwrap() :
  imageMode_(Image::BYATOM),
  RefParm_(0),
  orthogonal_(false),
  center_(false)
{ }

void Action_Unwrap::Help() const {
  mprintf("\t[center] [{bymol | byres | byatom}]\n"
          "\t[ %s ] [<mask>]\n", DataSetList::RefArgs);
  mprintf("  Reverse of 'image'; unwrap coordinates in <mask> according\n"
          "  to a reference structure.\n");
}

// Action_Unwrap::Init()
Action::RetType Action_Unwrap::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get Keywords
  center_ = actionArgs.hasKey("center");
  if (actionArgs.hasKey("bymol"))
    imageMode_ = Image::BYMOL;
  else if (actionArgs.hasKey("byres"))
    imageMode_ = Image::BYRES;
  else if (actionArgs.hasKey("byatom")) {
    imageMode_ = Image::BYATOM;
    // Unwrapping to center by atom makes no sense
    if (center_) center_ = false;
  } else
    imageMode_ = Image::BYATOM;
  // Get reference
  ReferenceFrame REF = init.DSL().GetReferenceFrame( actionArgs );
  if (REF.error()) return Action::ERR;
  if (!REF.empty()) {
    RefFrame_ = REF.Coord();
    // Get reference parm for frame
    RefParm_ = REF.ParmPtr();
  }

  // Get mask string
  maskExpression_ = actionArgs.GetMaskNext();

  mprintf("    UNWRAP: By %s", Image::ModeString(imageMode_));
  if (!maskExpression_.empty())
    mprintf(" using mask '%s'", maskExpression_.c_str());
  else
    mprintf(" using all atoms");
  if (imageMode_ != Image::BYATOM) {
    if (center_)
      mprintf(" based on center of mass.");
    else
      mprintf(" based on first atom position.");
  }
  mprintf("\n");
  if ( !REF.empty())
    mprintf("\tReference is %s", REF.refName());
  else
    mprintf("\tReference is first frame.");
  mprintf("\n");

  return Action::OK;
}

// Action_Unwrap::Setup()
Action::RetType Action_Unwrap::Setup(ActionSetup& setup) {
  // Ensure same number of atoms in current parm and ref parm
  if ( RefParm_!=0 ) {
    if ( setup.Top().Natom() != RefParm_->Natom() ) {
      mprinterr("Error: unwrap: # atoms in reference parm %s is not\n", RefParm_->c_str());
      mprinterr("Error:         equal to # atoms in parm %s\n", setup.Top().c_str());
      return Action::ERR;
    }
  }
  // Check box type
  if (setup.CoordInfo().TrajBox().Type()==Box::NOBOX) {
    mprintf("Error: unwrap: Parm %s does not contain box information.\n",
            setup.Top().c_str());
    return Action::ERR;
  }
  orthogonal_ = false;
  if (setup.CoordInfo().TrajBox().Type()==Box::ORTHO)
    orthogonal_ = true;

  // Setup atom pairs to be unwrapped.
  imageList_ = Image::CreatePairList(setup.Top(), imageMode_, maskExpression_);
  if (imageList_.empty()) {
    mprintf("Warning: Mask selects no atoms for topology '%s'.\n", setup.Top().c_str());
    return Action::SKIP;
  }
  mprintf("\tNumber of %ss to be unwrapped is %zu\n",
          Image::ModeString(imageMode_), imageList_.size()/2);

  // Use current parm as reference if not already set
  if (RefParm_ == 0)
    RefParm_ = setup.TopAddress();
  return Action::OK;
}

// Action_Unwrap::DoAction()
Action::RetType Action_Unwrap::DoAction(int frameNum, ActionFrame& frm) {
  Matrix_3x3 ucell, recip;
  // Set reference structure if not already set
  if (RefFrame_.empty()) {
    RefFrame_ = frm.Frm();
    return Action::OK;
  }
 
  if (orthogonal_)
    Image::UnwrapOrtho( frm.ModifyFrm(), RefFrame_, imageList_, center_, true );
  else {
    frm.Frm().BoxCrd().ToRecip( ucell, recip );
    Image::UnwrapNonortho( frm.ModifyFrm(), RefFrame_, imageList_, ucell, recip, center_, true );
  }

  return Action::MODIFY_COORDS;
}
