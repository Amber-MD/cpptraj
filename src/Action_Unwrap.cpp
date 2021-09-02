#include "Action_Unwrap.h"
#include "CpptrajStdio.h"
#include "ImageRoutines.h"
#include "Image_List.h"

// CONSTRUCTOR
Action_Unwrap::Action_Unwrap() :
  imageList_(0),
  imageMode_(Image::BYATOM),
  RefParm_(0),
  center_(false)
{ }

/** DESTRUCTOR */
Action_Unwrap::~Action_Unwrap() {
  if (imageList_ != 0) delete imageList_;
}

void Action_Unwrap::Help() const {
  mprintf("\t[center] [{bymol | byres | byatom}]\n"
          "\t[ %s ] [<mask>]\n", DataSetList::RefArgs);
  mprintf("  Reverse of 'image'; unwrap coordinates in <mask> according\n"
          "  to a reference structure.\n");
}

// Action_Unwrap::Init()
Action::RetType Action_Unwrap::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  if (init.TrajComm().Size() > 1) {
    mprinterr("Error: 'unwrap' action does not work with > 1 process (%i processes currently).\n",
              init.TrajComm().Size());
    return Action::ERR;
  }
# endif
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
  if (!setup.CoordInfo().TrajBox().HasBox()) {
    mprintf("Error: unwrap: Parm %s does not contain box information.\n",
            setup.Top().c_str());
    return Action::ERR;
  }

  // Setup atom pairs to be unwrapped. Always use CoM TODO why?
  if (imageList_ != 0) delete imageList_;
  imageList_ = Image::CreateImageList(setup.Top(), imageMode_, maskExpression_,
                                      true, center_);
  if (imageList_ == 0) {
    mprinterr("Internal Error: Could not allocate unwrap list.\n");
    return Action::ERR;
  }
  if (imageList_->nEntities() < 1) {
    mprintf("Warning: Mask selects no atoms for topology '%s'.\n", setup.Top().c_str());
    return Action::SKIP;
  }
  mprintf("\tNumber of %ss to be unwrapped is %u\n",
          Image::ModeString(imageMode_), imageList_->nEntities());
  // Get entities that need to be updated in reference
  allEntities_ = imageList_->AllEntities();

  // Use current parm as reference if not already set
  if (RefParm_ == 0)
    RefParm_ = setup.TopAddress();
  return Action::OK;
}

// Action_Unwrap::DoAction()
Action::RetType Action_Unwrap::DoAction(int frameNum, ActionFrame& frm) {
  if (RefFrame_.empty()) {
    // Set reference structure if not already set
    RefFrame_ = frm.Frm();
    return Action::OK;
  }
 
  if (frm.Frm().BoxCrd().Is_X_Aligned_Ortho())
    Image::UnwrapOrtho( frm.ModifyFrm(), RefFrame_, *imageList_, allEntities_ );
  else {
    Image::UnwrapNonortho( frm.ModifyFrm(), RefFrame_, *imageList_, allEntities_, frm.Frm().BoxCrd().UnitCell(), frm.Frm().BoxCrd().FracCell() );
  }

  return Action::MODIFY_COORDS;
}
