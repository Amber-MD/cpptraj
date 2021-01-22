// Action_Image 
#include "Action_Image.h"
#include "CpptrajStdio.h"
#include "ImageRoutines.h"
#include "Image_List.h"

// CONSTRUCTOR
Action_Image::Action_Image() :
  imageMode_(Image::BYMOL),
  ComMask_(0),
  offset_(0.0),
  origin_(false),
  center_(false),
  useMass_(true),
  truncoct_(false),
  triclinic_(OFF),
  debug_(0),
  imageList_(0)
{ }

/** DESTRUCTOR */
Action_Image::~Action_Image() {
  if (imageList_ != 0) delete imageList_;
  if (ComMask_!=0) delete ComMask_;
}

void Action_Image::Help() const {
  mprintf("\t[origin] [center] [triclinic | familiar [com <commask>]] <mask>\n" 
          "\t[ bymol | byres | byatom ] [xoffset <x>] [yoffset <y>] [zoffset <z>]\n"
          "  Image atoms in <mask> into the primary unit cell.\n"
          "    origin: center at 0.0, 0.0, 0.0, otherwise center at box center.\n"
          "    center: Use center of mass for imaging, otherwise use first atom.\n"
          "    triclinic: Force imaging with triclinic code.\n"
          "    familiar: Image with triclinic code and shape into familiar trunc. oct. shape.\n"
          "    com <commask>: If familiar, center based on COM of atoms in mask, otherwise use\n"
          "                   origin/box.\n"
          "    <mask>: Only image atoms in <mask>. If no mask given all atoms are imaged.\n");
}

// Action_Image::Init()
Action::RetType Action_Image::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  // Get keywords
  origin_ = actionArgs.hasKey("origin");
  center_ = actionArgs.hasKey("center");
  if (actionArgs.hasKey("familiar")) triclinic_ = FAMILIAR;
  if (actionArgs.hasKey("triclinic")) triclinic_ = FORCE;
  if (actionArgs.hasKey("bymol"))
    imageMode_ = Image::BYMOL;
  else if (actionArgs.hasKey("byres"))
    imageMode_ = Image::BYRES;
  else if (actionArgs.hasKey("byatom")) {
    imageMode_ = Image::BYATOM;
    // Imaging to center by atom makes no sense
    if (center_) center_ = false;
  } else
    imageMode_ = Image::BYMOL;
  offset_[0] = actionArgs.getKeyDouble("xoffset", 0.0);
  offset_[1] = actionArgs.getKeyDouble("yoffset", 0.0);
  offset_[2] = actionArgs.getKeyDouble("zoffset", 0.0);
  // Get Masks
  if (triclinic_ == FAMILIAR) {
    std::string maskexpr = actionArgs.GetStringKey("com");
    if (!maskexpr.empty()) {
      ComMask_ = new AtomMask();
      if (ComMask_->SetMaskString(maskexpr)) return Action::ERR;
    }
  }
  maskExpression_ = actionArgs.GetMaskNext();
  
  mprintf("    IMAGE: By %s to", Image::ModeString(imageMode_));
  if (origin_)
    mprintf(" origin");
  else
    mprintf(" box center");
  if (imageMode_ != Image::BYATOM) {
    if (center_)
      mprintf(" based on center of mass");
    else
      mprintf(" based on first atom position");
  }
  if (!maskExpression_.empty())
    mprintf(" using atoms in mask %s\n", maskExpression_.c_str());
  else
    mprintf(" using all atoms\n");
  if (triclinic_ == FORCE)
    mprintf( "           Triclinic On.\n");
  else if (triclinic_ == FAMILIAR) {
    mprintf( "           Triclinic On, familiar shape");
    if (ComMask_!=0) 
      mprintf( " centering on atoms in mask %s", ComMask_->MaskString());
    mprintf(".\n");
  }
  if (!offset_.IsZero())
    mprintf("\tOffsetting unit cells by factors X=%g, Y=%g, Z=%g\n",
            offset_[0], offset_[1], offset_[2]);

  return Action::OK;
}

// Action_Image::Setup()
/** Set Imaging up for this parmtop. Get masks etc.
  * currentParm is set in Action::Setup
  */
Action::RetType Action_Image::Setup(ActionSetup& setup) {
  // Check box type
  if (!setup.CoordInfo().TrajBox().HasBox()) {
    mprintf("Warning: Topology %s does not contain box information.\n",
            setup.Top().c_str());
    return Action::SKIP;
  }
  // Setup atom pairs to be unwrapped.
  if (imageList_ != 0) delete imageList_;
  imageList_ = Image::CreateImageList(setup.Top(), imageMode_, maskExpression_,
                                      useMass_, center_);
  if (imageList_ == 0) {
    mprinterr("Internal Error: Could not allocate image list.\n");
    return Action::ERR;
  }
  if (imageList_->nEntities() < 1) {
    mprintf("Warning: Nothing selected to image for topology '%s'.\n", setup.Top().c_str());
    return Action::SKIP;
  }
  mprintf("\tNumber of %ss to be imaged is %u\n",
          Image::ModeString(imageMode_), imageList_->nEntities());
  // DEBUG: Print all entities 
  if (debug_>0)
    imageList_->PrintEntities();
  // Setup for truncated octahedron
  if (triclinic_ == FAMILIAR) {
    if (ComMask_!=0) {
      if ( setup.Top().SetupIntegerMask( *ComMask_ ) ) return Action::ERR;
      if (ComMask_->None()) {
        mprintf("Warning: Mask for 'familiar com' contains no atoms.\n");
        return Action::SKIP;
      }
      mprintf("\tcom: mask [%s] contains %i atoms.\n",ComMask_->MaskString(),ComMask_->Nselected());
    }
  }

  // Truncoct flag
  truncoct_ = (triclinic_ == FAMILIAR);

  return Action::OK;
}

// Action_Image::DoAction()
Action::RetType Action_Image::DoAction(int frameNum, ActionFrame& frm) {
  if (frm.Frm().BoxCrd().Is_X_Aligned_Ortho() && triclinic_ == OFF) {
    Vec3 bp, bm;
    if (Image::SetupOrtho(frm.Frm().BoxCrd(), bp, bm, origin_)) {
      mprintf("Warning: Frame %i imaging failed, box lengths are zero.\n",frameNum+1);
      // TODO: Return OK for now so next frame is tried; eventually indicate SKIP?
      return Action::OK;
    }
    Image::Ortho(frm.ModifyFrm(), bp, bm, offset_, *imageList_);
  } else {
    Vec3 fcom;
    if (truncoct_)
      fcom = Image::SetupTruncoct( frm.Frm(), ComMask_, useMass_, origin_ );
    Image::Nonortho( frm.ModifyFrm(), origin_, fcom, offset_,
                     frm.Frm().BoxCrd().UnitCell(), frm.Frm().BoxCrd().FracCell(), truncoct_,
                     *imageList_);
  }
  return Action::MODIFY_COORDS;
}
