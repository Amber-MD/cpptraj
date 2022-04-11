#include <cmath> // modf
#include "Action_AutoImage.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "ImageRoutines.h"
#include "CharMask.h"
#include "Image_List_Unit.h"

// CONSTRUCTOR
Action_AutoImage::Action_AutoImage() :
  debug_(0),
  origin_(false),
  usecom_(true),
  truncoct_(false),
  useMass_(false),
  movingAnchor_(false),
  triclinic_(OFF),
  fixedList_(0),
  mobileList_(0)
{}

/** DESTRUCTOR */
Action_AutoImage::~Action_AutoImage() {
  if (fixedList_ != 0) delete fixedList_;
  if (mobileList_ != 0) delete mobileList_;
}

void Action_AutoImage::Help() const {
  mprintf("\t[<mask> | anchor <mask> [fixed <fmask>] [mobile <mmask>]]\n"
          "\t[origin] [firstatom] [familiar | triclinic] [moveanchor]\n"
          "  Automatically center and image periodic trajectory.\n"
          "  The \"anchor\" molecule (default the first molecule) will be centered;\n"
          "  all \"fixed\" molecules will be imaged only if imaging brings them closer\n"
          "  to the \"anchor\" molecule; default for \"fixed\" molecules is all\n"
          "  non-solvent non-ion molecules. All other molecules (referred to as\n"
          "  \"mobile\") will be imaged freely.\n"
          "  If 'moveanchor' is specified the anchor point will be set to the\n"
          "  previous \"fixed\" molecule; this is only expected to work well\n"
          "  when \"fixed\" molecules that are sequential are also geometrically\n"
          "  close.\n");
}

// Action_AutoImage::Init()
Action::RetType Action_AutoImage::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  // Get keywords
  origin_ = actionArgs.hasKey("origin");
  usecom_ = !actionArgs.hasKey("firstatom");
  movingAnchor_ = actionArgs.hasKey("moveanchor");
  if (actionArgs.hasKey("familiar")) triclinic_ = FAMILIAR;
  if (actionArgs.hasKey("triclinic")) triclinic_ = FORCE;
  anchor_ = actionArgs.GetStringKey("anchor");
  fixed_  = actionArgs.GetStringKey("fixed");
  mobile_ = actionArgs.GetStringKey("mobile");
  // Get mask expression for anchor if none yet specified
  if (anchor_.empty())  
    anchor_ = actionArgs.GetMaskNext();

  mprintf("    AUTOIMAGE: To");
  if (origin_)
    mprintf(" origin");
  else
    mprintf(" box center");
  mprintf(" based on");
  if (usecom_)
    mprintf(" center of mass");
  else
    mprintf(" first atom position");
  if (!anchor_.empty())
    mprintf(", anchor mask is [%s]\n", anchor_.c_str());
  else
    mprintf(", anchor is first molecule.\n");
  if (!fixed_.empty())
    mprintf("\tAtoms in mask [%s] will be fixed to anchor region.\n", fixed_.c_str());
  if (!mobile_.empty())
    mprintf("\tAtoms in mask [%s] will be imaged independently of anchor region.\n",
            mobile_.c_str());
  if (movingAnchor_)
    mprintf("\tWhen imaging fixed molecules anchor will be set to previous fixed molecule.\n");

  return Action::OK;
}

// Action_AutoImage::Setup()
Action::RetType Action_AutoImage::Setup(ActionSetup& setup) {
  bool fixedauto = false;
  bool mobileauto = false;

  if (setup.Top().Nmol() < 1) {
    mprintf("Warning: Topology %s does not contain molecule information\n", setup.Top().c_str());
    return Action::SKIP;
  }
  // Determine Box info
  if (!setup.CoordInfo().TrajBox().HasBox()) {
    mprintf("Warning: Topology %s does not contain box information.\n", setup.Top().c_str());
    return Action::SKIP;
  }
  // If box is originally truncated oct and not forcing triclinic, 
  // turn familiar on.
  if (triclinic_ != FORCE && triclinic_ != FAMILIAR && setup.CoordInfo().TrajBox().CellShape() == Box::OCTAHEDRAL)
  {
    mprintf("\tOriginal box is truncated octahedron, turning on 'familiar'.\n");
    triclinic_ = FAMILIAR;
  }

  // Set up anchor mask
  anchorMask_.ResetMask();
  int anchormolnum = -1;
  if (!anchor_.empty()) {
    // Anchor molecule/region specified
    mprintf("\tAnchoring on atoms selected by mask '%s'\n", anchor_.c_str());
    if (anchorMask_.SetMaskString( anchor_ )) return Action::ERR;
    if ( setup.Top().SetupIntegerMask( anchorMask_ ) ) return Action::ERR;
    anchorMask_.MaskInfo();
    if (anchorMask_.None()) {
      mprinterr("Error: No atoms selected for anchor.\n");
      return Action::ERR;
    }
    // If mask pertains to only one molecule, do not include that molecule
    // in the fixed region.
    std::vector<int> molnums = setup.Top().MolnumsSelectedBy( anchorMask_ );
    if (molnums.size() == 1)
      anchormolnum = molnums.front();
    if (anchormolnum != -1)
      mprintf("\tMask [%s] corresponds to molecule %i\n",
              anchorMask_.MaskString(), anchormolnum+1);
  } else {
    // No anchor specified. Use first molecule as anchor.
    anchormolnum = 0;
    mprintf("\tUsing first molecule as anchor.\n");
    anchorMask_.AddUnit( setup.Top().Mol(0).MolUnit() );
  }

  if (fixedList_ != 0) delete fixedList_;
  if (mobileList_ != 0) delete mobileList_;
  // Set up fixed region
  // NOTE: Always use molecule center when imaging fixed list
  if (!fixed_.empty())
    fixedList_ = (Image::List_Unit*)
                 Image::CreateImageList(setup.Top(), Image::BYMOL, fixed_, useMass_, true); 
  else { 
    fixedauto = true;
    fixedList_ = (Image::List_Unit*)
                 Image::CreateImageList(Image::BYMOL, useMass_, true);
  }
  if (fixedList_ == 0) {
    mprinterr("Internal Error: Could not allocate fixed list.\n");
    return Action::ERR;
  }
  // Set up mobile region
  if (!mobile_.empty())
    mobileList_ = (Image::List_Unit*)
                  Image::CreateImageList(setup.Top(), Image::BYMOL, mobile_, useMass_, usecom_);
  else {
    mobileauto = true;
    mobileList_ = (Image::List_Unit*)
                  Image::CreateImageList(Image::BYMOL, useMass_, usecom_);
  }
  if (mobileList_ == 0) {
    mprinterr("Internal Error: Could not allocate mobile list.\n");
    return Action::ERR;
  }
  // Automatic search through molecules for fixed/mobile
  if (fixedauto || mobileauto) {
    int molnum = 0;
    for (Topology::mol_iterator mol = setup.Top().MolStart();
                                mol != setup.Top().MolEnd(); mol++)
    {
      // Skip the anchor molecule
      if (molnum != anchormolnum) { 
        // Solvent and 1 atom molecules (prob. ions) go in mobile list,
        // everything else into fixed list.
        if ( mol->IsSolvent() || mol->NumAtoms() == 1 ) {
          if (mobileauto) {
            mobileList_->AddUnit( mol->MolUnit() );
          }
        } else {
          if (fixedauto) {
            fixedList_->AddUnit( mol->MolUnit() );
          }
        }
      }
      ++molnum;
    }
  }
  // Print fixed and mobile lists
  if (!fixedList_->empty()) {
    mprintf("\t%u molecules are fixed to anchor:", fixedList_->nEntities());
    for (Image::List_Unit::const_iterator it = fixedList_->begin();
                                it != fixedList_->end(); ++it)
      mprintf(" %i", setup.Top()[ it->Front() ].MolNum()+1 );
    mprintf("\n");
  }
  mprintf("\t%u molecules are mobile.\n", mobileList_->nEntities() );
  if (debug_ > 1) {
    mprintf("\tThe following molecules are mobile:\n");
    for (Image::List_Unit::const_iterator it = mobileList_->begin();
                                it != mobileList_->end(); ++it)
      mprintf(" %i\n", setup.Top()[ it->Front() ].MolNum()+1 );
    mprintf("\n");
  }

  truncoct_ = (triclinic_==FAMILIAR);

  return Action::OK;
}

// Action_AutoImage::DoAction()
Action::RetType Action_AutoImage::DoAction(int frameNum, ActionFrame& frm) {
  Vec3 fcom;
  Vec3 bp, bm, offset(0.0);
  Vec3 Trans, framecenter, imagedcenter, anchorcenter;

  Box const& box = frm.Frm().BoxCrd();
  bool is_ortho = frm.Frm().BoxCrd().Is_X_Aligned_Ortho();
  bool use_ortho = (is_ortho && triclinic_ == OFF);
  // Store anchor point in fcom for now.
  if (useMass_)
    fcom = frm.Frm().VCenterOfMass( anchorMask_ );
  else
    fcom = frm.Frm().VGeometricCenter( anchorMask_ );
  // Determine translation to anchor point, store in fcom.
  // Anchor center will be in anchorcenter.
  if (origin_) {
    // Center is coordinate origin (0,0,0)
    fcom.Neg();
    anchorcenter.Zero();
  } else {
    // Center on box center
    if (is_ortho || truncoct_)
      // Center is box xyz over 2
      anchorcenter = box.Center();
    else
      // Center in frac coords is (0.5,0.5,0.5)
      anchorcenter = box.UnitCell().TransposeMult(Vec3(0.5));
    fcom = anchorcenter - fcom;
  }
  frm.ModifyFrm().Translate(fcom);

  // Setup imaging, and image everything in current Frame
  // according to mobileList_.
  if (is_ortho) {
    if (Image::SetupOrtho(box, bp, bm, origin_)) {
      mprintf("Warning: Frame %i imaging failed, box lengths are zero.\n",frameNum+1);
      // TODO: Return OK for now so next frame is tried; eventually indicate SKIP?
      return Action::OK; // FIXME return MODIFY_COORDS instead?
    }
    Image::Ortho(frm.ModifyFrm(), bp, bm, offset, *mobileList_);
  } else {
    if (truncoct_)
      fcom = Image::SetupTruncoct( frm.Frm(), 0, useMass_, origin_ );
    Image::Nonortho(frm.ModifyFrm(), origin_, fcom, offset, box.UnitCell(), box.FracCell(), truncoct_,
                    *mobileList_);
  }

  if (movingAnchor_) {
    // TODO I think the way the translation is calculated here is robust and
    //      more efficient than the !movingAnchor_ case but more testing is
    //      needed.
    // Loop over fixed molecules
    for (unsigned int idx = 0; idx != fixedList_->nEntities(); ++idx)
    {
      framecenter = fixedList_->GetCoord(idx, frm.Frm());

      // Determine distance in terms of box lengths
      if (use_ortho) {
        // Determine direction from molecule to anchor
        Vec3 delta = anchorcenter - framecenter;
        //mprintf("DEBUG: anchorcenter - framecenter = %g %g %g\n", delta[0], delta[1], delta[2]);
        Vec3 minTrans( floor(delta[0]/box.Param(Box::X)+0.5)*box.Param(Box::X),
                       floor(delta[1]/box.Param(Box::Y)+0.5)*box.Param(Box::Y),
                       floor(delta[2]/box.Param(Box::Z)+0.5)*box.Param(Box::Z) );
        Vec3 minImage = framecenter + minTrans;
        //mprintf("DBG: %5i %3u %6i %6i {%8.2f %8.2f %8.2f}\n",
        //        frameNum, (atom1-fixedList_.begin())/2, firstAtom+1, lastAtom,
        //        minTrans[0], minTrans[1], minTrans[2]);
        // Move atoms closer to anchor. Update coords in currentFrame.
        fixedList_->DoTranslation(frm.ModifyFrm(), idx, minTrans);
        // New anchor is previous fixed mol
        anchorcenter = minImage;
      } else {
        Vec3 newAnchor = framecenter;
        Trans = Image::Nonortho(framecenter, truncoct_, origin_, box.UnitCell(), box.FracCell(), fcom, -1.0);
        // If molecule was imaged, determine whether imaged position is closer to anchor.
        if (Trans[0] != 0 || Trans[1] != 0 || Trans[2] != 0) {
          imagedcenter = framecenter + Trans;
          double framedist2 = DIST2_NoImage( anchorcenter, framecenter );
          double imageddist2 = DIST2_NoImage( anchorcenter, imagedcenter );
          //mprintf("DBG: %5i %3u %6i %6i {%8.2f %8.2f %8.2f}"
          //        " frame dist2=%6.2f, imaged dist2=%6.2f\n",
          //        frameNum, (atom1-fixedList_.begin())/2, firstAtom+1, lastAtom,
          //        Trans[0], Trans[1], Trans[2], sqrt(framedist2), sqrt(imageddist2));
          if (imageddist2 < framedist2) {
            // Imaging these atoms moved them closer to anchor. Update coords in currentFrame.
            fixedList_->DoTranslation(frm.ModifyFrm(), idx, Trans);
            newAnchor = imagedcenter;
          }
        }
        anchorcenter = newAnchor;
      }
    }
  } else {
    // For each molecule defined by atom pairs in fixedList, determine if the
    // imaged position is closer to anchor center than the current position.
    // Always use molecule center when imaging fixedList.
    for (unsigned int idx = 0; idx != fixedList_->nEntities(); ++idx)
    {
      framecenter = fixedList_->GetCoord(idx, frm.Frm());
      
      // Determine if molecule would be imaged.
      if (use_ortho)
        Trans = Image::Ortho(framecenter, bp, bm, box);
      else
        Trans = Image::Nonortho(framecenter, truncoct_, origin_, box.UnitCell(), box.FracCell(), fcom, -1.0);
      // If molecule was imaged, determine whether imaged position is closer to anchor.
      if (Trans[0] != 0 || Trans[1] != 0 || Trans[2] != 0) {
        imagedcenter = framecenter + Trans;
        double framedist2 = DIST2_NoImage( anchorcenter, framecenter );
        double imageddist2 = DIST2_NoImage( anchorcenter, imagedcenter );
//        mprintf("DBG: %5i %3u %6i %6i {%8.2f %8.2f %8.2f} frame dist2=%6.2f, imaged dist2=%6.2f\n",
//                frameNum, (atom1-fixedList_.begin())/2, firstAtom+1, lastAtom,
//                Trans[0], Trans[1], Trans[2], sqrt(framedist2), sqrt(imageddist2));
        if (imageddist2 < framedist2) {
          // Imaging these atoms moved them closer to anchor. Update coords in currentFrame.
          fixedList_->DoTranslation(frm.ModifyFrm(), idx, Trans);
        }
      }
    }
  }

  return Action::MODIFY_COORDS;
}
