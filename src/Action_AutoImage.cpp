#include "Action_AutoImage.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "ImageRoutines.h"

// CONSTRUCTOR
Action_AutoImage::Action_AutoImage() :
  debug_(0),
  origin_(false),
  ortho_(false),
  usecom_(true),
  truncoct_(false),
  useMass_(false),
  triclinic_(OFF)
{}

void Action_AutoImage::Help() const {
  mprintf("\t[<mask> | anchor <mask> [fixed <fmask>] [mobile <mmask>]]\n"
          "\t[origin] [firstatom] [familiar | triclinic]\n"
          "  Automatically center and image periodic trajectory.\n"
          "  The 'anchor' molecule (default the first molecule) will be centered;\n"
          "  all 'fixed' molecules will be imaged only if imaging brings them closer\n"
          "  to the 'anchor' molecule; default for 'fixed' molecules is all\n"
          "  non-solvent non-ion molecules. All other molecules (referred to as\n"
          "  'mobile') will be imaged freely.\n");
}

// Action_AutoImage::Init()
Action::RetType Action_AutoImage::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  // Get keywords
  origin_ = actionArgs.hasKey("origin");
  usecom_ = !actionArgs.hasKey("firstatom");
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

  return Action::OK;
}

// Action_AutoImage::SetupAtomRanges()
/** Based on the given atom mask expression determine what molecules are
  * selected by the mask. If a mask selects any part of a molecule the
  * entire molecule will be selected.
  * \return A list of atom pairs that mark the beginning and end of each
  *         selected molecule.
  */
Action_AutoImage::pairList
  Action_AutoImage::SetupAtomRanges(Topology const& currentParm, std::string const& maskexpr)
{
  pairList imageList;
  CharMask Mask1( maskexpr.c_str() );

  if (currentParm.SetupCharMask( Mask1 )) return imageList;
  if (Mask1.None()) return imageList;
  for (Topology::mol_iterator mol = currentParm.MolStart(); mol != currentParm.MolEnd(); mol++)
  {
    int firstAtom = mol->BeginAtom();
    int lastAtom = mol->EndAtom();
    bool rangeIsValid = false;
    // Check that any atom in the range is in Mask1
    for (int atom = firstAtom; atom < lastAtom; ++atom) {
      if (Mask1.AtomInCharMask(atom)) {
        rangeIsValid = true;
        break;
      }
    }
    if (rangeIsValid) {
      imageList.push_back( firstAtom );
      imageList.push_back( lastAtom );
    }
  }
  mprintf("\tMask [%s] corresponds to %zu molecules\n", Mask1.MaskString(), imageList.size()/2);
  return imageList;
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
  Box::BoxType boxType = setup.CoordInfo().TrajBox().Type();
  if (boxType == Box::NOBOX) {
    mprintf("Warning: Topology %s does not contain box information.\n", setup.Top().c_str());
    return Action::SKIP;
  }
  ortho_ = false;
  if (boxType == Box::ORTHO && triclinic_ == OFF) ortho_ = true;
  // If box is originally truncated oct and not forcing triclinic, 
  // turn familiar on.
  if (boxType == Box::TRUNCOCT && triclinic_ != FORCE && triclinic_ != FAMILIAR) {
    mprintf("\tOriginal box is truncated octahedron, turning on 'familiar'.\n");
    triclinic_ = FAMILIAR;
  }

  // Set up anchor mask
  anchorMask_.ResetMask();
  int anchormolnum = -1;
  if (!anchor_.empty()) {
    // Anchor molecule/region specified
    mprintf("\tAnchoring on atoms selected by mask '%s'\n", anchor_.c_str());
    anchorMask_.SetMaskString( anchor_ );
    if ( setup.Top().SetupIntegerMask( anchorMask_ ) ) return Action::ERR;
    anchorMask_.MaskInfo();
    if (anchorMask_.None()) {
      mprinterr("Error: No atoms selected for anchor.\n");
      return Action::ERR;
    }
    // If mask pertains to only one molecule, do not include that molecule
    // in the fixed region.
    AtomMask::const_iterator at = anchorMask_.begin();
    anchormolnum = setup.Top()[ *at ].MolNum();
    ++at;
    for (; at != anchorMask_.end(); ++at) {
      if ( setup.Top()[ *at ].MolNum() != anchormolnum ) {
        anchormolnum = 1;
        break;
      }
    }
    if (anchormolnum != -1)
      mprintf("\tMask [%s] corresponds to molecule %i\n",
              anchorMask_.MaskString(), anchormolnum+1);
  } else {
    // No anchor specified. Use first molecule as anchor.
    anchormolnum = 0;
    mprintf("\tUsing first molecule as anchor.\n");
    anchorMask_.AddAtomRange( setup.Top().Mol(0).BeginAtom(),
                              setup.Top().Mol(0).EndAtom()    );
  }

  // Set up fixed region
  if (!fixed_.empty()) 
    fixedList_ = SetupAtomRanges( setup.Top(), fixed_ );
  else { 
    fixedauto = true;
    fixedList_.clear();
  }
  // Set up mobile region
  if (!mobile_.empty())
    mobileList_ = SetupAtomRanges( setup.Top(), mobile_ );
  else {
    mobileauto = true;
    mobileList_.clear();
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
            mobileList_.push_back( mol->BeginAtom() );
            mobileList_.push_back( mol->EndAtom()   );
          }
        } else {
          if (fixedauto) {
            fixedList_.push_back( mol->BeginAtom() );
            fixedList_.push_back( mol->EndAtom()   );
          }
        }
      }
      ++molnum;
    }
  }
  // Print fixed and mobile lists
  if (!fixedList_.empty()) {
    mprintf("\t%zu molecules are fixed to anchor:", fixedList_.size() / 2);
    for (pairList::const_iterator atom = fixedList_.begin();
                                  atom != fixedList_.end(); atom += 2)
      mprintf(" %i", setup.Top()[ *atom ].MolNum()+1 );
    mprintf("\n");
  }
  mprintf("\t%zu molecules are mobile.\n", mobileList_.size() / 2 );
  if (debug_ > 1) {
    mprintf("\tThe following molecules are mobile:\n");
    for (pairList::const_iterator atom = mobileList_.begin();
                                  atom != mobileList_.end(); atom += 2)
      mprintf(" %i\n", setup.Top()[ *atom ].MolNum()+1 );
    mprintf("\n");
  }

  truncoct_ = (triclinic_==FAMILIAR);

  return Action::OK;
}

// Action_AutoImage::DoAction()
Action::RetType Action_AutoImage::DoAction(int frameNum, ActionFrame& frm) {
  Matrix_3x3 ucell, recip;
  Vec3 fcom;
  Vec3 bp, bm, offset(0.0);
  Vec3 Trans, framecenter, imagedcenter, anchorcenter;

  if (!ortho_) frm.Frm().BoxCrd().ToRecip(ucell, recip);
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
    if (ortho_ || truncoct_)
      // Center is box xyz over 2
      anchorcenter = frm.Frm().BoxCrd().Center();
    else
      // Center in frac coords is (0.5,0.5,0.5)
      anchorcenter = ucell.TransposeMult(Vec3(0.5));
    fcom = anchorcenter - fcom;
  }
  frm.ModifyFrm().Translate(fcom);

  // Setup imaging, and image everything in current Frame 
  // according to mobileList_. 
  if (ortho_) {
    if (Image::SetupOrtho(frm.Frm().BoxCrd(), bp, bm, origin_)) {
      mprintf("Warning: Frame %i imaging failed, box lengths are zero.\n",frameNum+1);
      // TODO: Return OK for now so next frame is tried; eventually indicate SKIP?
      return Action::OK; // FIXME return MODIFY_COORDS instead?
    }
    Image::Ortho(frm.ModifyFrm(), bp, bm, offset, usecom_, useMass_, mobileList_);
  } else {
    if (truncoct_)
      fcom = Image::SetupTruncoct( frm.Frm(), 0, useMass_, origin_ );
    Image::Nonortho(frm.ModifyFrm(), origin_, fcom, offset, ucell, recip, truncoct_,
                    usecom_, useMass_, mobileList_);
  }

  // For each molecule defined by atom pairs in fixedList, determine if the
  // imaged position is closer to anchor center than the current position.
  // Always use molecule center when imaging fixedList.
  for (pairList::const_iterator atom1 = fixedList_.begin();
                                atom1 != fixedList_.end(); atom1 += 2)
  {
    int firstAtom = *atom1;
    int lastAtom = *(atom1+1);
    if (useMass_) 
      framecenter = frm.Frm().VCenterOfMass(firstAtom, lastAtom);
    else
      framecenter = frm.Frm().VGeometricCenter(firstAtom, lastAtom);
    // Determine if molecule would be imaged.
    if (ortho_)
      Trans = Image::Ortho(framecenter, bp, bm, frm.Frm().BoxCrd());
    else
      Trans = Image::Nonortho(framecenter, truncoct_, origin_, ucell, recip, fcom, -1.0);
    // If molecule was imaged, determine whether imaged position is closer to anchor.
    if (Trans[0] != 0 || Trans[1] != 0 || Trans[2] != 0) {
      imagedcenter = framecenter + Trans;
      double framedist2 = DIST2_NoImage( anchorcenter, framecenter );
      double imageddist2 = DIST2_NoImage( anchorcenter, imagedcenter );
//      mprintf("DBG: [%5i] Fixed @%i-%i frame dist2=%12.4f, imaged dist2=%12.4f\n",
//              frameNum, firstAtom+1, lastAtom, framedist2, imageddist2);
      if (imageddist2 < framedist2) {
        // Imaging these atoms moved them closer to anchor. Update coords in currentFrame.
        frm.ModifyFrm().Translate(Trans, firstAtom, lastAtom);
      }
    }
  }
  return Action::MODIFY_COORDS;
}
