#include <cmath> // modf
#include "Action_AutoImage.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "Dist_Imaged.h"
#include "ImageRoutines.h"
#include "CharMask.h"
#include "Image_List_Unit.h"

// CONSTRUCTOR
Action_AutoImage::Action_AutoImage() :
  debug_(0),
  mode_(UNSPECIFIED),
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
  mprintf("\t[{<mask> | anchor <mask> [fixed <fmask>] [mobile <mmask>]}]\n"
          "\t[origin] [firstatom] [{familiar|triclinic}] [moveanchor]\n"
          "\t[mode {bydist|byvec}]\n"
          "  Automatically center and image periodic trajectory.\n"
          "  The \"anchor\" molecule (default the first molecule) will be centered;\n"
          "  all \"fixed\" molecules will be imaged only if imaging brings them closer\n"
          "  to the \"anchor\" molecule; default for \"fixed\" molecules is all\n"
          "  non-solvent non-ion molecules. All other molecules (referred to as\n"
          "  \"mobile\") will be imaged freely.\n"
          "  If 'moveanchor' is specified the anchor point will be set to the\n"
          "  previous \"fixed\" molecule; this is only expected to work well\n"
          "  when \"fixed\" molecules that are sequential are also geometrically\n"
          "  close.\n"
          "  The 'mode' keyword determines how \"fixed\" molecules will be treated.\n"
          "  If 'bydist' (the default), \"fixed\" molecules will use the image closest\n"
          "  to the \"anchor\" molecule. If 'byvec', \"fixed\" molecules will use the image\n"
          "  closest to their orientation with respect to the anchor in the first frame.\n"
         );
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
  std::string modestr = actionArgs.GetStringKey("mode");
  if (modestr.empty()) {
    // Default mode
    mode_ = BY_DISTANCE;
  } else {
    if (modestr == "bydist")
      mode_ = BY_DISTANCE;
    else if (modestr == "byvec") {
      mode_ = BY_VECTOR;
      mprintf("Warning: 'mode byvec' is still being tested. Check results carefully.\n");
    } else {
      mprinterr("Error: '%s' is not a recognized autoimage mode.\n", modestr.c_str());
      return Action::ERR;
    }
  }
  // Get mask expression for anchor if none yet specified
  if (anchor_.empty())  
    anchor_ = actionArgs.GetMaskNext();
  RefVecs_.clear();

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
  if (mode_ == BY_DISTANCE) {
    mprintf("\tAuto-imaging using molecule distances.\n");
    if (movingAnchor_)
      mprintf("\tWhen imaging fixed molecules anchor will be set to previous fixed molecule.\n");
  } else if (mode_ == BY_VECTOR) {
    mprintf("\tAuto-imaging fixed molecules using molecule vectors.\n");
    if (movingAnchor_)
      mprintf("\tFor the first frame the anchor will be set to previous fixed molecule\n"
              "\t  when imaging fixed molecules.\n");
  }

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

  if (mode_ == BY_VECTOR)
    RefVecs_.reserve( fixedList_->nEntities() );

  return Action::OK;
}

// Action_AutoImage::DoAction()
Action::RetType Action_AutoImage::DoAction(int frameNum, ActionFrame& frm) {
  Action::RetType ret = Action::ERR;
  switch (mode_) {
    case BY_DISTANCE : ret = autoimage_by_distance(frameNum, frm); break;
    case BY_VECTOR   : ret = autoimage_by_vector(frameNum, frm); break;
    case UNSPECIFIED : ret = Action::ERR;
  }
  return ret;
}

/** \return center of given unit. */
Vec3 Action_AutoImage::unit_center(Varray const& fracCoords, Unit const& unit)
{
  Vec3 vcenter(0.0);
  unsigned int natoms = 0;
  for (Unit::const_iterator seg = unit.segBegin(); seg != unit.segEnd(); ++seg)
  {
    for (int at = seg->Begin(); at != seg->End(); ++at) {
      vcenter += fracCoords[at];
      natoms++;
    }
  }
  vcenter /= (double)natoms;
  return vcenter;
}

/** Translate given unit. */
void Action_AutoImage::translate_unit(Varray& fracCoords, Vec3 const& trans, Unit const& unit)
{
  for (Unit::const_iterator seg = unit.segBegin(); seg != unit.segEnd(); ++seg)
  {
    for (int at = seg->Begin(); at != seg->End(); ++at) {
      fracCoords[at] += trans;
    }
  }
}

/** \return frac coord vector needed to properly image the unit. */
Vec3 Action_AutoImage::calc_frac_image_vec(Vec3 const& delta_frac, bool& need_to_move) {
  Vec3 ivec_frac(0.0);
  need_to_move = false;

  for (int idx = 0; idx != 3; idx++) {
    double abs_dval = delta_frac[idx];
    if (abs_dval < 0.0) abs_dval = -abs_dval;
    if (abs_dval > 0.5) {
      // Vector has shifted more than half a box length.
      need_to_move = true;
      double currentVal = delta_frac[idx];
      if (delta_frac[idx] > 0.0) {
        //increment = -1.0;
        while (currentVal > 0.5) {
          currentVal -= 1.0;
          ivec_frac[idx] -= 1.0;
        }
      } else {
        //increment = 1.0;
        while (currentVal < -0.5) {
          currentVal += 1.0;
          ivec_frac[idx] += 1.0;
        }
      }
    } // END more than half box length traveled
  } // END loop over xyz idx
  return ivec_frac;
}

/** Calculate vector needed to wrap molecule back in the primary cell.
  * \return vector needed to image the molecule in frac space
  * \param min minimum frac coords
  * \param max maximum frac coords
  * \param pos current molecule position in frac coords
  * \param need_to_move set to true if returned vector is not zero.
  */
Vec3 Action_AutoImage::wrap_frac(Vec3 const& min, Vec3 const& max, Vec3 const& pos, bool& need_to_move) {
  Vec3 ivec_frac(0.0);
  need_to_move = false;

  //Vec3 min = center - 0.5;
  //Vec3 max = center + 0.5;

  Vec3 current = pos;
  for (int idx = 0; idx != 3; idx++) {
    while (current[idx] < min[idx]) {
      current[idx] += 1.0;
      ivec_frac[idx] += 1.0;
      need_to_move = true;
    }
    while (current[idx] > max[idx]) {
      current[idx] -= 1.0;
      ivec_frac[idx] -= 1.0;
      need_to_move = true;
    }
  }
  return ivec_frac;
}

/** Center the anchor molecule. */
Vec3 Action_AutoImage::center_anchor_molecule(ActionFrame& frm, bool is_ortho, bool use_ortho) const {
  Box const& box = frm.Frm().BoxCrd();
  //bool is_ortho = frm.Frm().BoxCrd().Is_X_Aligned_Ortho();
  //bool use_ortho = (is_ortho && triclinic_ == OFF);
  // Store anchor point in fcom for now.
  Vec3 fcom;
  if (useMass_)
    fcom = frm.Frm().VCenterOfMass( anchorMask_ );
  else
    fcom = frm.Frm().VGeometricCenter( anchorMask_ );
  // Determine translation to anchor point, store in fcom.
  // Anchor center will be in anchorcenter.
  Vec3 anchorcenter;
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
  return anchorcenter;
}

/** Autoimage molecules using reference vectors to anchor. */
Action::RetType Action_AutoImage::autoimage_by_vector(int frameNum, ActionFrame& frm) {
  Frame const& frameIn = frm.Frm();
  Box const& box = frameIn.BoxCrd();
  bool is_ortho = box.Is_X_Aligned_Ortho();
  bool use_ortho = (is_ortho && triclinic_ == OFF);

//  mprintf("DEBUG: ---------- autoimage_by_vector %i\n", frameNum);
  // We need the first frame to be properly imaged. Always use by distance first.
  if (RefVecs_.empty()) {
    Action::RetType ret = autoimage_by_distance(frameNum, frm);
    if (ret != Action::MODIFY_COORDS) {
      mprinterr("Error: autoimage by vector initial re-image by distance failed.\n");
      return ret;
    }
  } else {
    // Just move the anchor to the center
    center_anchor_molecule(frm, is_ortho, use_ortho);
  }

  // Convert everything to fractional coords.
  Varray fracCoords_;
  fracCoords_.reserve( frameIn.Natom() );
  for (int at = 0; at != frameIn.Natom(); ++at)
    fracCoords_.push_back( frameIn.BoxCrd().FracCell() * Vec3(frameIn.XYZ(at)) );

  // Calculate anchor center in fractional space
  Vec3 anchor_center_frac( 0.0 );
  for (AtomMask::const_iterator at = anchorMask_.begin(); at != anchorMask_.end(); ++at)
    anchor_center_frac += fracCoords_[*at];
  anchor_center_frac /= (double)anchorMask_.Nselected();
//  anchor_center_frac.Print("anchor_center_frac"); // DEBUG

  if (RefVecs_.empty()) {
  // If this is the first frame, save reference vectors to center
//    mprintf("DEBUG: Populating reference vectors.\n");
    for (Image::List_Unit::const_iterator it = fixedList_->begin();
                                          it != fixedList_->end(); ++it)
    {
      // DEBUG
//      for (Unit::const_iterator seg = it->segBegin(); seg != it->segEnd(); ++seg) {
//        mprintf("DEBUG: Fixed unit %li segment %li : %i to %i\n",
//                it - fixedList_->begin(), seg - it->segBegin(),
//                seg->Begin(), seg->End());
//      }
      Vec3 fixed_unit_center_frac = unit_center(fracCoords_, *it);
//      fixed_unit_center_frac.Print("fixed_unit_center_frac"); // DEBUG
      // Vector from fixed unit back to the anchor
      RefVecs_.push_back( fixed_unit_center_frac - anchor_center_frac );
//      RefVecs_.back().Print("RefVec"); // DEBUG
    } // END loop over fixed entities
  } else {
    // Not the first frame. Compare reference vectors to center
//    mprintf("DEBUG: Comparing reference vectors.\n");
    Varray::const_iterator refVec = RefVecs_.begin();
    for (Image::List_Unit::const_iterator it = fixedList_->begin();
                                          it != fixedList_->end(); ++it, ++refVec)
    {
      Vec3 fixed_unit_center_frac = unit_center(fracCoords_, *it);
//      fixed_unit_center_frac.Print("fixed_unit_center_frac"); // DEBUG
      Vec3 anchor_to_fixed_frac = fixed_unit_center_frac - anchor_center_frac;
//      anchor_to_fixed_frac.Print(  "anchor_to_fixed_frac  "); // DEBUG
//      refVec->Print(               "currentref            "); // DEBUG
      Vec3 delta_frac = anchor_to_fixed_frac - *refVec;
//      delta_frac.Print(            "delta_frac            "); // DEBUG
      bool need_to_move;
      Vec3 image_vec = calc_frac_image_vec( delta_frac, need_to_move );
//      mprintf("\tNeed to move= %i\n", (int)need_to_move); // DEBUG
//      image_vec.Print(             "image_vec             "); // DEBUG
      // Test that image_vec would do a good job moving the fixed unit
//      Vec3 test_vec = fixed_unit_center_frac + image_vec; // DEBUG
//      test_vec.Print(              "test_vec              "); // DEBUG
      // Move the unit if needed
      if (need_to_move) {
        translate_unit(fracCoords_, image_vec, *it);
        // Test that the image worked
//        test_vec = unit_center(fracCoords_, *it); // DEBUG
//        test_vec.Print(            "after imaging         "); // DEBUG
      }
//      mprintf("\t--------------------\n"); // DEBUG
    }

    // Mobile molecules
    Vec3 minvec = anchor_center_frac - 0.5;
    Vec3 maxvec = anchor_center_frac + 0.5;
//    minvec.Print("MinVec"); // DEBUG
//    maxvec.Print("MaxVec"); // DEBUG
    for (Image::List_Unit::const_iterator it = mobileList_->begin();
                                          it != mobileList_->end(); ++it, ++refVec)
    {
      Vec3 mobile_unit_center_frac = unit_center(fracCoords_, *it);
      bool need_to_move;
      Vec3 image_vec = wrap_frac( minvec, maxvec, mobile_unit_center_frac, need_to_move );
      if (need_to_move) {
        translate_unit(fracCoords_, image_vec, *it);
      }
      // Convert to familiar truncated octahedron shape.
      if (truncoct_) {
        // Center of mobile unit after imaging
        Vec3 translated_coord_frac = mobile_unit_center_frac + image_vec;
        // Closest distance of mobile unit center to anchor center
        int ixyz[3];
        Cpptraj::Dist2_Imaged_Frac( translated_coord_frac, anchor_center_frac, box.UnitCell(), box.FracCell(), ixyz );
        if (ixyz[0] != 0 || ixyz[1] != 0 || ixyz[2] != 0) {
          // The reflection is closer to the center, so move the mobile unit.
          translate_unit(fracCoords_, Vec3(ixyz[0], ixyz[1], ixyz[2]), *it);
          //boxTransOut += ucell.TransposeMult( ixyz );
          //if (debug > 2)
          //  mprintf("  IMAGING, FAMILIAR OFFSETS ARE %i %i %i\n", ixyz[0], ixyz[1], ixyz[2]);
        }
      }
    }
    // Convert back to Cartesian
    for (int at = 0; at != frameIn.Natom(); ++at)
      frm.ModifyFrm().SetXYZ(at, box.UnitCell().TransposeMult( fracCoords_[at] ));
  }

  return MODIFY_COORDS;
}

/** Original autoimage algorithm by distance. */
Action::RetType Action_AutoImage::autoimage_by_distance(int frameNum, ActionFrame& frm) {
  Vec3 fcom;
  Vec3 bp, bm, offset(0.0);
  Vec3 Trans, framecenter, imagedcenter;

  Box const& box = frm.Frm().BoxCrd();
  bool is_ortho = frm.Frm().BoxCrd().Is_X_Aligned_Ortho();
  bool use_ortho = (is_ortho && triclinic_ == OFF);
  Vec3 anchorcenter = center_anchor_molecule(frm, is_ortho, use_ortho);
/*
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
  frm.ModifyFrm().Translate(fcom);*/

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
