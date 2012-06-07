#include "Action_AutoImage.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"

// CONSTRUCTOR
Action_AutoImage::Action_AutoImage() :
  origin_(false),
  ortho_(false),
  center_(false),
  triclinic_(OFF)
{}

// Action_AutoImage::init()
/** Usage: autoimage <mask> | anchor <mask> [fixed <fmask>] [mobile <mmask>]
  *                  [origin] [familiar | triclinic]
  */
int Action_AutoImage::init() {
  // Get keywords
  origin_ = actionArgs.hasKey("origin");
  center_ = actionArgs.hasKey("center");
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
  if (center_)
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

  return 0;
}

// Action_AutoImage::SetupAtomRanges()
/** Based on the given atom mask expression determine what molecules are
  * selected by the mask.
  * \return A list of atom pairs that mark the beginning and end of each
  *         selected molecule.
  */
Action_AutoImage::pairList Action_AutoImage::SetupAtomRanges( std::string const& maskexpr )
{
  pairList imageList;
  AtomMask Mask1( maskexpr.c_str() );

  if (currentParm->SetupCharMask( Mask1 )) return imageList;
  if (Mask1.None()) return imageList;
  for (Topology::mol_iterator mol = currentParm->MolStart();
                              mol != currentParm->MolEnd(); mol++)
  {
    int firstAtom = (*mol).BeginAtom();
    int lastAtom = (*mol).EndAtom();
    // Check that each atom in the range is in Mask1
    bool rangeIsValid = true;
    for (int atom = firstAtom; atom < lastAtom; ++atom) {
      if (!Mask1.AtomInCharMask(atom)) {
        rangeIsValid = false;
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

// Action_AutoImage::setup()
int Action_AutoImage::setup() {
  bool fixedauto = false;
  bool mobileauto = false;

  // Determine Box info
  if (currentParm->BoxType()==Box::NOBOX) {
    mprintf("Warning: Image::setup: Parm %s does not contain box information.\n",
            currentParm->c_str());
    return 1;
  }
  ortho_ = false;
  if (currentParm->BoxType()==Box::ORTHO && triclinic_==OFF) ortho_=true;
  // If box is originally truncated oct and not forcing triclinic, 
  // turn familiar on.
  if (currentParm->BoxType()==Box::TRUNCOCT && triclinic_!=FORCE && triclinic_!=FAMILIAR) {
    mprintf("\tOriginal box is truncated octahedron, turning on 'familiar'.\n");
    triclinic_=FAMILIAR;
  }

  // Set up anchor region
  if (!anchor_.empty()) {
    anchorList_ = SetupAtomRanges( anchor_ );
  } else {
    anchorList_.clear();
    anchorList_.push_back( currentParm->Mol(0).BeginAtom() );
    anchorList_.push_back( currentParm->Mol(0).EndAtom() );
  }
  if (anchorList_.empty() || anchorList_.size() > 2) {
    mprinterr("Error: Anchor mask [%s] corresponds to %zu mols, should only be 1.\n",
              anchor_.c_str(), anchorList_.size() / 2);
    return 1;
  }
  // Set up mask for centering anchor
  anchorMask_.AddAtomRange( anchorList_[0], anchorList_[1] );
  int anchormolnum = (*currentParm)[ anchorList_[0] ].Mol();
  mprintf("\tAnchor molecule is %i\n", anchormolnum+1);
  // Set up fixed region
  if (!fixed_.empty()) 
    fixedList_ = SetupAtomRanges( fixed_ );
  else 
    fixedauto = true;
  // Set up mobile region
  if (!mobile_.empty())
    mobileList_ = SetupAtomRanges( mobile_ );
  else
    mobileauto = true;
  // Automatic search through molecules for fixed/mobile
  if (fixedauto || mobileauto) {
    int molnum = 0;
    for (Topology::mol_iterator mol = currentParm->MolStart();
                                mol != currentParm->MolEnd(); mol++)
    {
      // Skip the anchor molecule
      if (molnum != anchormolnum) { 
        // Solvent and 1 atom molecules (prob. ions) go in mobile list,
        // everything else into fixed list.
        if ( (*mol).IsSolvent() || (*mol).NumAtoms() == 1 ) {
          if (mobileauto) {
            mobileList_.push_back( (*mol).BeginAtom() );
            mobileList_.push_back( (*mol).EndAtom()   );
          }
        } else {
          if (fixedauto) {
            fixedList_.push_back( (*mol).BeginAtom() );
            fixedList_.push_back( (*mol).EndAtom()   );
          }
        }
      }
      ++molnum;
    }
  }
  // DEBUG: Print fixed and mobile lists
  mprintf("\tThe following molecules are fixed to anchor:");
  for (pairList::iterator atom = fixedList_.begin(); 
                          atom != fixedList_.end(); atom += 2)
    mprintf(" %i", (*currentParm)[ *atom ].Mol()+1 );
  mprintf("\n\t%zu molecules are mobile.\n", mobileList_.size() / 2 );
  //mprintf("\tThe following molecules are mobile:\n");
  //for (pairList::iterator atom = mobileList_.begin(); 
  //                        atom != mobileList_.end(); atom += 2)
  //  mprintf("\t\t%i\n", (*currentParm)[ *atom ].Mol()+1 );

  return 0;
}

// Action_AutoImage::action()
int Action_AutoImage::action() {
  double center[3], ucell[9], recip[9], fixedcenter[3], framecenter[3];

  // Center w.r.t. anchor
  currentFrame->Center( anchorMask_, origin_, useMass_);
  // Determine whether anchor center is at box center or coordinate origin
  if (origin_) {
    center[0] = 0;
    center[1] = 0;
    center[2] = 0;
  } else {
    center[0] = currentFrame->BoxX() / 2;
    center[1] = currentFrame->BoxY() / 2;
    center[2] = currentFrame->BoxZ() / 2;
  }

  // Create copy of currentFrame
  // TODO: This is relatively expensive. Change so that potential
  //       imaging translation can be calcd and not actually applied.
  Frame fixedFrame( *currentFrame );

  // Image everything in fixedFrame according to fixedList, and everything
  // in currentFrame according to mobileList. Always use molecule COM when
  // imaging fixedList.
  if (ortho_) {
    fixedFrame.ImageOrtho(origin_, true, useMass_, fixedList_);
    currentFrame->ImageOrtho(origin_, center_, useMass_, mobileList_);
  } else {
    // TODO: Move BoxToRecip call out of ImageNonortho
    currentFrame->BoxToRecip(ucell, recip);
    fixedFrame.ImageNonortho(origin_, NULL, (triclinic_==FAMILIAR),
                             true, useMass_, fixedList_);
    currentFrame->ImageNonortho(origin_, NULL, (triclinic_==FAMILIAR),
                                center_, useMass_, mobileList_);
  }  

  // Determine if fixed position is closer to anchor center in currentFrame or fixedFrame
  for (pairList::iterator atom1 = fixedList_.begin();
                          atom1 != fixedList_.end(); ++atom1)
  {
    int firstAtom = *atom1;
    ++atom1;
    int lastAtom = *atom1;
    if (useMass_) {
      currentFrame->CenterOfMass(framecenter, firstAtom, lastAtom);
      fixedFrame.CenterOfMass(fixedcenter, firstAtom, lastAtom);
    } else {
      currentFrame->GeometricCenter(framecenter, firstAtom, lastAtom);
      fixedFrame.GeometricCenter(fixedcenter, firstAtom, lastAtom);
    }
    double framedist2 = DIST2_NoImage( center, framecenter );
    double fixeddist2 = DIST2_NoImage( center, fixedcenter );
    //mprintf("DBG: [%5i] Fixed @%i-%i frame dist2=%lf, imaged dist2=%lf\n", frameNum,
    //        firstAtom+1, lastAtom+1,
    //        framedist2, fixeddist2);
    if (fixeddist2 < framedist2) {
      // Imaging these atoms moved them closer to anchor. Update coords in currentFrame. 
      for (int idx = firstAtom*3; idx < lastAtom*3; ++idx)
        (*currentFrame)[idx] = fixedFrame[idx];
    }
  }
    
  return 0;
}

