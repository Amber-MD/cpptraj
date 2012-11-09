// Action_Image 
#include "Action_Image.h"
#include "CpptrajStdio.h"
#include "ImageRoutines.h"

// CONSTRUCTOR
Action_Image::Action_Image() :
  imageMode_(BYMOL),
  ComMask_(NULL),
  origin_(false),
  center_(false),
  ortho_(false),
  useMass_(true),
  truncoct_(false),
  triclinic_(OFF),
  debug_(0)
{ } 

void Action_Image::Help() {
  mprintf("Expected call: image [origin] [center] [triclinic | familiar [com <mask>]] <mask>\n"); 
  mprintf("                     [ bymol | byres | byatom ]\n");
  mprintf("- origin: center at 0.0, 0.0, 0.0, otherwise center at box center.\n");
  mprintf("- center: Use center of mass for imaging, otherwise use first atom.\n");
  mprintf("- triclinic: Force imaging with triclinic code.\n");
  mprintf("- familiar: Image with triclinic code and shape into familiar trunc. oct. shape.\n");
  mprintf("- com <mask>: If familiar, center based on COM of atoms in mask, otherwise use\n");
  mprintf("              origin/box.\n");
  mprintf("- <mask>: Only image atoms in <mask>. If no mask given all atoms are imaged.\n");
}

// DESTRUCTOR
Action_Image::~Action_Image() {
  if (ComMask_!=NULL) delete ComMask_;
}

const char Action_Image::ImageModeString[3][9] = {
  "molecule", "residue", "atom"
};

// Action_Image::init()
Action::RetType Action_Image::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  // Get keywords
  origin_ = actionArgs.hasKey("origin");
  center_ = actionArgs.hasKey("center");
  if (actionArgs.hasKey("familiar")) triclinic_ = FAMILIAR;
  if (actionArgs.hasKey("triclinic")) triclinic_ = FORCE;
  if (actionArgs.hasKey("bymol"))
    imageMode_ = BYMOL;
  else if (actionArgs.hasKey("byres"))
    imageMode_ = BYRES;
  else if (actionArgs.hasKey("byatom"))
    imageMode_ = BYATOM;

  // Get Masks
  if (triclinic_ == FAMILIAR) {
    std::string maskexpr = actionArgs.GetStringKey("com");
    if (!maskexpr.empty()) {
      ComMask_ = new AtomMask();
      ComMask_->SetMaskString(maskexpr);
    }
  }
  Mask1_.SetMaskString(actionArgs.GetMaskNext());
  
  mprintf("    IMAGE: By %s to", ImageModeString[imageMode_]);
  if (origin_)
    mprintf(" origin");
  else
    mprintf(" box center");
  mprintf(" based on");
  if (center_)
    mprintf(" center of mass");
  else
    mprintf(" first atom position");
  mprintf(" using atoms in mask %s\n",Mask1_.MaskString());
  if (triclinic_ == FORCE)
    mprintf( "           Triclinic On.\n");
  else if (triclinic_ == FAMILIAR) {
    mprintf( "           Triclinic On, familiar shape");
    if (ComMask_!=NULL) 
      mprintf( " centering on atoms in mask %s", ComMask_->MaskString());
    mprintf(".\n");
  }

  return Action::OK;
}

/** Check that at least 1 atom in the range is in Mask1 */
void Action_Image::CheckRange(int firstAtom, int lastAtom) {
  bool rangeIsValid = false;
  for (int atom = firstAtom; atom < lastAtom; ++atom) {
    if (Mask1_.AtomInCharMask(atom)) {
      rangeIsValid = true; 
      break;
    }
  }
  if (rangeIsValid) {
    imageList_.push_back( firstAtom );
    imageList_.push_back( lastAtom );
  }
}

// Action_Image::setup()
/** Set Imaging up for this parmtop. Get masks etc.
  * currentParm is set in Action::Setup
  */
Action::RetType Action_Image::Setup(Topology* currentParm, Topology** parmAddress) {
  if ( imageMode_ == BYMOL || imageMode_ == BYRES ) {
    if ( currentParm->SetupCharMask( Mask1_ ) ) return Action::ERR;
  } else { // BYATOM
    if ( currentParm->SetupIntegerMask( Mask1_ ) ) return Action::ERR;
  }
  if (Mask1_.None()) {
    mprintf("Warning: Image::setup: Mask contains 0 atoms.\n");
    return Action::ERR;
  }

  if (currentParm->BoxType()==Box::NOBOX) {
    mprintf("Warning: Image::setup: Parm %s does not contain box information.\n",
            currentParm->c_str());
    return Action::ERR;
  }

  ortho_ = false;  
  if (currentParm->BoxType()==Box::ORTHO && triclinic_==OFF) ortho_=true;

  // If box is originally truncated oct and not forcing triclinic, 
  // turn familiar on.
  if (currentParm->BoxType()==Box::TRUNCOCT && triclinic_!=FORCE && triclinic_!=FAMILIAR) {
    mprintf("\tOriginal box is truncated octahedron, turning on 'familiar'.\n");
    triclinic_=FAMILIAR;
  }

  if (triclinic_ == FAMILIAR) {
    if (ComMask_!=NULL) {
      if ( currentParm->SetupIntegerMask( *ComMask_ ) ) return Action::ERR;
      if (ComMask_->None()) {
        mprintf("Warning: Image::setup: Mask for 'familiar com' contains no atoms.\n");
        return Action::ERR;
      }
      mprintf("\tcom: mask [%s] contains %i atoms.\n",ComMask_->MaskString(),ComMask_->Nselected());
    }
  }

  // Set up atom range for each entity to be imaged. 
  // Currently imaging by molecule only, so each pair will be the first and
  // last atom of each molecule. Check that all atoms between first and last
  // are actually in the mask.
  imageList_.clear();

  switch (imageMode_) {
    case BYMOL:
      imageList_.reserve( currentParm->Nmol()*2 );
      for (Topology::mol_iterator mol = currentParm->MolStart();
                                  mol != currentParm->MolEnd(); ++mol)
        CheckRange( (*mol).BeginAtom(), (*mol).EndAtom());
     break;
    case BYRES:
      imageList_.reserve( currentParm->Nres()*2 );
      for (int resnum = 0; resnum < currentParm->Nres(); ++resnum)
        CheckRange( currentParm->ResFirstAtom( resnum ), currentParm->ResLastAtom( resnum ) );
      break;
    case BYATOM:
      imageList_.reserve( currentParm->Natom()*2 );
      for (AtomMask::const_iterator atom = Mask1_.begin(); atom != Mask1_.end(); ++atom) {
        imageList_.push_back( *atom );
        imageList_.push_back( (*atom)+1 );
      }
      break;
  }
  mprintf("\tNumber of %ss to be imaged is %zu based on mask [%s]\n", 
           ImageModeString[imageMode_], imageList_.size()/2, Mask1_.MaskString());
 
  // DEBUG: Print all pairs
  if (debug_>0) {
    for (std::vector<int>::iterator ap = imageList_.begin();
                                    ap != imageList_.end(); ap+=2)
      mprintf("\t\tFirst-Last atom#: %i - %i\n", (*ap)+1, *(ap+1) );
  }

  // Truncoct flag
  truncoct_ = (triclinic_==FAMILIAR);

  return Action::OK;  
}

// Action_Image::action()
Action::RetType Action_Image::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  // Ortho
  Vec3 bp, bm;
  // Nonortho
  double ucell[9], recip[9];
  Vec3 fcom;
  
  if (ortho_) {
    SetupImageOrtho(*currentFrame, bp, bm, origin_);
    ImageOrtho(*currentFrame, bp, bm, center_, useMass_, imageList_);
  } else {
    currentFrame->BoxToRecip( ucell, recip );
    if (truncoct_)
      fcom = SetupImageTruncoct( *currentFrame, ComMask_, useMass_, origin_ );
    ImageNonortho( *currentFrame, origin_, fcom, ucell, recip, truncoct_,
                                center_, useMass_, imageList_);
  }
  return Action::OK;
} 
