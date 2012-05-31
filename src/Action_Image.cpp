// Action_Image 
#include "Action_Image.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Image::Action_Image() :
  ComMask_(NULL),
  origin_(false),
  center_(false),
  ortho_(false),
  triclinic_(OFF)
{
  //fprintf(stderr,"Image Con\n");
  useMass_ = true;
} 

// DESTRUCTOR
Action_Image::~Action_Image() {
  if (ComMask_!=NULL) delete ComMask_;
}

// Action_Image::init()
/** Expected call: image [origin] [center] [triclinic | familiar [com <mask>]] <mask>  
  * - origin: center at 0.0, 0.0, 0.0, otherwise center at box center.
  * - center: Use center of mass for imaging, otherwise use first atom.
  * - triclinic: Force imaging with triclinic code.
  * - familiar: Image with triclinic code and shape into familiar trunc. oct. shape.
  * - com <mask>: If familiar, center based on COM of atoms in mask, otherwise use
  *               origin/box.
  * - <mask>: Only image atoms in <mask>. If no mask given all atoms are imaged.
  */
// Check order is:
//    1) Keywords
//    2) Masks
int Action_Image::init() {
  char *mask1;

  // Get keywords
  origin_ = actionArgs.hasKey("origin");
  center_ = actionArgs.hasKey("center");
  if (actionArgs.hasKey("familiar")) triclinic_ = FAMILIAR;
  if (actionArgs.hasKey("triclinic")) triclinic_ = FORCE;

  // Get Masks
  if (triclinic_ == FAMILIAR) {
    mask1 = actionArgs.getKeyString("com",NULL);
    if (mask1!=NULL) {
      ComMask_ = new AtomMask();
      ComMask_->SetMaskString(mask1);
    }
  }
  mask1 = actionArgs.getNextMask();
  Mask1_.SetMaskString(mask1);
  
  mprintf("    IMAGE: To");
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

  return 0;
}

// Action_Image::setup()
/** Set Imaging up for this parmtop. Get masks etc.
  * currentParm is set in Action::Setup
  */
int Action_Image::setup() {
  //atomPair apair;

  if ( currentParm->SetupCharMask( Mask1_ ) ) return 1;
  if (Mask1_.None()) {
    mprintf("Warning: Image::setup: Mask contains 0 atoms.\n");
    return 1;
  }

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

  if (triclinic_ == FAMILIAR) {
    if (ComMask_!=NULL) {
      if ( currentParm->SetupIntegerMask( *ComMask_ ) ) return 1;
      if (ComMask_->None()) {
        mprintf("Warning: Image::setup: Mask for 'familiar com' contains no atoms.\n");
        return 1;
      }
      mprintf("\tcom: mask [%s] contains %i atoms.\n",ComMask_->MaskString(),ComMask_->Nselected());
    }
  }

  // Set up atom range for each entity to be imaged. 
  // Currently imaging by molecule only, so each pair will be the first and
  // last atom of each molecule. Check that all atoms between first and last
  // are actually in the mask.
  imageList_.clear();
  imageList_.reserve( currentParm->Nmol() );
  for (Topology::mol_iterator mol = currentParm->MolStart();
                              mol != currentParm->MolEnd(); mol++)
  {
    int firstAtom = (*mol).BeginAtom();
    int lastAtom = (*mol).EndAtom();
    // Check that each atom in the range is in Mask1
    bool rangeIsValid = true;
    for (int atom = firstAtom; atom < lastAtom; atom++) {
      if (!Mask1_.AtomInCharMask(atom)) {
        rangeIsValid = false; 
        break;
      }
    }
    if (rangeIsValid) {
      imageList_.push_back( firstAtom );
      imageList_.push_back( lastAtom );
    }
  }
  mprintf("\tNumber of molecules to be imaged is %u based on mask [%s]\n", imageList_.size()/2,
           Mask1_.MaskString()); 
  // DEBUG: Print all pairs
  if (debug>0) {
    for (std::vector<int>::iterator ap = imageList_.begin();
                                    ap != imageList_.end(); ap+=2)
      mprintf("\t\tMol First-Last atom#: %i - %i\n", (*ap)+1, *(ap+1) );
  }

  return 0;  
}

// Action_Image::action()
int Action_Image::action() {

  if (ortho_)
    currentFrame->ImageOrtho(origin_, center_, useMass_, imageList_);
  else
    currentFrame->ImageNonortho(origin_, ComMask_, (triclinic_==FAMILIAR),
                                center_, useMass_, imageList_);
  return 0;
} 
