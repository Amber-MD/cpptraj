// Image 
#include "Action_Image.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Image::Image() {
  //fprintf(stderr,"Image Con\n");
  ComMask=NULL;
  origin = false;
  center = false;
  ortho = false;
  useMass_ = true;
  triclinic = OFF;
} 

// DESTRUCTOR
Image::~Image() {
  if (ComMask!=NULL) delete ComMask;
}

// Image::init()
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
int Image::init() {
  char *mask1;

  // Get keywords
  origin = actionArgs.hasKey("origin");
  center = actionArgs.hasKey("center");
  if (actionArgs.hasKey("familiar")) triclinic = FAMILIAR;
  if (actionArgs.hasKey("triclinic")) triclinic = FORCE;

  // Get Masks
  if (triclinic == FAMILIAR) {
    mask1 = actionArgs.getKeyString("com",NULL);
    if (mask1!=NULL) {
      ComMask = new AtomMask();
      ComMask->SetMaskString(mask1);
    }
  }
  mask1 = actionArgs.getNextMask();
  Mask1.SetMaskString(mask1);
  
  mprintf("    IMAGE: To");
  if (origin)
    mprintf(" origin");
  else
    mprintf(" box center");
  mprintf(" based on");
  if (center)
    mprintf(" center of mass");
  else
    mprintf(" first atom position");
  mprintf(" using atoms in mask %s\n",Mask1.MaskString());
  if (triclinic == FORCE)
    mprintf( "           Triclinic On.\n");
  else if (triclinic == FAMILIAR) {
    mprintf( "           Triclinic On, familiar shape");
    if (ComMask!=NULL) 
      mprintf( " centering on atoms in mask %s", ComMask->MaskString());
    mprintf(".\n");
  }

  return 0;
}

// Image::setup()
/** Set Imaging up for this parmtop. Get masks etc.
  * currentParm is set in Action::Setup
  */
int Image::setup() {
  //atomPair apair;

  if ( currentParm->SetupCharMask( Mask1 ) ) return 1;
  if (Mask1.None()) {
    mprintf("Warning: Image::setup: Mask contains 0 atoms.\n");
    return 1;
  }

  if (currentParm->BoxType()==Box::NOBOX) {
    mprintf("Warning: Image::setup: Parm %s does not contain box information.\n",
            currentParm->c_str());
    return 1;
  }

  ortho = false;  
  if (currentParm->BoxType()==Box::ORTHO && triclinic==OFF) ortho=true;

  // If box is originally truncated oct and not forcing triclinic, 
  // turn familiar on.
  if (currentParm->BoxType()==Box::TRUNCOCT && triclinic!=FORCE && triclinic!=FAMILIAR) {
    mprintf("\tOriginal box is truncated octahedron, turning on 'familiar'.\n");
    triclinic=FAMILIAR;
  }

  if (triclinic == FAMILIAR) {
    if (ComMask!=NULL) {
      if ( currentParm->SetupIntegerMask( *ComMask ) ) return 1;
      if (ComMask->None()) {
        mprintf("Warning: Image::setup: Mask for 'familiar com' contains no atoms.\n");
        return 1;
      }
      mprintf("\tcom: mask [%s] contains %i atoms.\n",ComMask->MaskString(),ComMask->Nselected());
    }
  }

  // Set up atom range for each entity to be imaged. 
  // Currently imaging by molecule only, so each pair will be the first and
  // last atom of each molecule. Check that all atoms between first and last
  // are actually in the mask.
  imageList.clear();
  imageList.reserve( currentParm->Nmol() );
  for (Topology::mol_iterator mol = currentParm->MolStart();
                              mol != currentParm->MolEnd(); mol++)
  {
    int firstAtom = (*mol).BeginAtom();
    int lastAtom = (*mol).EndAtom();
    // Check that each atom in the range is in Mask1
    bool rangeIsValid = true;
    for (int atom = firstAtom; atom < lastAtom; atom++) {
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
  mprintf("\tNumber of molecules to be imaged is %u based on mask [%s]\n", imageList.size()/2,
           Mask1.MaskString()); 
  // DEBUG: Print all pairs
  if (debug>0) {
    for (std::vector<int>::iterator ap = imageList.begin();
                                    ap != imageList.end(); ap+=2)
      mprintf("\t\tMol First-Last atom#: %i - %i\n", (*ap)+1, *(ap+1) );
  }

  return 0;  
}

// Image::action()
int Image::action() {

  if (ortho)
    currentFrame->ImageOrtho(origin, center, useMass_, imageList);
  else
    currentFrame->ImageNonortho(origin, ComMask, (triclinic==FAMILIAR),
                                center, useMass_, imageList);

  return 0;
} 


