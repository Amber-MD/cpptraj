// Action_Strip
#include "Action_Strip.h"
#include "CpptrajStdio.h"
#include "ParmFile.h"

// CONSTRUCTOR
Action_Strip::Action_Strip() :
  oldParm_(NULL),
  newParm_(NULL),
  removeBoxInfo_(false)
{
  //fprintf(stderr,"Strip Con\n");
} 

// DESTRUCTOR
Action_Strip::~Action_Strip() {
  //fprintf(stderr,"Strip Des\n");
  // oldParm does not need dealloc because it is passed in
  if (newParm_!=NULL) delete newParm_;
}

// Action_Strip::init()
/** Expected call: strip <mask1> [outprefix <name>] [nobox] */
int Action_Strip::init( ) {
  // Get output stripped parm filename
  prefix_ = actionArgs.GetStringKey("outprefix");
  removeBoxInfo_ = actionArgs.hasKey("nobox");

  // Get mask of atoms to be stripped
  ArgList::ConstArg mask1 = actionArgs.getNextMask();
  //mprintf("    Mask 1: %s\n",mask1);
  if (mask1==NULL) {
    mprinterr("Error: strip: Requires atom mask.\n");
    return 1;
  }
  M1_.SetMaskString(mask1);
  // We want to strip the atoms inside the mask and keep those outside
  // the mask. Since modifyStateByMask needs to know the kept atoms,
  // invert the mask selection.
  M1_.InvertMask();

  mprintf("    STRIP: Stripping atoms in mask [%s]\n",M1_.MaskString());
  if (!prefix_.empty()) 
    mprintf("           Stripped topology will be output with prefix %s\n",prefix_.c_str());
  if (removeBoxInfo_)
    mprintf("           Any existing box information will be removed.\n");

  return 0;
}

// Action_Strip::Setup()
/** Attempt to create a new stripped down version of the input parmtop
  */
int Action_Strip::setup() {
  if (currentParm->SetupIntegerMask( M1_ )) return 1;
  //mprintf("    STRIP: Mask %s contains %i atoms\n",mask1,m1atoms);
  if (M1_.None()) {
    mprintf("Warning: strip: Mask [%s] has no atoms.\n",M1_.MaskString());
    return 1;
  }
  mprintf("\tStripping %i atoms.\n",currentParm->Natom() - M1_.Nselected());

  // Store old parm
  oldParm_ = currentParm;

  // Attempt to create new parmtop based on mask
  if (newParm_!=NULL) delete newParm_;
  newParm_ = currentParm->modifyStateByMask(M1_);
  if (newParm_==NULL) {
    mprinterr("Error: strip: Could not create new parmtop.\n");
    return 1;
  }
  // Remove box information if asked
  if (removeBoxInfo_)
    newParm_->ParmBox().SetNoBox(); 

  newParm_->Summary();

  // Allocate space for new frame
  newFrame_.SetupFrameM(newParm_->Atoms());

  // If prefix given then output stripped parm
  if (!prefix_.empty()) {
    std::string newfilename(prefix_);
    newfilename += ".";
    newfilename += oldParm_->OriginalFilename();
    mprintf("\tWriting out amber topology file %s to %s\n",newParm_->c_str(),newfilename.c_str());
    ParmFile pfile;
    pfile.SetDebug( debug );
    if ( pfile.Write( *newParm_, newfilename, ParmFile::AMBERPARM ) ) {
      mprinterr("Error: STRIP: Could not write out stripped parm file %s\n",
                newParm_->c_str());
    }
  }

  // Set parm
  currentParm = newParm_;

  return 0;  
}

// Action_Strip::action()
/** Modify the coordinate frame to reflect stripped parmtop. */
int Action_Strip::action() {

  newFrame_.SetFrame(*currentFrame, M1_);

  // Set frame
  currentFrame = &newFrame_;

  return 0;
} 

