// Strip
#include <cstdio> //sprintf 
#include <cstring>
#include "Action_Strip.h"
#include "CpptrajStdio.h"
#include "ParmFile.h"

// CONSTRUCTOR
Strip::Strip() {
  //fprintf(stderr,"Strip Con\n");
  // oldParm does not need dealloc because it is passed in
  oldParm=NULL;
  newParm=NULL;
  prefix=NULL;
  removeBoxInfo=false;
} 

// DESTRUCTOR
Strip::~Strip() {
  //fprintf(stderr,"Strip Des\n");
  if (newParm!=NULL) delete newParm;
}

// Strip::init()
/** Expected call: strip <mask1> [outprefix <name>] [nobox] */
int Strip::init( ) {
  char *mask1;

  // Get output stripped parm filename
  prefix = actionArgs.getKeyString("outprefix",NULL);
  removeBoxInfo = actionArgs.hasKey("nobox");

  // Get mask of atoms to be stripped
  mask1 = actionArgs.getNextMask();
  //mprintf("    Mask 1: %s\n",mask1);
  if (mask1==NULL) {
    mprinterr("Error: Strip::init: Requires atom mask.\n");
    return 1;
  }
  M1.SetMaskString(mask1);
  // We want to strip the atoms inside the mask and keep those outside
  // the mask. Since modifyStateByMask needs to know the kept atoms,
  // invert the mask selection.
  M1.InvertMask();

  mprintf("    STRIP: Stripping atoms in mask [%s]\n",M1.MaskString());
  if (prefix!=NULL) 
    mprintf("           Stripped topology will be output with prefix %s\n",prefix);
  if (removeBoxInfo)
    mprintf("           Any existing box information will be removed.\n");

  return 0;
}

// Strip::Setup()
/** Attempt to create a new stripped down version of the input parmtop
  */
int Strip::setup() {
  if (currentParm->SetupIntegerMask( M1 )) return 1;
  //mprintf("    STRIP: Mask %s contains %i atoms\n",mask1,m1atoms);
  if (M1.None()) {
    mprintf("Warning: Strip::setup: Mask [%s] has no atoms.\n",M1.MaskString());
    return 1;
  }
  mprintf("\tStripping %i atoms.\n",currentParm->Natom() - M1.Nselected());

  // Store old parm
  oldParm = currentParm;

  // Attempt to create new parmtop based on mask
  if (newParm!=NULL) delete newParm;
  newParm = currentParm->modifyStateByMask(M1, prefix);
  if (newParm==NULL) {
    mprinterr("Error: Strip::setup: Could not create new parmtop.\n");
    return 1;
  }
  // Remove box information if asked
  if (removeBoxInfo)
    newParm->SetNoBox(); 

  newParm->Summary();

  // Allocate space for new frame
  newFrame.SetupFrame(newParm->Natom(), newParm->Mass());

  // If prefix given then output stripped parm
  if (prefix!=NULL) {
    mprintf("\tWriting out amber topology file %s\n",newParm->c_str());
    ParmFile pfile;
    pfile.SetDebug( debug );
    if ( pfile.Write( *newParm, newParm->ParmName(), ParmFile::AMBERPARM ) ) {
    //if ( newParm->WriteTopology(newParm->parmName) ) {
      mprinterr("Error: STRIP: Could not write out stripped parm file %s\n",
                newParm->c_str());
    }
  }

  // Set parm
  currentParm = newParm;

  return 0;  
}

// Strip::action()
/** Modify the coordinate frame to reflect stripped parmtop. */
int Strip::action() {

  newFrame.SetFrame(*currentFrame,M1);

  // Set frame
  currentFrame = &newFrame;

  return 0;
} 

