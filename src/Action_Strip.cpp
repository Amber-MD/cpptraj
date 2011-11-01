// Strip
#include <cstdio> //sprintf 
#include <cstdlib>
#include <cstring>
#include "Action_Strip.h"
#include "CpptrajStdio.h"

Strip::Strip() {
  //fprintf(stderr,"Strip Con\n");
  // oldParm does not need dealloc because it is passed in
  oldParm=NULL;
  newParm=NULL;
  prefix=NULL;
} 

Strip::~Strip() {
  //fprintf(stderr,"Strip Des\n");
  if (newParm!=NULL) delete newParm;
}

/*
 * Strip::init()
 * Expected call: strip <mask1> [outprefix <name>]
 */
int Strip::init( ) {
  char *mask1;

  // Get mask of atoms to be stripped
  mask1 = actionArgs.getNextMask();
  //mprintf("    Mask 1: %s\n",mask1);
  if (mask1==NULL) {
    mprinterr("    Error: Strip::init: Requires atom mask.\n");
    return 1;
  }
  M1.SetMaskString(mask1);

  // Get output stripped parm filename
  prefix = actionArgs.getKeyString("outprefix",NULL);

  mprintf("    STRIP: Stripping atoms in mask [%s]\n",M1.maskString);
  if (prefix!=NULL) 
    mprintf("           Stripped topology will be output with prefix %s\n",prefix);

  return 0;
}

/*
 * Strip::Setup()
 * Attempt to create a new stripped down version of the input parmtop
 */
int Strip::setup() {
  M1.invertMask=true;   // Want to keep atoms outside the mask selection
                        // Atoms in the mask will be stripped
                        // Selected atoms will be kept.
  M1.SetupMask(P,activeReference,debug);
  //mprintf("    STRIP: Mask %s contains %i atoms\n",mask1,m1atoms);
  if (M1.None()) {
    mprinterr("      Error: Strip::setup: Mask has no atoms.\n");
    return 1;
  }
  mprintf("      STRIP: Stripping %i atoms.\n",P->natom - M1.Nselected);

  // Store old parm
  oldParm = P;

  // Attempt to create new parmtop based on mask
  if (newParm!=NULL) delete newParm;
  newParm = P->modifyStateByMask(M1.Selected, M1.Nselected);
  if (newParm==NULL) {
    mprinterr("      Error: Strip::setup: Could not create new parmtop.\n");
    return 1;
  }
  newParm->Summary();

  // Allocate space for new frame
  newFrame.SetupFrame(newParm->natom, newParm->mass);

  // If prefix given then output stripped parm
  if (prefix!=NULL && newParm->parmName==NULL) {
    newParm->parmName=(char*)malloc((strlen(oldParm->parmName)+strlen(prefix)+2)*sizeof(char));
    sprintf(newParm->parmName,"%s.%s",prefix,oldParm->parmName);
    mprintf("             Writing out amber topology file %s\n",newParm->parmName);
    if ( newParm->WriteAmberParm(newParm->parmName) ) {
      mprinterr("      Error: STRIP: Could not write out stripped parm file %s\n",
              newParm->parmName);
    }

  // Otherwise Set stripped parm name only, default prefix strip
  } else if ( newParm->parmName==NULL ) {
    newParm->parmName=(char*)malloc((strlen(oldParm->parmName)+7)*sizeof(char));
    sprintf(newParm->parmName,"strip.%s",oldParm->parmName);
  }

  // Set parm
  P = newParm;

  return 0;  
}

/*
 * Strip::action()
 * Modify the coordinate frame to reflect stripped parmtop.
 */
int Strip::action() {

  newFrame.SetFrameFromMask(F, &M1);

  // Set frame
  F = &newFrame;

  return 0;
} 


