// Strip 
#include <cstdlib>
#include <cstring>
#include "Action_Strip.h"

Strip::Strip() {
  //fprintf(stderr,"Strip Con\n");
  // oldParm does not need dealloc because it is passed in
  oldParm=NULL;
  newParm=NULL;
  // mask1 does not need deallocation since it is a passed in Arg
  //mask1=NULL;
  //M1=NULL;
  //currentType=STRIP;
  newFrame=NULL;
  prefix=NULL;
} 

Strip::~Strip() {
  //fprintf(stderr,"Strip Des\n");
  //if (M1!=NULL) free(M1);
  if (newParm!=NULL) delete newParm;
  if (newFrame!=NULL) delete newFrame;
}

/*
 * Strip::init()
 * Expected call: strip <mask1> [outprefix <name>]
 */
int Strip::init( ) {
  char *mask1;

  // Get mask of atoms to be stripped
  mask1 = A->getNextMask();
  fprintf(stdout,"    Mask 1: %s\n",mask1);
  if (mask1==NULL) {
    fprintf(stdout,"    Error: Strip::init: Requires atom mask.\n");
    return 1;
  }
  M1.SetMaskString(mask1);

  // Get output stripped parm filename
  prefix = A->getKeyString("outprefix",NULL);

  fprintf(stdout,"    STRIP: %s",M1.maskString);
  if (prefix!=NULL) 
    fprintf(stdout,", stripped parm will be output with prefix %s",prefix);
  fprintf(stdout,"\n");

  return 0;
}

/*
 * Strip::Setup()
 * Attempt to create a new stripped down version of the input parmtop
 */
int Strip::setup() {
  fprintf(stdout,"    STRIP: ");
  M1.invertMask=true;   // Want to keep atoms outside the mask selection
                        // Atoms in the mask will be stripped
                        // Selected atoms will be kept.
  M1.SetupMask(P,debug);
  //fprintf(stdout,"    STRIP: Mask %s contains %i atoms\n",mask1,m1atoms);

  if (M1.None()) {
    fprintf(stdout,"    Error: Strip::setup: Mask has no atoms.\n");
    return 1;
  }

  // Store old parm
  oldParm = P;

  // Attempt to create new parmtop based on mask
  if (newParm!=NULL) delete newParm;
  newParm = P->modifyStateByMask(M1.Selected, M1.Nselected);
  if (newParm==NULL) {
    fprintf(stdout,"    Error: Strip::setup: Could not create new parmtop.\n");
    return 1;
  }

  fprintf(stdout,"    STRIP: %s was stripped of %i atoms.\n",oldParm->parmName,
          oldParm->natom - M1.Nselected);
  fprintf(stdout,"           New parmtop contains %i atoms.\n",newParm->natom);
  fprintf(stdout,"                                %i residues.\n",newParm->nres);
  fprintf(stdout,"                                %i molecules.\n",newParm->molecules);
  fprintf(stdout,"                                %i solvent molcules.\n",newParm->solventMolecules);

  // Allocate space for new frame
  if (newFrame!=NULL) delete newFrame;
  newFrame = new Frame(newParm->natom, NULL);

  // If prefix given then output stripped parm
  if (prefix!=NULL && newParm->parmName==NULL) {
    newParm->parmName=(char*)malloc((strlen(oldParm->parmName)+strlen(prefix)+2)*sizeof(char));
    sprintf(newParm->parmName,"%s.%s",prefix,oldParm->parmName);
    if ( newParm->WriteAmberParm() ) {
      fprintf(stdout,"ERROR: STRIP: Could not write out stripped parm file %s\n",
              newParm->parmName);
    }
    //free(outfilename);

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
  //int i;//,j,natom3;

  newFrame->SetFrameFromMask(F, &M1);

  // Set frame
  F = newFrame;

  return 0;
} 


