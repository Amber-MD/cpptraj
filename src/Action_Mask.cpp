// ActionMask
#include <cstdlib>
#include "Action_Mask.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
ActionMask::ActionMask() {
  //fprintf(stderr,"ActionMask Con\n");
} 

// DESTRUCTOR
ActionMask::~ActionMask() {
  //fprintf(stderr,"ActionMask Destructor.\n");
}

/*
 * ActionMask::init()
 * Expected call: mask <mask1> [out filename] 
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int ActionMask::init( ) {
  char *mask1;
  char *maskFilename;

  // Get Keywords
  maskFilename = A->getKeyString("out",NULL);

  // Get Mask
  mask1 = A->getNextMask();
  //mprintf("    Mask 1: %s\n",mask1);
  Mask1.SetMaskString(mask1);

  mprintf("    ActionMask: Information on atoms in mask %s will be printed",
          Mask1.maskString);
  if (maskFilename!=NULL)
    mprintf(" to file %s",maskFilename);
  mprintf(".\n");

  // Open output file
  // NOTE: Should this be a buffer? Output at end?
  if ( outfile.SetupFile(maskFilename, WRITE, UNKNOWN_FORMAT, UNKNOWN_TYPE, debug) )
    return 1;
  if ( outfile.OpenFile() ) return 1;
  // Header
  outfile.IO->Printf("%-8s %8s %4s %8s %4s %8s","#Frame","AtomNum","Atom","ResNum","Res",
                     "MolNum");
  outfile.IO->Printf("\n");

  return 0;
}

int ActionMask::setup() {

  //if ( Mask1.SetupMask(P,debug) ) return 1;

  //if (Mask1.None()) {
  //  mprintf("    Error: ActionMask::setup: Mask has no atoms.\n");
  //  return 1;
  //}

  return 0;  
}

/*
 * ActionMask::action()
 */
int ActionMask::action() {
  int atom, res;
  char *mask;

  mask = P->mask(Mask1.maskString, F->X);
  if (mask==NULL) {
    mprintf("Warning: ActionMask::action: Atom mask %s is NULL.\n",Mask1.maskString);
    return 1;
  }
  for (atom=0; atom < P->natom; atom++) {
    if (mask[atom]=='F') continue;
    res = P->atomToResidue(atom);
    outfile.IO->Printf("%8i %8i %4s %8i %4s %8i",currentFrame,atom+1, P->names[atom], res+1,
            P->resnames[res], P->atomToMolecule(atom)+1);
    /*mprintf(" Type=%4s",P->types[atom]);
    mprintf(" Charge=%lf",P->charge[atom]);
    mprintf(" Mass=%lf",P->mass[atom]);*/
    outfile.IO->Printf("\n");
  }
  free(mask); 
  return 0;
} 

void ActionMask::print() {

  outfile.CloseFile();

}

