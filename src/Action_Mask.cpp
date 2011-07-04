// ActionMask
#include <cstdlib>
#include "Action_Mask.h"
#include "CpptrajStdio.h"
#include "TrajectoryFile.h"

// CONSTRUCTOR
ActionMask::ActionMask() {
  //fprintf(stderr,"ActionMask Con\n");
  maskpdb=NULL;
} 

// DESTRUCTOR
ActionMask::~ActionMask() {
  //fprintf(stderr,"ActionMask Destructor.\n");
}

/* ActionMask::init()
 * Expected call: mask <mask1> [maskout filename] 
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int ActionMask::init( ) {
  char *mask1;
  char *maskFilename;

  // Get Keywords
  maskFilename = A->getKeyString("maskout",NULL);
  maskpdb = A->getKeyString("maskpdb",NULL);

  // Get Mask
  mask1 = A->getNextMask();
  //mprintf("    Mask 1: %s\n",mask1);
  Mask1.SetMaskString(mask1);

  mprintf("    ActionMask: Information on atoms in mask %s will be printed",
          Mask1.maskString);
  if (maskFilename!=NULL)
    mprintf(" to file %s",maskFilename);
  mprintf(".\n");
  if (maskpdb!=NULL) {
    mprintf("                PDBs of atoms in mask will be written to %s.X\n",maskpdb);
    // Set pdb output options: multi so that 1 file per frame is written; dumpq
    // so that charges are written out. 
    // NOTE: Could also split the arglist at maskpdb and make it so any type of 
    //       file can be written out.
    maskpdbarg.Add((char*)"multi\0");
    maskpdbarg.Add((char*)"dumpq\0");
  }

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

/* ActionMask::setup()
 * No setup performed here, all done in action so that the coords can be
 * passed to the mask parser.
 */
int ActionMask::setup() {
  return 0;  
}

/* ActionMask::action()
 */
int ActionMask::action() {
  int atom, res;
  TrajectoryFile pdbout;

  // Get atom selection
  if ( Mask1.SetupCharMask(P, F->X, debug) ) {
    mprintf("Warning: ActionMask::action: Could not set up atom mask [%s]\n",Mask1.maskString);
    return 1;
  }

  // Print out information for every atom in the mask
  for (atom=0; atom < P->natom; atom++) {
    if (Mask1.AtomInCharMask(atom)) {
      res = P->atomToResidue(atom);
      outfile.IO->Printf("%8i %8i %4s %8i %4s %8i",
                         currentFrame+OUTPUTFRAMESHIFT,atom+1, P->names[atom], res+1,
                         P->ResidueName(res), P->atomToMolecule(atom)+1);
      /*mprintf(" Type=%4s",P->types[atom]);
      mprintf(" Charge=%lf",P->charge[atom]);
      mprintf(" Mass=%lf",P->mass[atom]);*/
      outfile.IO->Printf("\n");
    }
  }

  // Optional PDB write out of selected atoms for the frame.
  if (maskpdb!=NULL) {
    AtomMask *Mask2 = new AtomMask();
    for (atom=0; atom < P->natom; atom++) 
      if (Mask1.AtomInCharMask(atom)) Mask2->AddAtom(atom);
    // Create new parm and frame based on atoms in Mask
    AmberParm *pdbParm = P->modifyStateByMask(Mask2->Selected, Mask2->Nselected);
    //pdbParm->Summary(); // DEBUG
    Frame *pdbFrame = new Frame(Mask2->Nselected,NULL);
    // Set only coords
    pdbFrame->SetFrameCoordsFromMask(F->X, Mask2);
    // Set up output file. Reset the arg list before every setup so that
    // any args are again available to SetupWrite
    maskpdbarg.ResetAll();
    pdbout.SetDebug(debug);
    if (pdbout.SetupWrite(maskpdb,&maskpdbarg,pdbParm,PDBFILE)) {
      mprinterr("Error: Action_Mask: maskpdb %s: Could not set up for write of frame %i.\n",
                maskpdb,currentFrame);
    } else {
      pdbout.PrintInfo(0);
      pdbout.WriteFrame(currentFrame,pdbParm,pdbFrame->X,pdbFrame->box,pdbFrame->T);
      pdbout.EndTraj();
    }
    delete pdbParm;
    delete pdbFrame;
    delete Mask2;
  }

  return 0;
} 

/* ActionMask::print()
 * Close the output file.
 */
void ActionMask::print() {
  outfile.CloseFile();
}

