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

// ActionMask::init()
/** Expected call: mask <mask1> [maskout <filename>] [maskpdb <filename>] 
  */
// NOTE: Could also split the arglist at maskpdb and make it so any type of 
//       file can be written out.
int ActionMask::init( ) {
  char *mask1;
  char *maskFilename;

  // Get Keywords
  maskFilename = actionArgs.getKeyString("maskout",NULL);
  maskpdb = actionArgs.getKeyString("maskpdb",NULL);

  // Get Mask
  mask1 = actionArgs.getNextMask();
  //mprintf("    Mask 1: %s\n",mask1);
  Mask1.SetMaskString(mask1);

  mprintf("    ActionMask: Information on atoms in mask %s will be printed",
          Mask1.MaskString());
  if (maskFilename!=NULL)
    mprintf(" to file %s",maskFilename);
  mprintf(".\n");
  if (maskpdb!=NULL) {
    mprintf("                PDBs of atoms in mask will be written to %s.X\n",maskpdb);

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

// ActionMask::setup()
/** No setup performed here, all done in action so that the coords can be
  * passed to the mask parser.
  */
int ActionMask::setup() {
  return 0;  
}

// ActionMask::action()
int ActionMask::action() {
  int atom, res;
  TrajectoryFile pdbout;

  // Get atom selection
  if ( currentParm->SetupCharMask(Mask1, currentFrame->X) ) {
    mprintf("Warning: ActionMask::action: Could not set up atom mask [%s]\n",Mask1.MaskString());
    return 1;
  }

  // Print out information for every atom in the mask
  for (atom=0; atom < currentParm->natom; atom++) {
    if (Mask1.AtomInCharMask(atom)) {
      res = currentParm->atomToResidue(atom);
      outfile.IO->Printf("%8i %8i %4s %8i %4s %8i",
                         frameNum+OUTPUTFRAMESHIFT,atom+1, currentParm->AtomName(atom), res+1,
                         currentParm->ResidueName(res), currentParm->atomToMolecule(atom)+1);
      /*mprintf(" Type=%4s",currentParm->types[atom]);
      mprintf(" Charge=%lf",currentParm->charge[atom]);
      mprintf(" Mass=%lf",currentParm->mass[atom]);*/
      outfile.IO->Printf("\n");
    }
  }

  // Optional PDB write out of selected atoms for the frame.
  if (maskpdb!=NULL) {
    AtomMask *Mask2 = new AtomMask();
    for (atom=0; atom < currentParm->natom; atom++) 
      if (Mask1.AtomInCharMask(atom)) Mask2->AddAtom(atom);
    // Create new parm and frame based on atoms in Mask
    AmberParm *pdbParm = currentParm->modifyStateByMask(Mask2->Selected, NULL);
    //pdbParm->Summary(); // DEBUG
    Frame *pdbFrame = new Frame();
    pdbFrame->SetupFrame(Mask2->Nselected,NULL);
    // Set only coords
    pdbFrame->SetFrameCoordsFromMask(currentFrame->X, Mask2);
    // Set up output file. 
    pdbout.SetDebug(debug);
    // Set pdb output options: multi so that 1 file per frame is written; dumpq
    // so that charges are written out. 
    if (pdbout.SetupWriteWithArgs(maskpdb,"multi dumpq",pdbParm,PDBFILE)) {
      mprinterr("Error: Action_Mask: maskpdb %s: Could not set up for write of frame %i.\n",
                maskpdb,frameNum);
    } else {
      pdbout.PrintInfo(0);
      pdbout.WriteFrame(frameNum,pdbParm,*pdbFrame);
      pdbout.EndTraj();
    }
    delete pdbParm;
    delete pdbFrame;
    delete Mask2;
  }

  return 0;
} 

// ActionMask::print()
/** Close the output file. */
void ActionMask::print() {
  outfile.CloseFile();
}

