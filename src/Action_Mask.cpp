// Action_Mask
#include "Action_Mask.h"
#include "CpptrajStdio.h"
#include "Trajout.h"

// CONSTRUCTOR
Action_Mask::Action_Mask() 
{ 
  //fprintf(stderr,"ActionMask Con\n");
} 

// Action_Mask::init()
/** Expected call: mask <mask1> [maskout <filename>] [maskpdb <filename>] 
  */
// NOTE: Could also split the arglist at maskpdb and make it so any type of 
//       file can be written out.
int Action_Mask::init( ) {
  // Get Keywords
  std::string maskFilename = actionArgs.GetStringKey("maskout");
  maskpdb_ = actionArgs.GetStringKey("maskpdb");

  // Get Mask
  Mask1_.SetMaskString( actionArgs.getNextMask() );

  mprintf("    ActionMask: Information on atoms in mask %s will be printed",
          Mask1_.MaskString());
  if (!maskFilename.empty())
    mprintf(" to file %s",maskFilename.c_str());
  mprintf(".\n");
  if (!maskpdb_.empty()) 
    mprintf("\tPDBs of atoms in mask will be written to %s.X\n",maskpdb_.c_str());

  // Open output file
  // TODO: Buffer write out
  if ( outfile_.OpenWrite( maskFilename ) )
    return 1;
  // Header
  outfile_.Printf("%-8s %8s %4s %8s %4s %8s\n","#Frame","AtomNum","Atom",
                  "ResNum","Res", "MolNum");

  return 0;
}

// Action_Mask::action()
int Action_Mask::action() {
  Trajout pdbout;

  // Get atom selection
  if ( currentParm->SetupCharMask(Mask1_, *currentFrame) ) {
    mprintf("Warning: ActionMask::action: Could not set up atom mask [%s]\n",
            Mask1_.MaskString());
    return 1;
  }

  // Print out information for every atom in the mask
  for (int atom=0; atom < currentParm->Natom(); atom++) {
    if (Mask1_.AtomInCharMask(atom)) {
      int res = (*currentParm)[atom].ResNum();
      outfile_.Printf("%8i %8i %4s %8i %4s %8i\n", frameNum+OUTPUTFRAMESHIFT,
                      atom+1, (*currentParm)[atom].c_str(), res+1,
                      currentParm->ResidueName(res), (*currentParm)[atom].Mol()+1);
      /*mprintf(" Type=%4s",currentParm->types[atom]);
      mprintf(" Charge=%lf",currentParm->charge[atom]);
      mprintf(" Mass=%lf",currentParm->mass[atom]);
      outfile.Printf("\n");*/
    }
  }

  // Optional PDB write out of selected atoms for the frame.
  if (!maskpdb_.empty()) {
    // Convert Mask1 to an integer mask for use in parm/frame functions
    AtomMask Mask2 = Mask1_;
    Mask2.ConvertMaskType();
    // Create new parm and frame based on atoms in Mask
    Topology* pdbParm = currentParm->modifyStateByMask(Mask2);
    //pdbParm->Summary(); // DEBUG
    Frame pdbFrame(*currentFrame, Mask2);
    // Set up output file. 
    pdbout.SetDebug(debug);
    // Set pdb output options: multi so that 1 file per frame is written; dumpq
    // so that charges are written out. 
    if (pdbout.SetupTrajWriteWithArgs(maskpdb_,"multi dumpq",pdbParm,TrajectoryFile::PDBFILE)) 
    {
      mprinterr("Error: Action_Mask: maskpdb %s: Could not set up for write of frame %i.\n",
                maskpdb_.c_str(),frameNum);
    } else {
      pdbout.PrintInfo(0);
      pdbout.WriteFrame(frameNum,pdbParm,pdbFrame);
      pdbout.EndTraj();
    }
    delete pdbParm;
  }

  return 0;
} 

// Action_Mask::print()
/** Close the output file. */
void Action_Mask::print() {
  outfile_.CloseFile();
}

