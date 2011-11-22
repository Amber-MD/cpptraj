// DihedralScan
#include <cmath>
#include <cstdlib> // rand, srand
#include <ctime> // time
#include "Action_DihedralScan.h"
#include "CpptrajStdio.h"
#include "Constants.h" // DEGRAD

// CONSTRUCTOR
DihedralScan::DihedralScan() {
  //fprintf(stderr,"DihedralScan Con\n");
  random_angle=false;
  interval = 60.0;
  outfilename=NULL;
  outfmt=NULL;
} 

// DESTRUCTOR
DihedralScan::~DihedralScan() {
  //fprintf(stderr,"DihedralScan Destructor.\n");
}

// DihedralScan::init()
/// Expected call: dihedralscan <mask> [<interval> | random ] outtraj <traj>
// Dataset name will be the last arg checked for. Check order is:
//    1) Keywords
//    2) Masks
//    3) Dataset name
int DihedralScan::init( ) {
  int rseed_arg;
  unsigned int rseed;
  // Get mask
  Mask1.SetMaskString( actionArgs.getNextMask() );

  // Get Keywords
  random_angle = actionArgs.hasKey("random");
  if (!random_angle) 
    interval = actionArgs.getNextDouble(60.0);
  outfilename = actionArgs.getKeyString("outtraj",NULL);
  outfmt = actionArgs.getKeyString("outfmt",NULL);
  rseed_arg = actionArgs.getKeyInt("rseed",-1);

  mprintf("    DIHEDRALSCAN: Dihedrals in mask [%s]\n",Mask1.maskString);
  if (random_angle) {
    mprintf("                  Dihedrals will be rotated to random values.\n");
    if (rseed_arg==-1)
      mprintf("                  Random number generator will be seeded using time.\n");
    else
      mprintf("                  Random number generator will be seeded using %i\n",rseed_arg);
  } else
    mprintf("                  Dihedrals will be rotated at intervals of %.2lf degrees.\n",
            interval);
  if (outfilename!=NULL)
    mprintf("                  Coordinates output to %s, format %s\n",outfilename, outfmt);
  
  // Seed random number generator
  if (random_angle) {
    if (rseed_arg == -1)
      rseed = (unsigned int) time(NULL);
    else
      rseed = (unsigned int) rseed_arg;
    mprintf("                  Seeding random number generator: [%u]\n",rseed);
    srand( rseed );
  }

  // Initialize CheckStructure
  checkStructure.SeparateInit(-1,-1,debug);

  return 0;
}

// DihedralScan::setup()
/** Determine from selected mask atoms which dihedrals will be rotated.
  */
int DihedralScan::setup() {
  DihedralScanType dst;
  BondInfo mol;

  if ( Mask1.SetupCharMask(currentParm,activeReference,debug) ) return 1;
  if (Mask1.None()) {
    mprintf("    Error: DihedralScan::setup: Mask has no atoms.\n");
    return 1;
  }

  // Need bonding info to determine whether two atoms are actually
  // bonded, so that the bond can be rotated around. Set up for all
  // atoms. 
  if (currentParm->SetupBondInfo()) return 1;

  // For now just focus on backbone phi/psi dihedrals:
  //   C-N-CA-C  N-CA-C-N
  for (int atom=0; atom < currentParm->natom; atom++) {
    if (Mask1.AtomInCharMask(atom)) {
      int atom2 = -1;
      // PHI: C-N-CA-C
      if (currentParm->AtomNameIs(atom,"N   ")) {
        atom2 = currentParm->GetBondedAtomIdx(atom,"CA  ");
      // PSI: N-CA-C-N
      } else if (currentParm->AtomNameIs(atom,"CA  ")) {
        atom2 = currentParm->GetBondedAtomIdx(atom,"C   ");
      }
      // If a second atom was found dihedral is defined and in mask, store it
      if (atom2 != -1 && Mask1.AtomInCharMask(atom2)) {
        // Set up mask of atoms that will move upon rotation of dihedral
        int tmpNmask = 0;
        int *tmpMask = currentParm->MaskOfAtomsAroundBond(atom, atom2, &tmpNmask);
        dst.Rmask.ResetMask();
        dst.Rmask.AddAtoms(tmpMask,tmpNmask);
        if (tmpMask!=NULL) delete[] tmpMask;
        dst.atom1 = atom;
        dst.atom2 = atom2;
        dst.currentVal = 0;
        dst.interval = interval;
        dst.maxVal = (int) (360.0 / interval);
        dst.isRandom = random_angle;
        BB_dihedrals.push_back(dst);
      }
    } // END if atom in char mask
  }

  // DEBUG: List defined dihedrals
  if (debug>0) {
    mprintf("DEBUG: Dihedrals (central 2 atoms only):\n");
    for (unsigned int dih = 0; dih < BB_dihedrals.size(); dih++) {
      mprintf("\t%8i%4s %8i%4s\n",BB_dihedrals[dih].atom1, 
              currentParm->names[BB_dihedrals[dih].atom1], BB_dihedrals[dih].atom2,
              currentParm->names[BB_dihedrals[dih].atom2]);
      if (debug>1) {
        mprintf("\tRmask:");
        BB_dihedrals[dih].Rmask.PrintMaskAtoms();
        mprintf("\n");
      }
    }
  }
  // Set up CheckStructure for this parm
  checkStructure.Setup(&currentParm); 

  return 0;  
}

// DihedralScan::action()
int DihedralScan::action() {
  double axisOfRotation[3], theta_in_degrees, theta_in_radians;

  for (std::vector<DihedralScanType>::iterator dih = BB_dihedrals.begin();
                                               dih != BB_dihedrals.end();
                                               dih++)
  {
    // Set axis of rotation
    currentFrame->SetAxisOfRotation(axisOfRotation,(*dih).atom1,(*dih).atom2);
    // Generate random value to rotate by in radians
    // Guaranteed to rotate by at least 1 degree.
    // NOTE: could potentially rotate 360 - prevent?
    theta_in_degrees = (rand() % 360) + 1;
    theta_in_radians = theta_in_degrees * DEGRAD;
    //mprintf("DEBUG:\tRotating around %8i%4s %8i%4s, %.2lf degrees.\n",(*dih).atom1,
    //        currentParm->names[(*dih).atom1],(*dih).atom2,
    //        currentParm->names[(*dih).atom2],theta_in_degrees);
    // Rotate around axis
    currentFrame->RotateAroundAxis(axisOfRotation, theta_in_radians, (*dih).Rmask);
  }
  // Check the resulting structure
  int n_problems = checkStructure.SeparateAction( currentFrame );
  mprintf("%i\tResulting structure has %i problems.\n",frameNum,n_problems);
 
  return 0;
} 

