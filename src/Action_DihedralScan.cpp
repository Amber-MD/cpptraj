// DihedralScan
#include <cmath>
#include <cstdlib> // rand, srand
#include <ctime> // time
#include "Action_DihedralScan.h"
#include "CpptrajStdio.h"
#include "Constants.h" // DEGRAD
#include "vectormath.h" // calcRotationMatrix

// CONSTRUCTOR
DihedralScan::DihedralScan() {
  //fprintf(stderr,"DihedralScan Con\n");
  random_angle=false;
  interval = 60.0;
  outfilename=NULL;
  outfmt=NULL;
  check_for_clashes = false;
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
  check_for_clashes = actionArgs.hasKey("check");

  mprintf("    DIHEDRALSCAN: Dihedrals in mask [%s]\n",Mask1.MaskString());
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
  if (check_for_clashes)
    mprintf("                  Will attempt to recover from bad steric clashes.\n");
  
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
  int Nres;
  ResidueCheckType rct;
  std::vector<int> tmpMask;

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
        if (currentParm->MaskOfAtomsAroundBond(atom, atom2, tmpMask)) return 1;
        dst.Rmask.ResetMask();
        dst.Rmask.AddAtoms(tmpMask);
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
        mprintf("DEBUG:\t");
        BB_dihedrals[dih].Rmask.PrintMaskAtoms("Rmask:");
      }
    }
  }
  // Set up CheckStructure for this parm
  checkStructure.Setup(&currentParm);

  // CheckStructure can take quite a long time. Set up an alternative
  // structure check. First step is coarse; check distances between a 
  // certain atom in each residue (first, COM, CA, some other atom?) 
  // to see if residues are in each others neighborhood. Second step
  // is to check the atoms in each close residue.
  if (check_for_clashes) {
    if (currentParm->finalSoluteRes>0)
      Nres = currentParm->finalSoluteRes;
    else
      Nres = currentParm->nres;
    for (int res = 0; res < Nres; res++) {
      rct.resnum = res;
      rct.checkatom = currentParm->resnums[res];
      rct.start = currentParm->resnums[res];
      rct.stop = currentParm->resnums[res+1];
      ResCheck.push_back(rct);
    }
  }

  return 0;  
}

// DihedralScan::CheckResidues()
int DihedralScan::CheckResidues( Frame *frameIn, int second_atom ) {
  int Nproblems = 0;
  double atomD2;
  // Check residues in resulting structure
  std::vector<ResidueCheckType>::iterator endres1 = ResCheck.end();
  endres1--;
  for (std::vector<ResidueCheckType>::iterator res1 = ResCheck.begin();
                                               res1 != endres1;
                                               res1++)
  {
    // Check self
    for (int resatom1 = (*res1).start; resatom1 < (*res1).stop - 1; resatom1++) {
      for (int resatom2 = resatom1 + 1; resatom2 < (*res1).stop; resatom2++) {
        atomD2 = frameIn->DIST2( resatom1, resatom2 );
        if (atomD2 < 0.64) { // 0.8^2 Ang
          //mprintf("\tWarning: Atoms %i@%s and %i@%s are close (%.2lf)\n",
          mprintf("\tWarning: Res %i Atoms %i@%s and %i@%s are close\n",
                  (*res1).resnum+1, resatom1+1, currentParm->names[resatom1], 
                  resatom2+1, currentParm->names[resatom2]);//,
                  //sqrt(atomD2));
          Nproblems++;
          // If the second atom of this pair is lower than the input second
          // atom (which is the second atom of the next dihedral) need to calc
          // a new dihedral, exit now.
          if (resatom2 < second_atom) return -1;
        }    
      }
    }
    // Check other residues
    std::vector<ResidueCheckType>::iterator res2 = res1;
    res2++;
    for (; res2 != ResCheck.end(); res2++) {
      double D2 = frameIn->DIST2( (*res1).checkatom, (*res2).checkatom );
      if (D2 < 100.0) { // 10^2 Ang
        for (int resatom1 = (*res1).start; resatom1 < (*res1).stop; resatom1++) {
          for (int resatom2 = (*res2).start; resatom2 < (*res2).stop; resatom2++) {
            atomD2 = frameIn->DIST2( resatom1, resatom2 );
            if (atomD2 < 0.64) { // 0.8^2 Ang
              //mprintf("\tWarning: Atoms %i@%s and %i@%s are close (%.2lf)\n",
              mprintf("\tWarning: Atoms %i:%i@%s and %i:%i@%s are close\n",
                      (*res1).resnum+1, resatom1+1, currentParm->names[resatom1], 
                      (*res2).resnum+1, resatom2+1, currentParm->names[resatom2]);//,
                      //sqrt(atomD2));
              Nproblems++;
              // If the second atom of this pair is lower than the input second
              // atom (which is the second atom of the next dihedral) need to calc
              // a new dihedral, exit now.
              if (resatom2 < second_atom) return -1;
            }
          }
        } 
      }
    }
  } 
  return Nproblems;
}

// DihedralScan::action()
int DihedralScan::action() {
  double axisOfRotation[3], rotationMatrix[9], theta_in_degrees, theta_in_radians;
  int n_problems, second_atom;

  std::vector<DihedralScanType>::iterator next_dih = BB_dihedrals.begin();
  next_dih++;
  for (std::vector<DihedralScanType>::iterator dih = BB_dihedrals.begin();
                                               dih != BB_dihedrals.end();
                                               dih++)
  {
    // Get the second atom of the next dihdedral. If the second atom in any 
    // of the problem pairs of atoms is less than this there will be no way
    // to correct any problems, so this dihedral will need to be rotated
    // again.
    if (next_dih!=BB_dihedrals.end()) 
      second_atom = (*next_dih).atom2;
    else
      second_atom = 0;

    // Set axis of rotation
    currentFrame->SetAxisOfRotation(axisOfRotation,(*dih).atom1,(*dih).atom2);
    // Generate random value to rotate by in radians
    // Guaranteed to rotate by at least 1 degree.
    // NOTE: could potentially rotate 360 - prevent?
    theta_in_degrees = (rand() % 360) + 1;
    theta_in_radians = theta_in_degrees * DEGRAD;
    // Calculate rotation matrix for random theta
    calcRotationMatrix(rotationMatrix, axisOfRotation, theta_in_radians);
    bool rotate_dihedral = true;
    unsigned int loop_count = 0;
    while (rotate_dihedral) {
      //mprintf("DEBUG:\tRotating around %8i%4s %8i%4s, %.2lf degrees (%u).\n",(*dih).atom1,
      //        currentParm->names[(*dih).atom1],(*dih).atom2,
      //        currentParm->names[(*dih).atom2],theta_in_degrees,loop_count);
      // Rotate around axis
      currentFrame->RotateAroundAxis(rotationMatrix, theta_in_radians, (*dih).Rmask);
      // If we dont care about sterics exit here
      if (!check_for_clashes) break;
      // Check resulting structure for issues
      n_problems = CheckResidues( currentFrame, second_atom );
      if (n_problems > -1) {
        mprintf("%i\tCheckResidues: %i problems.\n",frameNum,n_problems);
        rotate_dihedral = false;
      } else if (loop_count==0) {
        mprintf("\tTrying dihedral increments of +1\n");
        // Instead of a new random dihedral, try increments of 1
        theta_in_degrees = 1.0;
        theta_in_radians = theta_in_degrees * DEGRAD;
        // Calculate rotation matrix for new theta
        calcRotationMatrix(rotationMatrix, axisOfRotation, theta_in_radians);
      }
      ++loop_count;
      if (loop_count == 360) {
        mprintf("360 iterations! Bailing!\n");
        return 1;
      }
    }
    next_dih++;
  }
  // Check the resulting structure
  n_problems = checkStructure.SeparateAction( currentFrame );
  mprintf("%i\tResulting structure has %i problems.\n",frameNum,n_problems);

  return 0;
} 

