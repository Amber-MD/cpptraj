// Action_DihedralScan
#include <cmath>
#include <ctime> // time
#include "Action_DihedralScan.h"
#include "CpptrajStdio.h"
#include "Constants.h" // DEGRAD
#include "vectormath.h" // calcRotationMatrix
#include "DistRoutines.h"

// Activate DEBUG info
//#define DEBUG_DIHEDRALSCAN
#ifdef DEBUG_DIHEDRALSCAN
#include "TrajectoryFile.h"
#endif

// CONSTRUCTOR
Action_DihedralScan::Action_DihedralScan() {
  //fprintf(stderr,"DihedralScan Con\n");
  random_angle=false;
  interval = 60.0;
  check_for_clashes = false;
  max_rotations = 0;
  cutoff = 0.64; // 0.8^2
  rescutoff = 100.0; // 10.0^2
  backtrack = 5;
  increment = 1;
  max_increment = 360;
  max_factor = 2;
  debug_ = 0;
  CurrentParm_ = 0;
} 

void Action_DihedralScan::Help() {
  mprintf("dihedralscan <mask> [<interval> | random ] [rseed <rseed>] outtraj <traj>\n");
  mprintf("             [ check [cutoff <cutoff>] [rescutoff <rescutoff>]\n");
  mprintf("             [backtrack <backtrack>] [increment <increment>] [maxfactor <max_factor>] ]\n");
}

// Action_DihedralScan::init()
Action::RetType Action_DihedralScan::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  // Get mask
  Mask1.SetMaskString( actionArgs.GetMaskNext() );

  // Get Keywords
  random_angle = actionArgs.hasKey("random");
  if (!random_angle) 
    interval = actionArgs.getNextDouble(60.0);
  outfilename_ = actionArgs.GetStringKey("outtraj");
  outfmt_ = actionArgs.GetStringKey("outfmt");
  check_for_clashes = actionArgs.hasKey("check");
  cutoff = actionArgs.getKeyDouble("cutoff",0.8);
  rescutoff = actionArgs.getKeyDouble("rescutoff",10.0);
  backtrack = actionArgs.getKeyInt("backtrack",4);
  increment = actionArgs.getKeyInt("increment",1);
  max_factor = actionArgs.getKeyInt("maxfactor",2);
  std::string problemFile = actionArgs.GetStringKey("out");
  // Seed random number gen
  int rseed_arg = actionArgs.getKeyInt("rseed",-1);
  RN_.rn_set( rseed_arg );

  // Check validity of args
  if (cutoff < SMALL) {
    mprinterr("Error: cutoff too small.\n");
    return Action::ERR;
  }
  if (rescutoff < SMALL) {
    mprinterr("Error: rescutoff too small.\n");
    return Action::ERR;
  }
  if (backtrack < 0) {
    mprinterr("Error: backtrack value must be >= 0\n");
    return Action::ERR;
  }
  if ( increment<1 || (360 % increment)!=0 ) {
    mprinterr("Error: increment must be a factor of 360.\n");
    return Action::ERR;
  }

  // Dataset to store number of problems
  number_of_problems = DSL->AddSet(DataSet::INT, actionArgs.GetStringNext(),"Nprob");
  if (number_of_problems==NULL) return Action::ERR;
  // Add dataset to data file list
  DFL->AddSetToFile(problemFile,number_of_problems);

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
  if (!outfilename_.empty())
    mprintf("\tCoordinates output to %s, format %s\n",outfilename_.c_str(), outfmt_.c_str());
  if (check_for_clashes) {
    mprintf("\tWill attempt to recover from bad steric clashes.\n");
    mprintf("\tAtom cutoff %.2lf, residue cutoff %.2lf, backtrack = %i\n",
            cutoff, rescutoff, backtrack);
    mprintf("\tWhen clashes occur dihedral will be incremented by %i\n",increment);
    mprintf("\tMax # attempted rotations = %i times number dihedrals.\n",
            max_factor);
  }
  
  // Square cutoffs to compare to dist^2 instead of dist
  cutoff *= cutoff;
  rescutoff *= rescutoff;

  // Increment backtrack by 1 since we need to skip over current res
  ++backtrack;

  // Calculate max increment
  max_increment = 360 / increment;

  // Initialize CheckStructure
  checkStructure.SeparateInit(-1,-1,debug_);

  return Action::OK;
}

static int GetBondedAtomIdx(Topology const& topIn, int atom, NameType nameIn) {
  for (Atom::bond_iterator bndatm = topIn[atom].bondbegin();
                           bndatm != topIn[atom].bondend(); ++bndatm)
  {
    if ( topIn[*bndatm].Name() == nameIn ) return *bndatm;
  }
  return -1;
}

static void VisitAtom( Topology const& topIn, int atm, std::vector<char>& Visited )
{
  // If this atom has already been visited return
  if (Visited[atm] == 'T') return;
  // Mark this atom as visited
  Visited[atm] = 'T';
  // Visit each atom bonded to this atom
  for (Atom::bond_iterator bondedatom = topIn[atm].bondbegin();
                           bondedatom != topIn[atm].bondend(); ++bondedatom)
    VisitAtom(topIn, *bondedatom, Visited);
}

// Action_DihedralScan::setup()
/** Determine from selected mask atoms which dihedrals will be rotated.
  */
Action::RetType Action_DihedralScan::Setup(Topology* currentParm, Topology** parmAddress) {
  DihedralScanType dst;

  if ( currentParm->SetupCharMask( Mask1 ) ) return Action::ERR;
  if (Mask1.None()) {
    mprintf("    Error: DihedralScan::setup: Mask has no atoms.\n");
    return Action::ERR;
  }

  // For now just focus on backbone phi/psi dihedrals:
  //   C-N-CA-C  N-CA-C-N
  std::vector<char> Visited( currentParm->Natom(), 'F' );
  for (int atom=0; atom < currentParm->Natom(); atom++) {
    if (Mask1.AtomInCharMask(atom)) {
      int atom2 = -1;
      // PHI: C-N-CA-C
      if ((*currentParm)[atom].Name() == "N   ") {
        atom2 = GetBondedAtomIdx(*currentParm, atom, "CA  ");
      // PSI: N-CA-C-N
      } else if ((*currentParm)[atom].Name() == "CA  ") {
        atom2 = GetBondedAtomIdx(*currentParm, atom, "C   ");
      }
      // If a second atom was found dihedral is defined and in mask, store it
      if (atom2 != -1 && Mask1.AtomInCharMask(atom2)) {
        // Set up mask of atoms that will move upon rotation of dihedral.
        // Also set up mask of atoms in this residue that will not move
        // upon rotation of dihedral, including atom2
        //if (currentParm->MaskOfAtomsAroundBond(atom, atom2, tmpMask)) return Action::ERR;
        dst.Rmask.ResetMask();
        Visited.assign( currentParm->Natom(), 'F' );
        for (Atom::bond_iterator bndatm = (*currentParm)[atom2].bondbegin();
                                 bndatm != (*currentParm)[atom2].bondend(); ++bndatm)
        {
          if ( *bndatm != atom )
            VisitAtom( *currentParm, *bndatm, Visited );
        }
        dst.checkAtoms.clear();
        int a1res = (*currentParm)[atom].ResNum();
        int a1res_start = currentParm->ResFirstAtom( a1res );
        int a1res_stop = currentParm->ResLastAtom( a1res );
        for (int maskatom = 0; maskatom < (int)Visited.size(); maskatom++) {
          if (Visited[maskatom]=='T')
            dst.Rmask.AddAtom(maskatom);
          else {
            // If this atom is in the same residue but will not move, it needs
            // to be checked for clashes since further rotations will not
            // help it.
            if (maskatom >= a1res_start && maskatom < a1res_stop   )
              dst.checkAtoms.push_back( maskatom );
          }
        }
        dst.checkAtoms.push_back( atom2 ); // atom already included from tmpMask
        dst.atom1 = atom;
        dst.atom2 = atom2;
        // Since only the second atom and atoms it is
        // bonded to move during rotation, base the check
        // on the residue of the second atom.
        dst.resnum = (*currentParm)[atom2].ResNum();
        dst.currentVal = 0;
        dst.interval = interval;
        dst.maxVal = (int) (360.0 / interval);
        dst.isRandom = random_angle;
        BB_dihedrals.push_back(dst);
      }
    } // END if atom in char mask
  }

  // DEBUG: List defined dihedrals
  if (debug_>0) {
    mprintf("DEBUG: Dihedrals (central 2 atoms only):\n");
    for (unsigned int dih = 0; dih < BB_dihedrals.size(); dih++) {
      mprintf("\t%8i%4s %8i%4s %8i%4s\n",
              BB_dihedrals[dih].atom1, (*currentParm)[BB_dihedrals[dih].atom1].c_str(), 
              BB_dihedrals[dih].atom2, (*currentParm)[BB_dihedrals[dih].atom2].c_str(), 
              BB_dihedrals[dih].resnum, currentParm->Res(BB_dihedrals[dih].resnum).c_str());
      if (debug_>1) {
        mprintf("\t\tCheckAtoms=");
        for (std::vector<int>::iterator ca = BB_dihedrals[dih].checkAtoms.begin();
                                        ca != BB_dihedrals[dih].checkAtoms.end(); ca++)
        {
          mprintf(" %i",*ca);
        }
        mprintf("\n");
      }
      if (debug_>2) {
        mprintf("\t\t");
        BB_dihedrals[dih].Rmask.PrintMaskAtoms("Rmask:");
      }
    }
  }
  // Set up CheckStructure for this parm
  checkStructure.Setup(currentParm, &currentParm);

  // Set the overall max number of rotations to try
  max_rotations = (int) BB_dihedrals.size();
  max_rotations *= max_factor;

  // CheckStructure can take quite a long time. Set up an alternative
  // structure check. First step is coarse; check distances between a 
  // certain atom in each residue (first, COM, CA, some other atom?) 
  // to see if residues are in each others neighborhood. Second step
  // is to check the atoms in each close residue.
  if (check_for_clashes) {
    ResidueCheckType rct;
    int Nres = currentParm->FinalSoluteRes();
    for (int res = 0; res < Nres; res++) {
      rct.resnum = res;
      rct.start = currentParm->ResFirstAtom( res );
      rct.stop = currentParm->ResLastAtom( res );
      rct.checkatom = rct.start;
      ResCheck.push_back(rct);
    }
  }
  CurrentParm_ = currentParm;
  return Action::OK;  
}

/*
// Action_DihedralScan::CheckResidues()
int Action_DihedralScan::CheckResidues( Frame *frameIn, int second_atom ) {
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
*/

// Action_DihedralScan::CheckResidue()
/** \return 1 if a new dihedral should be tried, 0 if no clashes, -1 if
  * \return further rotations will not help.
  */
int Action_DihedralScan::CheckResidue( Frame *FrameIn, DihedralScanType &dih, int nextres,
                                double *clash ) 
{
  int resnumIn = dih.resnum;
  int rstart = ResCheck[ resnumIn ].start;
  int rstop = ResCheck[ resnumIn ].stop;
  int rcheck = ResCheck[ resnumIn ].checkatom;
  // Check for clashes with self
#ifdef DEBUG_DIHEDRALSCAN
  mprintf("\tChecking residue %i\n",resnumIn+1);
  mprintf("\tATOMS %i to %i\n",rstart+1,rstop);
#endif
  for (int atom1 = rstart; atom1 < rstop - 1; atom1++) {
    for (int atom2 = atom1 + 1; atom2 < rstop; atom2++) {
      double atomD2 = DIST2_NoImage(FrameIn->XYZ(atom1), FrameIn->XYZ(atom2));
      if (atomD2 < cutoff) {
#ifdef DEBUG_DIHEDRALSCAN 
        mprintf("\t\tRes %i Atoms %i@%s and %i@%s are close (%.3lf)\n", resnumIn+1, 
                atom1+1, currentParm->AtomName(atom1),
                atom2+1, currentParm->AtomName(atom2), sqrt(atomD2));
#endif
        *clash = atomD2;
        return 1;
      }
    }
  }
  // Check for clashes with previous residues, as well as clashes up to and
  // including the next residue in which a dihedral will be rotated.
  for (int res = 0; res <= nextres; res++) {
    if (res == resnumIn) continue;
    int rstart2 = ResCheck[ res ].start;
    int rstop2 = ResCheck[ res ].stop;
    int rcheck2 = ResCheck[ res ].checkatom;
    double resD2 = DIST2_NoImage(FrameIn->XYZ(rcheck), FrameIn->XYZ(rcheck2));
    // If residues are close enough check each atom
    if (resD2 < rescutoff) { 
#ifdef DEBUG_DIHEDRALSCAN
      mprintf("\tRES %i ATOMS %i to %i\n",res+1,rstart2+2,rstop2);
#endif
      for (int atom1 = rstart; atom1 < rstop; atom1++) {
        for (int atom2 = rstart2; atom2 < rstop2; atom2++) {
          double D2 = DIST2_NoImage(FrameIn->XYZ(atom1), FrameIn->XYZ(atom2));
          if (D2 < cutoff) {
#ifdef DEBUG_DIHEDRALSCAN
            mprintf("\t\tRes %i atom %i@%s and res %i atom %i@%s are close (%.3lf)\n", resnumIn+1,
                    atom1+1, currentParm->AtomName(atom1), res+1,
                    atom2+1, currentParm->AtomName(atom2), sqrt(D2));
#endif
            *clash = D2;
            // If the clash involves any atom that will not be moved by further
            // rotation, indicate it is not possible to resolve clash by
            // more rotation by returning -1.
            //if (atom1 == dih.atom2 || atom1 == dih.atom1) return -1;
            for (std::vector<int>::iterator ca = dih.checkAtoms.begin();
                                            ca != dih.checkAtoms.end(); ca++) 
            {
              if (atom1 == *ca) return -1;
            }
            return 1;
          }
        }
      }
    }
  }
  return 0;
} 

// Action_DihedralScan::action()
Action::RetType Action_DihedralScan::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  double axisOfRotation[3], rotationMatrix[9], theta_in_degrees, theta_in_radians;
  int n_problems, next_resnum;
  double bestClash, clash; // The largest distance observed during clashes
  int bestLoop = 0;
  int number_of_rotations=0;
#ifdef DEBUG_DIHEDRALSCAN
  // DEBUG
  int debugframenum=0;
  TrajectoryFile DebugTraj;
  DebugTraj.SetupWrite((char*)"debugtraj.nc\0",NULL,currentParm,AMBERNETCDF);
  DebugTraj.WriteFrame(debugframenum++,currentParm,*currentFrame);
#endif

  std::vector<DihedralScanType>::iterator next_dih = BB_dihedrals.begin();
  next_dih++;
  for (std::vector<DihedralScanType>::iterator dih = BB_dihedrals.begin();
                                               dih != BB_dihedrals.end();
                                               dih++)
  {
    ++number_of_rotations;
    // Get the residue atom of the next dihedral. Residues up to and
    // including this residue will be checked for bad clashes 
    if (next_dih!=BB_dihedrals.end()) 
      next_resnum = (*next_dih).resnum;
    else
      next_resnum = (*dih).resnum-1;

    // Set axis of rotation
    currentFrame->SetAxisOfRotation(axisOfRotation,(*dih).atom1,(*dih).atom2);
    // Generate random value to rotate by in radians
    // Guaranteed to rotate by at least 1 degree.
    // NOTE: could potentially rotate 360 - prevent?
    // FIXME: Just use 2PI and rn_gen, get everything in radians
    if ((*dih).isRandom)
      theta_in_degrees = ((int)(RN_.rn_gen()*100000) % 360) + 1;
    else
      theta_in_degrees = (*dih).interval;
    theta_in_radians = theta_in_degrees * DEGRAD;
    // Calculate rotation matrix for random theta
    calcRotationMatrix(rotationMatrix, axisOfRotation, theta_in_radians);
    bool rotate_dihedral = (*dih).isRandom;
    int loop_count = 0;
    clash = 0;
    bestClash = 0;
    if (debug_>0) mprintf("DEBUG: Rotating res %8i:\n",(*dih).resnum+1);
    while (rotate_dihedral) {
      if (debug_>0) {
        mprintf("\t%8i %8i%4s %8i%4s, +%.2lf degrees (%i).\n",(*dih).resnum+1,
                (*dih).atom1+1, (*CurrentParm_)[(*dih).atom1].c_str(),
                (*dih).atom2+1, (*CurrentParm_)[(*dih).atom2].c_str(),
                theta_in_degrees,loop_count);
      }
      // Rotate around axis
      currentFrame->Rotate(rotationMatrix, (*dih).Rmask);
#ifdef DEBUG_DIHEDRALSCAN
      // DEBUG
      DebugTraj.WriteFrame(debugframenum++,currentParm,*currentFrame);
#endif
      // If we dont care about sterics exit here
      if (!check_for_clashes) break;
      // Check resulting structure for issues
      int checkresidue = CheckResidue(currentFrame, *dih, next_resnum, &clash);
      if (checkresidue==0)
        rotate_dihedral = false;
      else if (checkresidue==-1) {
        dih--; //  0
        dih--; // -1
        next_dih = dih;
        next_dih++;
        if (debug_>0)
          mprintf("\tCannot resolve clash with further rotations, trying previous again.\n");
        break;
      }
      if (clash > bestClash) {bestClash = clash; bestLoop = loop_count;}
      //n_problems = CheckResidues( currentFrame, second_atom );
      //if (n_problems > -1) {
      //  mprintf("%i\tCheckResidues: %i problems.\n",frameNum,n_problems);
      //  rotate_dihedral = false;
      //} else if (loop_count==0) {
      if (loop_count==0 && rotate_dihedral) {
        if (debug_>0)
          mprintf("\tTrying dihedral increments of +%i\n",increment);
        // Instead of a new random dihedral, try increments
        theta_in_degrees = (double)increment;
        theta_in_radians = theta_in_degrees * DEGRAD;
        // Calculate rotation matrix for new theta
        calcRotationMatrix(rotationMatrix, axisOfRotation, theta_in_radians);
      }
      ++loop_count;
      if (loop_count == max_increment) {
        if (debug_>0)
          mprintf("%i iterations! Best clash= %.3lf at %i\n",max_increment,
                  sqrt(bestClash),bestLoop);
        for (int bt = 0; bt < backtrack; bt++)
          dih--;
        next_dih = dih;
        next_dih++;
        if (debug_>0)
          mprintf("\tCannot resolve clash with further rotations, trying previous %i again.\n",
                  backtrack - 1);
        break;
        // Calculate how much to rotate back in order to get to best clash
        /*int num_back = bestLoop - 359;
        theta_in_degrees = (double) num_back;
        theta_in_radians = theta_in_degrees * DEGRAD;
        // Calculate rotation matrix for theta
        calcRotationMatrix(rotationMatrix, axisOfRotation, theta_in_radians);
        // Rotate back to best clash
        currentFrame->RotateAroundAxis(rotationMatrix, theta_in_radians, (*dih).Rmask);
        // DEBUG
        DebugTraj.WriteFrame(debugframenum++,currentParm,*currentFrame);
        // Sanity check
        CheckResidue(currentFrame, *dih, second_atom, &clash);
        rotate_dihedral=false;*/
        //DebugTraj.EndTraj();
        //return 1;
      }
    } // End dihedral rotation loop
    next_dih++;
    // Safety valve - number of defined dihedrals times 2
    if (number_of_rotations > max_rotations) {
      mprinterr("DihedralScan: [%i] # of rotations (%i) exceeds max rotations (%i), exiting.\n",
                frameNum+OUTPUTFRAMESHIFT,number_of_rotations, max_rotations);
//#ifdef DEBUG_DIHEDRALSCAN
//      DebugTraj.EndTraj();
//#endif
      // Return gracefully for now
      break;
      //return 1;
    }
  } // End loop over dihedrals
#ifdef DEBUG_DIHEDRALSCAN
  DebugTraj.EndTraj();
  mprintf("\tNumber of rotations %i, expected %u\n",number_of_rotations,BB_dihedrals.size());
#endif
  // Check the resulting structure
  n_problems = checkStructure.SeparateAction( currentFrame );
  //mprintf("%i\tResulting structure has %i problems.\n",frameNum,n_problems);
  number_of_problems->Add(frameNum, &n_problems);

  return Action::OK;
} 

