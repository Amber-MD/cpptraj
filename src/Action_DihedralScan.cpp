// Action_DihedralScan
#include <cmath>
#include <ctime> // time
#include "Action_DihedralScan.h"
#include "CpptrajStdio.h"
#include "Constants.h" // DEGRAD
#include "DistRoutines.h"
#include "TorsionRoutines.h"

// Activate DEBUG info
//#define DEBUG_DIHEDRALSCAN
#ifdef DEBUG_DIHEDRALSCAN
#include "TrajectoryFile.h"
#endif

// CONSTRUCTOR
Action_DihedralScan::Action_DihedralScan() :
  mode_(IMPOSE),
  outframe_(0),
  interval_(60.0),
  maxVal_(0),
  check_for_clashes_(false),
  max_rotations_(0),
  max_factor_(2),
  cutoff_(0.64), // 0.8^2
  rescutoff_(100.0), // 10.0^2
  backtrack_(5),
  increment_(1),
  max_increment_(360),
  debug_(0),
  CurrentParm_(0),
  number_of_problems_(0)
{} 

Action_DihedralScan::~Action_DihedralScan() {
  outtraj_.EndTraj();
}

void Action_DihedralScan::Help() {
  mprintf("dihedralscan <mask> [{interval|random|impose*}] [all] [phi] [psi]]\n");
  mprintf("\t'*' denotes default.\n");
  mprintf("\tOptions for 'random': [rseed <rseed>]\n");
  mprintf("\t\t[ check [cutoff <cutoff>] [rescutoff <rescutoff>]\n");
  mprintf("\t\t  [backtrack <backtrack>] [increment <increment>] [maxfactor <max_factor>] ]\n");
  mprintf("\tOptions for 'interval': <interval deg> [outtraj <filename> [<outfmt>]]\n");
  mprintf("\tOptions for 'impose': <impose deg>\n");
}

// Action_DihedralScan::init()
Action::RetType Action_DihedralScan::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  TrajectoryFile::TrajFormatType outfmt = TrajectoryFile::UNKNOWN_TRAJ;
  Topology* outtop = 0;
  int iseed = -1;

  debug_ = debugIn;
  // Get Keywords - first determine mode
  if (actionArgs.hasKey("random"))
    mode_ = RANDOM;
  else if (actionArgs.hasKey("interval"))
    mode_ = INTERVAL;
  else if (actionArgs.hasKey("impose"))
    mode_ = IMPOSE;
  else
    mode_ = IMPOSE;
  // Determine which angles to search for
  SearchFor_[PHI] = actionArgs.hasKey("phi");
  SearchFor_[PSI] = actionArgs.hasKey("psi");
  // If nothing is enabled or 'all' specified, enable all
  bool nothing_enabled = true;
  for (int i = 0; i < N_ANGLE; i++) {
    if (SearchFor_[i]) {
      nothing_enabled = false;
      break;
    }
  }
  if (nothing_enabled || actionArgs.hasKey("all"))
  for (int i = 0; i < N_ANGLE; i++)
    SearchFor_[i] = true;
  // Get mask
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );
  // For impose/interval, get value to impose/shift by
  if ( mode_ != RANDOM ) 
    interval_ = actionArgs.getNextDouble(60.0);
  // For interval, set max rotations and set up output trajectory
  if (mode_ == INTERVAL) {
    maxVal_ = (int) (360.0 / interval_);
    if (maxVal_ < 0) maxVal_ = -maxVal_;
    outfilename_ = actionArgs.GetStringKey("outtraj");
    if (!outfilename_.empty()) {
      outfmt = TrajectoryFile::GetFormatFromArg( actionArgs );
      outtop = PFL->GetParm( actionArgs );
      if (outtop == 0) {
        mprinterr("Error: dihedralscan: No topology for output traj.\n");
        return Action::ERR;
      }
    }
  }
  // Get 'random' args
  if (mode_ == RANDOM) {
    check_for_clashes_ = actionArgs.hasKey("check");
    cutoff_ = actionArgs.getKeyDouble("cutoff",0.8);
    rescutoff_ = actionArgs.getKeyDouble("rescutoff",10.0);
    backtrack_ = actionArgs.getKeyInt("backtrack",4);
    increment_ = actionArgs.getKeyInt("increment",1);
    max_factor_ = actionArgs.getKeyInt("maxfactor",2);
    // Check validity of args
    if (cutoff_ < SMALL) {
      mprinterr("Error: cutoff too small.\n");
      return Action::ERR;
    }
    if (rescutoff_ < SMALL) {
      mprinterr("Error: rescutoff too small.\n");
      return Action::ERR;
    }
    if (backtrack_ < 0) {
      mprinterr("Error: backtrack value must be >= 0\n");
      return Action::ERR;
    }
    if ( increment_<1 || (360 % increment_)!=0 ) {
      mprinterr("Error: increment must be a factor of 360.\n");
      return Action::ERR;
    }
    // Calculate max increment
    max_increment_ = 360 / increment_;
    // Seed random number gen
    iseed = actionArgs.getKeyInt("rseed",-1);
    RN_.rn_set( iseed );
  }
  // Output file for # of problems
  std::string problemFile = actionArgs.GetStringKey("out");
  // Dataset to store number of problems
  number_of_problems_ = DSL->AddSet(DataSet::INT, actionArgs.GetStringNext(),"Nprob");
  if (number_of_problems_==0) return Action::ERR;
  // Add dataset to data file list
  DFL->AddSetToFile(problemFile,number_of_problems_);

  mprintf("    DIHEDRALSCAN: Dihedrals in mask [%s]\n",Mask1_.MaskString());
  switch (mode_) {
    case RANDOM:
      mprintf("\tDihedrals will be rotated to random values.\n");
      if (iseed==-1)
        mprintf("\tRandom number generator will be seeded using time.\n");
      else
        mprintf("\tRandom number generator will be seeded using %i\n",iseed);
      if (check_for_clashes_) {
        mprintf("\tWill attempt to recover from bad steric clashes.\n");
        mprintf("\tAtom cutoff %.2f, residue cutoff %.2f, backtrack = %i\n",
                cutoff_, rescutoff_, backtrack_);
        mprintf("\tWhen clashes occur dihedral will be incremented by %i\n",increment_);
        mprintf("\tMax # attempted rotations = %i times number dihedrals.\n",
                max_factor_);
      }
      break;
    case INTERVAL:
      mprintf("\tDihedrals will be rotated at intervals of %.2f degrees.\n",
              interval_);
      if (!outfilename_.empty())
        mprintf("\tCoordinates output to %s, format %s\n",outfilename_.c_str(), 
                TrajectoryFile::FormatString(outfmt));
      break;
    case IMPOSE:
      mprintf("\tA value of %.2f degrees will be imposed on selected dihedrals.\n",
               interval_);
  }
  // Setup output trajectory
  if (!outfilename_.empty()) {
    if (outtraj_.SetupTrajWrite(outfilename_, 0, outtop, outfmt)) return Action::ERR;
    outframe_ = 0;
  } 
  // Square cutoffs to compare to dist^2 instead of dist
  cutoff_ *= cutoff_;
  rescutoff_ *= rescutoff_;
  // Increment backtrack by 1 since we need to skip over current res
  ++backtrack_;
  // Initialize CheckStructure
  ArgList cs_args("noimage nobondcheck");
  if (checkStructure_.Init( cs_args, PFL, FL, DSL, DFL, debug_) != Action::OK) {
    mprinterr("Error: Could not set up structure check for DIHEDRALSCAN.\n");
    return Action::ERR;
  }
  return Action::OK;
}

// GetBondedAtomIdx()
static int GetBondedAtomIdx(Topology const& topIn, int atom, NameType const& nameIn) {
  for (Atom::bond_iterator bndatm = topIn[atom].bondbegin();
                           bndatm != topIn[atom].bondend(); ++bndatm)
  {
    if ( topIn[*bndatm].Name() == nameIn ) return *bndatm;
  }
  return -1;
}

/** \return array containing atom indices involved with specified dihedral.
  * Always check that an atom is in Mask1.
  * \param topIn Topology
  * \param atom1 Index of 2nd atom in dihedral
  * \param A0 name of 1st atom in dihedral
  * \param A2 name of 3rd atom in dihedral
  * \param A3 name of 4th atom in dihedral
  */
// TODO: Should atom0 and atom3 be checked for in mask?
// TODO: Should atom0 and atom3 always be looked for?
int Action_DihedralScan::GetDihedralIdxs(int* idxs, Topology const& topIn, 
                                         int atom1, NameType const& A0, 
                                         NameType const& A2, NameType const& A3)
{
  idxs[0] = -1; // atom0
//idxs[1] = -1; // atom2
  idxs[2] = -1; // atom3
  if (mode_ == IMPOSE) { 
    // Get index of atom0
    idxs[0] = GetBondedAtomIdx(topIn, atom1, A0);
    if (idxs[0] == -1) return 1;
    //if (!Mask1_.AtomInCharMask(idxs[0])) return 1;
  }
  // Get index of atom2
  idxs[1] = GetBondedAtomIdx(topIn, atom1, A2);
  if (idxs[1] == -1) return 1;
  if (!Mask1_.AtomInCharMask(idxs[1])) return 1;
  if (mode_ == IMPOSE) {
    // Get index of atom 3
    idxs[2] = GetBondedAtomIdx(topIn, idxs[1], A3);
    if (idxs[2] == -1) return 1;
    //if (!Mask1_.AtomInCharMask(idxs[2])) return 1;
  }
  return 0;
}

// VisitAtom()
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
/** Determine from selected mask atoms which dihedrals will be rotated. */
Action::RetType Action_DihedralScan::Setup(Topology* currentParm, Topology** parmAddress) {
  DihedralScanType dst;
  // Set up Character mask
  if ( currentParm->SetupCharMask( Mask1_ ) ) return Action::ERR;
  Mask1_.MaskInfo();
  if (Mask1_.None()) {
    mprintf("    Error: DihedralScan::setup: Mask has no atoms.\n");
    return Action::ERR;
  }
  // For now just focus on backbone phi/psi dihedrals:
  //   C-N-CA-C  N-CA-C-N
  int idxs[3];
  std::vector<char> Visited( currentParm->Natom(), 'F' );
  for (int atom1 = 0; atom1 < currentParm->Natom(); atom1++) {
    if (Mask1_.AtomInCharMask(atom1)) {
      // PHI: C-N-CA-C
      if ((*currentParm)[atom1].Name() == "N   ") {
        if (!SearchFor_[PHI]) continue;
        if (GetDihedralIdxs(idxs, *currentParm, atom1, "C", "CA", "C")) continue;
      // PSI: N-CA-C-N
      } else if ((*currentParm)[atom1].Name() == "CA  ") {
        if (!SearchFor_[PSI]) continue;
        if (GetDihedralIdxs(idxs, *currentParm, atom1, "N", "C", "N")) continue;
      } else
        continue;
      // Set up mask of atoms that will move upon rotation of dihedral.
      // Also set up mask of atoms in this residue that will not move
      // upon rotation of dihedral, including atom2
      dst.Rmask.ResetMask();
      Visited.assign( currentParm->Natom(), 'F' );
      // Mark atom1 as already visited
      Visited[atom1] = 'T';
      for (Atom::bond_iterator bndatm = (*currentParm)[idxs[1]].bondbegin();
                               bndatm != (*currentParm)[idxs[1]].bondend(); ++bndatm)
      {
        if ( *bndatm != atom1 )
          VisitAtom( *currentParm, *bndatm, Visited );
      }
      dst.checkAtoms.clear();
      int a1res = (*currentParm)[atom1].ResNum();
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
      dst.checkAtoms.push_back( idxs[1] ); // atom1 already included from tmpMask
      dst.atom0 = idxs[0];
      dst.atom1 = atom1;
      dst.atom2 = idxs[1];
      dst.atom3 = idxs[2];
      // Since only the second atom and atoms it is bonded to move during 
      // rotation, base the check on the residue of the second atom.
      dst.resnum = (*currentParm)[idxs[1]].ResNum();
      BB_dihedrals_.push_back(dst);
    } // END if atom in char mask
  }

  // DEBUG: List defined dihedrals
  if (debug_>0) {
    mprintf("DEBUG: Dihedrals:");
    if (mode_ != IMPOSE)
      mprintf(" (central two atoms only)");
    mprintf("\n");
    for (unsigned int dih = 0; dih < BB_dihedrals_.size(); dih++) {
      mprintf("\t");
      std::string aname = currentParm->TruncResAtomName(BB_dihedrals_[dih].atom0);
      if (!aname.empty()) mprintf("%s-", aname.c_str());
      aname = currentParm->TruncResAtomName(BB_dihedrals_[dih].atom1);
      mprintf("%s-", aname.c_str());
      aname = currentParm->TruncResAtomName(BB_dihedrals_[dih].atom2);
      mprintf("%s", aname.c_str());
      aname = currentParm->TruncResAtomName(BB_dihedrals_[dih].atom3);
      if (!aname.empty()) mprintf("-%s", aname.c_str());
      mprintf("\n");
      if (debug_>1) {
        mprintf("\t\tCheckAtoms=");
        for (std::vector<int>::iterator ca = BB_dihedrals_[dih].checkAtoms.begin();
                                        ca != BB_dihedrals_[dih].checkAtoms.end(); ca++)
        {
          mprintf(" %i",*ca + 1);
        }
        mprintf("\n");
      }
      if (debug_>2) {
        mprintf("\t\t");
        BB_dihedrals_[dih].Rmask.PrintMaskAtoms("Rmask:");
      }
    }
  }
  // Set up CheckStructure for this parm
  if (checkStructure_.Setup(currentParm, &currentParm) != Action::OK)
    return Action::ERR;

  // Set the overall max number of rotations to try
  max_rotations_ = (int) BB_dihedrals_.size();
  max_rotations_ *= max_factor_;

  // Set up simple structure check. First step is coarse; check distances 
  // between a certain atom in each residue (first, COM, CA, some other atom?)
  // to see if residues are in each others neighborhood. Second step is to 
  // check the atoms in each close residue.
  if (check_for_clashes_) {
    ResidueCheckType rct;
    int Nres = currentParm->FinalSoluteRes();
    for (int res = 0; res < Nres; res++) {
      rct.resnum = res;
      rct.start = currentParm->ResFirstAtom( res );
      rct.stop = currentParm->ResLastAtom( res );
      rct.checkatom = rct.start;
      ResCheck_.push_back(rct);
    }
  }
  CurrentParm_ = currentParm;
  return Action::OK;  
}

// Action_DihedralScan::CheckResidue()
/** \return 1 if a new dihedral should be tried, 0 if no clashes, -1 if
  * \return further rotations will not help.
  */
int Action_DihedralScan::CheckResidue( Frame const& FrameIn, DihedralScanType &dih, 
                                       int nextres, double *clash ) 
{
  int resnumIn = dih.resnum;
  int rstart = ResCheck_[ resnumIn ].start;
  int rstop = ResCheck_[ resnumIn ].stop;
  int rcheck = ResCheck_[ resnumIn ].checkatom;
  // Check for clashes with self
#ifdef DEBUG_DIHEDRALSCAN
  mprintf("\tChecking residue %i\n",resnumIn+1);
  mprintf("\tATOMS %i to %i\n",rstart+1,rstop);
#endif
  for (int atom1 = rstart; atom1 < rstop - 1; atom1++) {
    for (int atom2 = atom1 + 1; atom2 < rstop; atom2++) {
      double atomD2 = DIST2_NoImage(FrameIn.XYZ(atom1), FrameIn.XYZ(atom2));
      if (atomD2 < cutoff_) {
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
    int rstart2 = ResCheck_[ res ].start;
    int rstop2 = ResCheck_[ res ].stop;
    int rcheck2 = ResCheck_[ res ].checkatom;
    double resD2 = DIST2_NoImage(FrameIn.XYZ(rcheck), FrameIn.XYZ(rcheck2));
    // If residues are close enough check each atom
    if (resD2 < rescutoff_) { 
#ifdef DEBUG_DIHEDRALSCAN
      mprintf("\tRES %i ATOMS %i to %i\n",res+1,rstart2+2,rstop2);
#endif
      for (int atom1 = rstart; atom1 < rstop; atom1++) {
        for (int atom2 = rstart2; atom2 < rstop2; atom2++) {
          double D2 = DIST2_NoImage(FrameIn.XYZ(atom1), FrameIn.XYZ(atom2));
          if (D2 < cutoff_) {
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

// Action_DihedralScan::RandomizeAngles()
void Action_DihedralScan::RandomizeAngles(Frame& currentFrame) {
  Matrix_3x3 rotationMatrix;
#ifdef DEBUG_DIHEDRALSCAN
  // DEBUG
  int debugframenum=0;
  Trajout DebugTraj;
  DebugTraj.SetupTrajWrite("debugtraj.nc",0,currentParm,TrajectoryFile::AMBERNETCDF);
  DebugTraj.WriteFrame(debugframenum++,currentParm,currentFrame);
#endif
  int next_resnum;
  int bestLoop = 0;
  int number_of_rotations = 0;

  std::vector<DihedralScanType>::iterator next_dih = BB_dihedrals_.begin();
  next_dih++;
  for (std::vector<DihedralScanType>::iterator dih = BB_dihedrals_.begin();
                                               dih != BB_dihedrals_.end();
                                               dih++)
  {
    ++number_of_rotations;
    // Get the residue atom of the next dihedral. Residues up to and
    // including this residue will be checked for bad clashes 
    if (next_dih!=BB_dihedrals_.end()) 
      next_resnum = (*next_dih).resnum;
    else
      next_resnum = (*dih).resnum-1;
    // Set axis of rotation
    Vec3 axisOfRotation = currentFrame.SetAxisOfRotation((*dih).atom1, (*dih).atom2);
    // Generate random value to rotate by in radians
    // Guaranteed to rotate by at least 1 degree.
    // NOTE: could potentially rotate 360 - prevent?
    // FIXME: Just use 2PI and rn_gen, get everything in radians
    double theta_in_degrees = ((int)(RN_.rn_gen()*100000) % 360) + 1;
    double theta_in_radians = theta_in_degrees * DEGRAD;
    // Calculate rotation matrix for random theta
    rotationMatrix.CalcRotationMatrix(axisOfRotation, theta_in_radians);
    int loop_count = 0;
    double clash = 0;
    double bestClash = 0;
    if (debug_>0) mprintf("DEBUG: Rotating res %8i:\n",(*dih).resnum+1);
    bool rotate_dihedral = true;
    while (rotate_dihedral) {
      if (debug_>0) {
        mprintf("\t%8i %8i%4s %8i%4s, +%.2lf degrees (%i).\n",(*dih).resnum+1,
                (*dih).atom1+1, (*CurrentParm_)[(*dih).atom1].c_str(),
                (*dih).atom2+1, (*CurrentParm_)[(*dih).atom2].c_str(),
                theta_in_degrees,loop_count);
      }
      // Rotate around axis
      currentFrame.Rotate(rotationMatrix, (*dih).Rmask);
#ifdef DEBUG_DIHEDRALSCAN
      // DEBUG
      DebugTraj.WriteFrame(debugframenum++,currentParm,*currentFrame);
#endif
      // If we dont care about sterics exit here
      if (!check_for_clashes_) break;
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
          mprintf("\tTrying dihedral increments of +%i\n",increment_);
        // Instead of a new random dihedral, try increments
        theta_in_degrees = (double)increment_;
        theta_in_radians = theta_in_degrees * DEGRAD;
        // Calculate rotation matrix for new theta
        rotationMatrix.CalcRotationMatrix(axisOfRotation, theta_in_radians);
      }
      ++loop_count;
      if (loop_count == max_increment_) {
        if (debug_>0)
          mprintf("%i iterations! Best clash= %.3lf at %i\n",max_increment_,
                  sqrt(bestClash),bestLoop);
        for (int bt = 0; bt < backtrack_; bt++)
          dih--;
        next_dih = dih;
        next_dih++;
        if (debug_>0)
          mprintf("\tCannot resolve clash with further rotations, trying previous %i again.\n",
                  backtrack_ - 1);
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
    if (number_of_rotations > max_rotations_) {
      mprinterr("Error: DihedralScan: # of rotations (%i) exceeds max rotations (%i), exiting.\n",
                number_of_rotations, max_rotations_);
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
  mprintf("\tNumber of rotations %i, expected %u\n",number_of_rotations,BB_dihedrals_.size());
#endif
}

// Action_DihedralScan::IntervalAngles()
void Action_DihedralScan::IntervalAngles(Frame& currentFrame) {
  Matrix_3x3 rotationMatrix;
  double theta_in_radians = interval_ * DEGRAD;
  // Write original frame
  if (!outfilename_.empty())
    outtraj_.WriteFrame(outframe_++, CurrentParm_, currentFrame);
  for (std::vector<DihedralScanType>::iterator dih = BB_dihedrals_.begin();
                                               dih != BB_dihedrals_.end();
                                               dih++)
  {
    // Set axis of rotation
    Vec3 axisOfRotation = currentFrame.SetAxisOfRotation((*dih).atom1, (*dih).atom2);
    // Calculate rotation matrix for interval 
    rotationMatrix.CalcRotationMatrix(axisOfRotation, theta_in_radians);
    if (debug_ > 0) {
      std::string a1name = CurrentParm_->TruncResAtomName( (*dih).atom1 );
      std::string a2name = CurrentParm_->TruncResAtomName( (*dih).atom2 );
      mprintf("\tRotating Dih %s-%s by %.2f deg %i times.\n",
               a1name.c_str(), a2name.c_str(), interval_, maxVal_); 
    }
    for (int rot = 0; rot < maxVal_; ++rot) {
      // Rotate around axis
      currentFrame.Rotate(rotationMatrix, (*dih).Rmask);
      // Write output trajectory
      if (outtraj_.TrajIsOpen())
        outtraj_.WriteFrame(outframe_++, CurrentParm_, currentFrame);
    }
  }
}

// Action_DihedralScan::ImposeAngles()
void Action_DihedralScan::ImposeAngles(Frame& currentFrame) {
  Matrix_3x3 rotationMatrix;
  double theta_in_radians = interval_ * DEGRAD;
  for (std::vector<DihedralScanType>::iterator dih = BB_dihedrals_.begin();
                                               dih != BB_dihedrals_.end();
                                               dih++)
  {
    // Calculate current value of dihedral
    double torsion = Torsion( currentFrame.XYZ( (*dih).atom0 ),
                              currentFrame.XYZ( (*dih).atom1 ),
                              currentFrame.XYZ( (*dih).atom2 ),
                              currentFrame.XYZ( (*dih).atom3 ) );
    // Calculate delta needed to get to theta
    double delta = theta_in_radians - torsion;
    // Set axis of rotation
    Vec3 axisOfRotation = currentFrame.SetAxisOfRotation((*dih).atom1, (*dih).atom2);
    // Calculate rotation matrix for delta 
    rotationMatrix.CalcRotationMatrix(axisOfRotation, delta);
    if (debug_ > 0) {
      std::string a0name = CurrentParm_->TruncResAtomName( (*dih).atom0 );
      std::string a1name = CurrentParm_->TruncResAtomName( (*dih).atom1 );
      std::string a2name = CurrentParm_->TruncResAtomName( (*dih).atom2 );
      std::string a3name = CurrentParm_->TruncResAtomName( (*dih).atom3 );
      mprintf("\tRotating Dih %s-%s-%s-%s (@%.2f) by %.2f deg to get to %.2f.\n",
               a0name.c_str(), a1name.c_str(), a2name.c_str(), a3name.c_str(),
               torsion*RADDEG, delta*RADDEG, interval_); 
    }
    // Rotate around axis
    currentFrame.Rotate(rotationMatrix, (*dih).Rmask);
  }
}

// Action_DihedralScan::action()
Action::RetType Action_DihedralScan::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  switch (mode_) {
    case RANDOM: RandomizeAngles(*currentFrame); break;
    case INTERVAL: IntervalAngles(*currentFrame); break;
    case IMPOSE: ImposeAngles(*currentFrame); break;
  }
  // Check the resulting structure
  int n_problems = checkStructure_.CheckFrame( frameNum+1, *currentFrame );
  //mprintf("%i\tResulting structure has %i problems.\n",frameNum,n_problems);
  number_of_problems_->Add(frameNum, &n_problems);

  return Action::OK;
} 

