#include <cmath>
#include "Exec_PermuteDihedrals.h"
#include "CpptrajStdio.h"
#include "DihedralSearch.h"
#include "Constants.h" // SMALL, DEGRAD
#include "DistRoutines.h"
#include "TorsionRoutines.h"
// Activate DEBUG info
//#define DEBUG_PERMUTEDIHEDRALS

// CONSTRUCTOR
Exec_PermuteDihedrals::Exec_PermuteDihedrals() : Exec(GENERAL),
  mode_(INTERVAL),
  debug_(0),
  outframe_(0),
  crdout_(0),
  check_for_clashes_(false),
  checkAllResidues_(false),
  max_factor_(2),
  cutoff_(0.64), // 0.8^2
  rescutoff_(100.0), // 10.0^2
  backtrack_(5),
  increment_(1),
  max_increment_(360),
  number_of_problems_(0)
{} 

// Exec_PermuteDihedrals::Help()
void Exec_PermuteDihedrals::Help() const {
  mprintf("\tcrdset <COORDS set> resrange <range> [{interval | random}]\n"
          "\t[outtraj <filename> [<outfmt>]] [crdout <output COORDS>] [<dihedral types>]\n"
          "  Options for 'random':\n"
          "\t[rseed <rseed>] [out <# problems file> [<set name>]]\n"
          "\t[ check [cutoff <cutoff>] [rescutoff <rescutoff>] [checkallresidues]\n"
          "\t  [backtrack <backtrack>] [increment <increment>] [maxfactor <max_factor>] ]\n"
          "  Options for 'interval':\n"
          "\t<interval deg>\n"
          "  <dihedral types> = ");
  DihedralSearch::ListKnownTypes();
  mprintf("  Rotate specified dihedral(s) in given COORDS set by specific interval\n"
          "  or to random values.\n");
}

// Exec_PermuteDihedrals::Execute()
Exec::RetType Exec_PermuteDihedrals::Execute(CpptrajState& State, ArgList& argIn) {
  mode_ = INTERVAL;
  // Get Keywords - first determine mode
  if (argIn.hasKey("random"))
    mode_ = RANDOM;
  else if (argIn.hasKey("interval"))
    mode_ = INTERVAL;
  // Get input COORDS set
  std::string setname = argIn.GetStringKey("crdset");
  if (setname.empty()) {
    mprinterr("Error: Specify COORDS dataset name with 'crdset'.\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL().FindCoordsSet( setname );
  if (CRD == 0) {
    mprinterr("Error: Could not find COORDS set '%s'\n", setname.c_str());
    return CpptrajState::ERR;
  }
  mprintf("    PERMUTEDIHEDRALS: Using COORDS '%s'\n", CRD->legend());

  // Get residue range
  Range resRange;
  resRange.SetRange(argIn.GetStringKey("resrange"));
  if (!resRange.Empty())
    resRange.ShiftBy(-1); // User res args start from 1
  mprintf("\tPermutating dihedrals in");
  if (resRange.Empty())
    mprintf(" all solute residues.\n");
  else
    mprintf(" residue range [%s]\n", resRange.RangeArg());

  // Determine which angles to search for
  DihedralSearch dihSearch;
  dihSearch.SearchForArgs(argIn);
  // If nothing is enabled, enable all 
  dihSearch.SearchForAll();
  mprintf("\tSearching for types:");
  dihSearch.PrintTypes();
  mprintf("\n");

  // Setup output trajectory
  outframe_ = 0; 
  std::string outfilename = argIn.GetStringKey("outtraj");
  if (!outfilename.empty()) {
    mprintf("\tCoordinates output to '%s'\n", outfilename.c_str());
    Topology* outtop = State.DSL().GetTopology( argIn );
    if (outtop == 0) {
      mprinterr("Error: No topology for output traj.\n");
      return CpptrajState::ERR;
    }
    // Setup output trajectory FIXME: Correct frames for # of rotations
    if (outtraj_.PrepareTrajWrite(outfilename, argIn, CRD->TopPtr(), CRD->CoordsInfo(),
                                  CRD->Size(), TrajectoryFile::UNKNOWN_TRAJ))
      return CpptrajState::ERR;
  }

  // Setup output coords
  outfilename = argIn.GetStringKey("crdout");
  if (!outfilename.empty()) {
    mprintf("\tCoordinates saved to set '%s'\n", outfilename.c_str());
    crdout_ = (DataSet_Coords_CRD*)State.DSL().AddSet(DataSet::COORDS, outfilename);
    if (crdout_ == 0) return CpptrajState::ERR;
    crdout_->CoordsSetup( CRD->Top(), CRD->CoordsInfo() );
  }

  // Get specific mode options.
  double interval_in_deg = 60.0;
  if ( mode_ == INTERVAL ) {
    interval_in_deg = argIn.getNextDouble(60.0);
    mprintf("\tDihedrals will be rotated at intervals of %.2f degrees.\n", interval_in_deg);
  } else if (mode_ == RANDOM) {
    check_for_clashes_ = argIn.hasKey("check");
    checkAllResidues_ = argIn.hasKey("checkallresidues");
    cutoff_ = argIn.getKeyDouble("cutoff",0.8);
    rescutoff_ = argIn.getKeyDouble("rescutoff",10.0);
    backtrack_ = argIn.getKeyInt("backtrack",4);
    increment_ = argIn.getKeyInt("increment",1);
    max_factor_ = argIn.getKeyInt("maxfactor",2);
    int iseed = argIn.getKeyInt("rseed",-1);
    // Output file for # of problems
    DataFile* problemFile = State.DFL().AddDataFile(argIn.GetStringKey("out"), argIn);
    // Dataset to store number of problems
    number_of_problems_ = State.DSL().AddSet(DataSet::INTEGER, argIn.GetStringNext(),"Nprob");
    if (number_of_problems_==0) return CpptrajState::ERR;
   // Add dataset to data file list
    if (problemFile != 0) problemFile->AddDataSet(number_of_problems_);
    // Check validity of args
    if (cutoff_ < Constants::SMALL) {
      mprinterr("Error: cutoff too small.\n");
      return CpptrajState::ERR;
    }
    if (rescutoff_ < Constants::SMALL) {
      mprinterr("Error: rescutoff too small.\n");
      return CpptrajState::ERR;
    }
    if (backtrack_ < 0) {
      mprinterr("Error: backtrack value must be >= 0\n");
      return CpptrajState::ERR;
    }
    if ( increment_<1 || (360 % increment_)!=0 ) {
      mprinterr("Error: increment must be a factor of 360.\n");
      return CpptrajState::ERR;
    }
    // Calculate max increment
    max_increment_ = 360 / increment_;
    // Seed random number gen
    RN_.rn_set( iseed );
    // Print info
    mprintf("\tDihedrals will be rotated to random values.\n");
    if (iseed==-1)
      mprintf("\tRandom number generator will be seeded using time.\n");
    else
      mprintf("\tRandom number generator will be seeded using %i\n",iseed);
    if (check_for_clashes_) {
      mprintf("\tWill attempt to recover from bad steric clashes.\n");
      if (checkAllResidues_)
        mprintf("\tAll residues will be checked.\n");
      else
        mprintf("\tResidues up to the currenly rotating dihedral will be checked.\n");
      mprintf("\tAtom cutoff %.2f, residue cutoff %.2f, backtrack = %i\n",
              cutoff_, rescutoff_, backtrack_);
      mprintf("\tWhen clashes occur dihedral will be incremented by %i\n",increment_);
      mprintf("\tMax # attempted rotations = %i times number dihedrals.\n",
              max_factor_);
    }
    // Square cutoffs to compare to dist^2 instead of dist
    cutoff_ *= cutoff_;
    rescutoff_ *= rescutoff_;
    // Increment backtrack by 1 since we need to skip over current res
    ++backtrack_;
    // Initialize CheckStructure
    if (checkStructure_.SeparateInit( false, "*", "", "", 0.8, 1.15, false, State.DFL() )) {
      mprinterr("Error: Could not set up structure check.\n");
      return CpptrajState::ERR;
    }
    // Set up CheckStructure for this parm (false = nobondcheck)
    if (checkStructure_.SeparateSetup(CRD->Top(),
                                      CRD->CoordsInfo().TrajBox().Type(), false) != Action::OK)
      return CpptrajState::ERR;
  }

  // Determine from selected mask atoms which dihedrals will be rotated.
  PermuteDihedralsType dst;
  // If range is empty (i.e. no resrange arg given) look through all 
  // solute residues.
  Range actualRange;
  if (resRange.Empty())
    actualRange = CRD->Top().SoluteResidues();
  else 
    actualRange = resRange;
  // Search for dihedrals
  if (dihSearch.FindDihedrals(CRD->Top(), actualRange))
    return CpptrajState::ERR;
  // For each found dihedral, set up mask of atoms that will move upon 
  // rotation. Also set up mask of atoms in this residue that will not
  // move, including atom2.
  if (debug_>0)
    mprintf("DEBUG: Dihedrals:\n");
  for (DihedralSearch::mask_it dih = dihSearch.begin();
                               dih != dihSearch.end(); ++dih)
  {
    dst.checkAtoms.clear();
    // Set mask of atoms that will move during dihedral rotation.
    dst.Rmask = DihedralSearch::MovingAtoms(CRD->Top(), dih->A1(), dih->A2());
    // If randomly rotating angles, check for atoms that are in the same
    // residue as A1 but will not move. They need to be checked for clashes
    // since further rotations will not help them.
    if (mode_ == RANDOM && check_for_clashes_) {
      CharMask cMask( dst.Rmask.ConvertToCharMask(), dst.Rmask.Nselected() );
      int a1res = CRD->Top()[dih->A1()].ResNum();
      for (int maskatom = CRD->Top().Res(a1res).FirstAtom();
               maskatom < CRD->Top().Res(a1res).LastAtom(); ++maskatom)
        if (!cMask.AtomInCharMask(maskatom))
          dst.checkAtoms.push_back( maskatom );
      dst.checkAtoms.push_back(dih->A1()); // TODO: Does this need to be added first?
      // Since only the second atom and atoms it is bonded to move during 
      // rotation, base the check on the residue of the second atom.
      dst.resnum = a1res;
    }
    dst.atom0 = dih->A0(); // FIXME: This duplicates info
    dst.atom1 = dih->A1();
    dst.atom2 = dih->A2();
    dst.atom3 = dih->A3();
    BB_dihedrals_.push_back(dst);
    // DEBUG: List dihedral info.
    if (debug_ > 0) {
      mprintf("\t%s-%s-%s-%s\n", 
              CRD->Top().TruncResAtomName(dih->A0()).c_str(),
              CRD->Top().TruncResAtomName(dih->A1()).c_str(),
              CRD->Top().TruncResAtomName(dih->A2()).c_str(),
              CRD->Top().TruncResAtomName(dih->A3()).c_str() );
      if (debug_ > 1 && mode_ == RANDOM && check_for_clashes_) {
        mprintf("\t\tCheckAtoms=");
        for (std::vector<int>::const_iterator ca = dst.checkAtoms.begin();
                                              ca != dst.checkAtoms.end(); ++ca)
          mprintf(" %i", *ca + 1);
        mprintf("\n");
      }
      if (debug_ > 2) {
        mprintf("\t\t");
        dst.Rmask.PrintMaskAtoms("Rmask:");
      }
    }
  }

  // Set up simple structure check. First step is coarse; check distances 
  // between a certain atom in each residue (first, COM, CA, some other atom?)
  // to see if residues are in each others neighborhood. Second step is to 
  // check the atoms in each close residue.
  if (check_for_clashes_) {
    ResidueCheckType rct;
    int res = 0;
    for (Topology::res_iterator residue = CRD->Top().ResStart();
                                residue != CRD->Top().ResEnd(); ++residue)
    {
      rct.resnum = res++;
      rct.start = residue->FirstAtom();
      rct.stop = residue->LastAtom();
      rct.checkatom = rct.start;
      ResCheck_.push_back(rct);
    }
  }

  // Perform dihedral permute
  Frame currentFrame = CRD->AllocateFrame();
  for (unsigned int set = 0; set != CRD->Size(); set++)
  {
    CRD->GetFrame(set, currentFrame);
    int n_problems = 0;
    switch (mode_) {
      case RANDOM:
        RandomizeAngles(currentFrame, CRD->Top());
        // Check the resulting structure
        n_problems = checkStructure_.CheckOverlap( set+1, currentFrame, CRD->Top() );
        //mprintf("%i\tResulting structure has %i problems.\n",frameNum,n_problems);
        number_of_problems_->Add(set, &n_problems);
        if (outtraj_.IsInitialized()) outtraj_.WriteSingle(outframe_++, currentFrame);
        if (crdout_ != 0) crdout_->AddFrame( currentFrame );
        break;
      case INTERVAL: IntervalAngles(currentFrame, CRD->Top(), interval_in_deg); break;
    }
  }
  if (outtraj_.IsInitialized()) outtraj_.EndTraj();
  return CpptrajState::OK;
}

// -----------------------------------------------------------------------------
// Exec_PermuteDihedrals::IntervalAngles()
void Exec_PermuteDihedrals::IntervalAngles(Frame const& frameIn, Topology const& topIn,
                                           double interval_in_deg)
{
  Matrix_3x3 rotationMatrix;
  double theta_in_radians = interval_in_deg * Constants::DEGRAD;
  int maxVal = (int) (360.0 / interval_in_deg);
  if (maxVal < 0) maxVal = -maxVal;
  // Write original frame
  if (outtraj_.IsInitialized())
    outtraj_.WriteSingle(outframe_++, frameIn);
  if (crdout_ != 0)
    crdout_->AddFrame( frameIn );
  Frame currentFrame = frameIn;
  for (std::vector<PermuteDihedralsType>::const_iterator dih = BB_dihedrals_.begin();
                                                     dih != BB_dihedrals_.end();
                                                     ++dih)
  {
    // Set axis of rotation
    Vec3 axisOfRotation = currentFrame.SetAxisOfRotation(dih->atom1, dih->atom2);
    // Calculate rotation matrix for interval 
    rotationMatrix.CalcRotationMatrix(axisOfRotation, theta_in_radians);
    if (debug_ > 0) {
      mprintf("\tRotating Dih %s-%s by %.2f deg %i times.\n",
               topIn.TruncResAtomName( dih->atom1 ).c_str(), 
               topIn.TruncResAtomName( dih->atom2 ).c_str(), interval_in_deg, maxVal); 
    }
    for (int rot = 0; rot != maxVal; ++rot) {
      // Rotate around axis
      currentFrame.Rotate(rotationMatrix, dih->Rmask);
      // Write output trajectory
      if (outtraj_.IsInitialized())
        outtraj_.WriteSingle(outframe_++, currentFrame);
      if (crdout_ != 0)
        crdout_->AddFrame( currentFrame );
    }
  }
}

// -----------------------------------------------------------------------------
// Exec_PermuteDihedrals::CheckResidue()
/** \return 1 if a new dihedral should be tried, 0 if no clashes
  * \return -1 if further rotations will not help.
  */
int Exec_PermuteDihedrals::CheckResidue( Frame const& FrameIn, Topology const& topIn,
                                         PermuteDihedralsType const& dih, 
                                         int nextres, double& clash ) 
{
  int resnumIn = dih.resnum;
  int rstart = ResCheck_[ resnumIn ].start;
  int rstop = ResCheck_[ resnumIn ].stop;
  int rcheck = ResCheck_[ resnumIn ].checkatom;
  // Check for clashes with self
# ifdef DEBUG_PERMUTEDIHEDRALS
  mprintf("\tChecking residue %i\n",resnumIn+1);
  mprintf("\tATOMS %i to %i\n",rstart+1,rstop);
# endif
  for (int atom1 = rstart; atom1 < rstop - 1; atom1++) {
    for (int atom2 = atom1 + 1; atom2 < rstop; atom2++) {
      // Skip bonded atoms
      bool isBonded = false;
      for (Atom::bond_iterator bndatm = topIn[atom1].bondbegin();
                               bndatm != topIn[atom1].bondend(); ++bndatm)
        if (*bndatm == atom2) {
          isBonded = true;
          break;
        }
      if (!isBonded) {
        double atomD2 = DIST2_NoImage(FrameIn.XYZ(atom1), FrameIn.XYZ(atom2));
        if (atomD2 < cutoff_) {
#         ifdef DEBUG_PERMUTEDIHEDRALS 
          mprintf("\t\tCurrent Res %i Atoms %s and %s are close (%.3lf)\n", resnumIn+1, 
                  topIn.AtomMaskName(atom1).c_str(),
                  topIn.AtomMaskName(atom2).c_str(), sqrt(atomD2));
#         endif
          clash = atomD2;
          return 1;
        }
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
#     ifdef DEBUG_PERMUTEDIHEDRALS
      mprintf("\tRES %i ATOMS %i to %i\n",res+1,rstart2+2,rstop2);
#     endif
      for (int atom1 = rstart; atom1 < rstop; atom1++) {
        for (int atom2 = rstart2; atom2 < rstop2; atom2++) {
          double D2 = DIST2_NoImage(FrameIn.XYZ(atom1), FrameIn.XYZ(atom2));
          if (D2 < cutoff_) {
#           ifdef DEBUG_PERMUTEDIHEDRALS
            mprintf("\t\tResCheck %i Atoms %s and %s are close (%.3lf)\n", res+1,
                    topIn.TruncResAtomName(atom1).c_str(),
                    topIn.TruncResAtomName(atom2).c_str(), sqrt(D2));
#           endif
            clash = D2;
            // If the clash involves any atom that will not be moved by further
            // rotation, indicate it is not possible to resolve clash by
            // more rotation by returning -1.
            //if (atom1 == dih.atom2 || atom1 == dih.atom1) return -1;
            for (std::vector<int>::const_iterator ca = dih.checkAtoms.begin();
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

// Exec_PermuteDihedrals::RandomizeAngles()
void Exec_PermuteDihedrals::RandomizeAngles(Frame& currentFrame, Topology const& topIn) {
  Matrix_3x3 rotationMatrix;
# ifdef DEBUG_PERMUTEDIHEDRALS
  // DEBUG
  int debugframenum=0;
  Trajout_Single DebugTraj;
  DebugTraj.PrepareTrajWrite("debugtraj.nc",ArgList(),(Topology*)&topIn,
                             topIn.ParmCoordInfo(), topIn.Nframes(),
                             TrajectoryFile::AMBERNETCDF);
  DebugTraj.WriteSingle(debugframenum++,currentFrame);
# endif
  int next_resnum;
  int bestLoop = 0;
  int number_of_rotations = 0;
  // Set max number of rotations to try.
  int max_rotations = (int)BB_dihedrals_.size();
  max_rotations *= max_factor_;

  // Loop over all dihedrals
  std::vector<PermuteDihedralsType>::const_iterator next_dih = BB_dihedrals_.begin();
  next_dih++;
  for (std::vector<PermuteDihedralsType>::const_iterator dih = BB_dihedrals_.begin();
                                                     dih != BB_dihedrals_.end(); 
                                                     ++dih, ++next_dih)
  {
    ++number_of_rotations;
    // Get the residue atom of the next dihedral. Residues up to and
    // including this residue will be checked for bad clashes 
    if (next_dih != BB_dihedrals_.end()) 
      next_resnum = next_dih->resnum;
    else
      next_resnum = dih->resnum - 1;
    // Set axis of rotation
    Vec3 axisOfRotation = currentFrame.SetAxisOfRotation(dih->atom1, dih->atom2);
    // Generate random value to rotate by in radians
    // Guaranteed to rotate by at least 1 degree.
    // NOTE: could potentially rotate 360 - prevent?
    // FIXME: Just use 2PI and rn_gen, get everything in radians
    double theta_in_degrees = ((int)(RN_.rn_gen()*100000) % 360) + 1;
    double theta_in_radians = theta_in_degrees * Constants::DEGRAD;
    // Calculate rotation matrix for random theta
    rotationMatrix.CalcRotationMatrix(axisOfRotation, theta_in_radians);
    int loop_count = 0;
    double clash = 0;
    double bestClash = 0;
    if (debug_>0) mprintf("DEBUG: Rotating dihedral %zu res %8i:\n", dih - BB_dihedrals_.begin(),
                          dih->resnum+1);
    bool rotate_dihedral = true;
    while (rotate_dihedral) {
      if (debug_>0) {
        mprintf("\t%8i %12s %12s, +%.2lf degrees (%i).\n",dih->resnum+1,
                topIn.AtomMaskName(dih->atom1).c_str(),
                topIn.AtomMaskName(dih->atom2).c_str(),
                theta_in_degrees,loop_count);
      }
      // Rotate around axis
      currentFrame.Rotate(rotationMatrix, dih->Rmask);
#     ifdef DEBUG_PERMUTEDIHEDRALS
      // DEBUG
      DebugTraj.WriteSingle(debugframenum++,currentFrame);
#     endif
      // If we dont care about sterics exit here
      if (!check_for_clashes_) break;
      // Check resulting structure for issues
      int checkresidue;
      if (!checkAllResidues_)
        checkresidue = CheckResidue(currentFrame, topIn, *dih, next_resnum, clash);
      else
        checkresidue = CheckResidue(currentFrame, topIn, *dih, topIn.Nres(), clash);
      if (checkresidue==0)
        rotate_dihedral = false;
      else if (checkresidue==-1) {
        if (dih - BB_dihedrals_.begin() < 2) {
          mprinterr("Error: Cannot backtrack; initial structure already has clashes.\n");
          number_of_rotations = max_rotations + 1;
        } else {
          dih--; //  0
          dih--; // -1
          next_dih = dih;
          next_dih++;
          if (debug_>0)
            mprintf("\tCannot resolve clash with further rotations, trying previous again.\n");
        }
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
        theta_in_radians = theta_in_degrees * Constants::DEGRAD;
        // Calculate rotation matrix for new theta
        rotationMatrix.CalcRotationMatrix(axisOfRotation, theta_in_radians);
      }
      ++loop_count;
      if (loop_count == max_increment_) {
        if (debug_>0)
          mprintf("%i iterations! Best clash= %.3lf at %i\n",max_increment_,
                  sqrt(bestClash),bestLoop);
        if (dih - BB_dihedrals_.begin() < backtrack_) {
          mprinterr("Error: Cannot backtrack; initial structure already has clashes.\n");
          number_of_rotations = max_rotations + 1;
        } else { 
          for (int bt = 0; bt < backtrack_; bt++)
            dih--;
          next_dih = dih;
          next_dih++;
          if (debug_>0)
            mprintf("\tCannot resolve clash with further rotations, trying previous %i again.\n",
                    backtrack_ - 1);
        }
        break;
        // Calculate how much to rotate back in order to get to best clash
        /*int num_back = bestLoop - 359;
        theta_in_degrees = (double) num_back;
        theta_in_radians = theta_in_degrees * Constants::DEGRAD;
        // Calculate rotation matrix for theta
        calcRotationMatrix(rotationMatrix, axisOfRotation, theta_in_radians);
        // Rotate back to best clash
        frm.Frm().RotateAroundAxis(rotationMatrix, theta_in_radians, dih->Rmask);
        // DEBUG
        DebugTraj.WriteFrame(debugframenum++,currentParm,*currentFrame);
        // Sanity check
        CheckResidue(currentFrame, *dih, second_atom, &clash);
        rotate_dihedral=false;*/
        //DebugTraj.EndTraj();
        //return 1;
      }
    } // End dihedral rotation loop
    // Safety valve - number of defined dihedrals times * maxfactor
    if (number_of_rotations > max_rotations) {
      mprinterr("Error: # of rotations (%i) exceeds max rotations (%i), exiting.\n",
                number_of_rotations, max_rotations);
//#     ifdef DEBUG_PERMUTEDIHEDRALS
//      DebugTraj.EndTraj();
//#     endif
      // Return gracefully for now
      break;
      //return 1;
    }
  } // End loop over dihedrals
# ifdef DEBUG_PERMUTEDIHEDRALS
  DebugTraj.EndTraj();
  mprintf("\tNumber of rotations %i, expected %u\n",number_of_rotations,BB_dihedrals_.size());
# endif
}
