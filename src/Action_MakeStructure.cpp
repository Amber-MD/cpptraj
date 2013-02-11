#include "Action_MakeStructure.h"
#include "CpptrajStdio.h"
#include "Constants.h" // DEGRAD
#include "TorsionRoutines.h"

struct SS_TYPE {
  double phi;
  double psi;
  double phi2;
  double psi2;
  int isTurn;
  const char* name;
};

static SS_TYPE SS[] = {
  {  -57.8,  -47.0,    0.0,   0.0, 0, "alpha"    },
  {   57.8,   47.0,    0.0,   0.0, 0, "left"     },
  {  -75.0,  145.0,    0.0,   0.0, 0, "pp2"      },
  { -100.0,  130.0,    0.0,   0.0, 0, "hairpin"  },
  { -150.0,  155.0,    0.0,   0.0, 0, "extended" },
  {  -60.0,  -30.0,  -90.0,   0.0, 1, "typeI"    }, 	
  {  -60.0,  120.0,   80.0,   0.0, 1, "typeII"   },	
  {  -60.0,  -30.0, -120.0, 120.0, 1, "typeVIII" }, 	
  {   60.0,   30.0,   90.0,   0.0, 1, "typeI'"   }, 	
  {   60.0, -120.0,  -80.0,   0.0, 1, "typeII'"  }, 	
  {  -60.0,  120.0,  -90.0,   0.0, 1, "typeVIa1" }, // 2nd res must be cis-PRO 	
  { -120.0,  120.0,  -60.0,   0.0, 1, "typeVIa2" }, // 2nd res must be cis-PRO 	
  { -135.0,  135.0,  -75.0, 160.0, 1, "typeVIb"  }  // 2nd res must be cis-PRO
};

// CONSTRUCTOR
Action_MakeStructure::Action_MakeStructure() { }

void Action_MakeStructure::Help() {
  mprintf(" <ss type>:<res range> [<ss type>:<res range> ...]\n");
}

// Action_MakeStructure::Init()
Action::RetType Action_MakeStructure::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  SecStructHolder ss_holder;
  secstruct_.clear();
  // Get all keywords
  std::string ss_expr = actionArgs.GetStringNext();
  while ( !ss_expr.empty() ) {
    ArgList ss_arg(ss_expr, ":");
    if (ss_arg.Nargs() < 2) {
      mprinterr("Error: Malformed SS arg, expected <ss type>,<res range>\n");
      return Action::ERR;
    }
    // Find SS type
    int ss = 0;
    for (; ss < (int)NSS; ss++) {
      if (ss_arg[0].compare(SS[ss].name) == 0) { 
        ss_holder.type = (ssType)ss;
        ss_holder.name.assign( SS[ss].name );
        ss_holder.resRange.SetRange(ss_arg[1]);
        // Since user range args start from 1, shift range
        ss_holder.resRange.ShiftBy(-1);
        ss_holder.dihSearch_.Clear();
        ss_holder.dihSearch_.SearchFor(DihedralSearch::PHI);
        ss_holder.dihSearch_.SearchFor(DihedralSearch::PSI);
        secstruct_.push_back( ss_holder );
        break;
      }
    }
    if (ss == (int)NSS) {
      mprinterr("Error: SS type [%s] not found.\n", ss_arg[0].c_str());
      return Action::ERR;
    }
    ss_expr = actionArgs.GetStringNext();
  }
  if (secstruct_.empty()) {
    mprinterr("Error: No SS types defined.\n");
    return Action::ERR;
  }

  return Action::OK;
}

// VisitAtom()
static void VisitAtom( Topology const& topIn, int atm, std::vector<bool>& Visited )
{
  // If this atom has already been visited return
  if (Visited[atm]) return;
  // Mark this atom as visited
  Visited[atm] = true;
  // Visit each atom bonded to this atom
  for (Atom::bond_iterator bondedatom = topIn[atm].bondbegin();
                           bondedatom != topIn[atm].bondend(); ++bondedatom)
    VisitAtom(topIn, *bondedatom, Visited);
}

/// Set up mask of atoms that will move upon rotation.
static AtomMask MovingAtoms(Topology const& topIn, std::vector<bool>& Visited, 
                            int atom0, int atom1) 
{
  Visited.assign( topIn.Natom(), false );
  // Mark atom0 as already visited
  Visited[atom0] = true;
  for (Atom::bond_iterator bndatm = topIn[atom1].bondbegin();
                           bndatm != topIn[atom1].bondend(); ++bndatm)
  {
    if ( *bndatm != atom0 )
      VisitAtom( topIn, *bndatm, Visited );
  }
  // Everything marked T will move.
  AtomMask Rmask;
  for (int maskatom = 0; maskatom < (int)Visited.size(); maskatom++) {
    if (Visited[maskatom])
      Rmask.AddAtom(maskatom);
  }
  return Rmask;
}

// Action_MakeStructure::Setup()
Action::RetType Action_MakeStructure::Setup(Topology* currentParm, Topology** parmAddress) {
  std::vector<bool> Visited( currentParm->Natom(), false );
  // Set up each SS type
  for (std::vector<SecStructHolder>::iterator ss = secstruct_.begin();
                                              ss != secstruct_.end(); ++ss)
  {
    if ((*ss).dihSearch_.FindDihedrals(*currentParm, (*ss).resRange))
      return Action::ERR;
    mprintf("\t[%s] %s %i dih", (*ss).resRange.RangeArg(), (*ss).name.c_str(),
            (*ss).dihSearch_.Ndihedrals());
    // Set up found dihedrals 
    // TODO: Check that # dihedrals is multiple of 2?
    ssType myType = (*ss).type;
    (*ss).Rmasks_.clear();
    (*ss).thetas_.clear();
    if (SS[myType].isTurn) {
      // Require phi/psi residue pairs; must be multiple of 4
      if ( ((*ss).dihSearch_.Ndihedrals() % 4) != 0) {
        mprintf("Error: Assigning turn SS requires residue phi/psi pairs.\n");
        return Action::ERR;
      }
      // Turns require that each pair of residues is consecutive
      for (DihedralSearch::mask_it dih = (*ss).dihSearch_.begin();
                                   dih != (*ss).dihSearch_.end(); ++dih)
      {
        // First has to be phi
        mprintf(" %i:%s", (*dih).ResNum()+1, (*dih).Name().c_str());
        if ( (*dih).Name() != "phi" ) {
          mprinterr("Error: Assigning turn SS requires 1st dihedral be phi.\n");
          return Action::ERR;
        }
        int res1num = (*dih).ResNum();
        (*ss).thetas_.push_back((float)(SS[myType].phi * DEGRAD));
        (*ss).Rmasks_.push_back( MovingAtoms(*currentParm, Visited, (*dih).A1(), (*dih).A2()) );
        ++dih;
        // Second has to be psi and +0
        mprintf("-%i:%s", (*dih).ResNum()+1, (*dih).Name().c_str());
        if ((*dih).Name() != "psi" || ((*dih).ResNum() != res1num)) {
          mprinterr("Error: Assigning turn SS requires 2nd dihedral be psi and consecutive.\n");
          return Action::ERR;
        }
        (*ss).thetas_.push_back((float)(SS[myType].psi * DEGRAD));
        (*ss).Rmasks_.push_back( MovingAtoms(*currentParm, Visited, (*dih).A1(), (*dih).A2()) );
        ++dih;
        // Third has to be phi and +1
        mprintf("-%i:%s", (*dih).ResNum()+1, (*dih).Name().c_str());
        if ((*dih).Name() != "phi" || ((*dih).ResNum() - res1num) != 1) {
          mprinterr("Error: Assigning turn SS requires 3rd dihedral be phi and consecutive.\n");
          return Action::ERR;
        }
        (*ss).thetas_.push_back((float)(SS[myType].phi2 * DEGRAD));
        (*ss).Rmasks_.push_back( MovingAtoms(*currentParm, Visited, (*dih).A1(), (*dih).A2()) );
        ++dih;
        // Fourth has to be phi and +1
        mprintf("-%i:%s", (*dih).ResNum()+1, (*dih).Name().c_str());
        if ((*dih).Name() != "psi" || ((*dih).ResNum() - res1num) != 1) {
          mprinterr("Error: Assigning turn SS requires 4th dihedral be psi and consecutive.\n");
          return Action::ERR;
        }
        (*ss).thetas_.push_back((float)(SS[myType].psi2 * DEGRAD));
        (*ss).Rmasks_.push_back( MovingAtoms(*currentParm, Visited, (*dih).A1(), (*dih).A2()) );
      }
    } else {
      // Not a turn. Assign SS.
      for (DihedralSearch::mask_it dih = (*ss).dihSearch_.begin();
                                   dih != (*ss).dihSearch_.end(); ++dih)
      {
        mprintf(" %i:%s", (*dih).ResNum()+1, (*dih).Name().c_str());
        if ((*dih).Name() == "phi") {
          (*ss).thetas_.push_back((float)(SS[myType].phi * DEGRAD));
          (*ss).Rmasks_.push_back( MovingAtoms(*currentParm, Visited, (*dih).A1(), (*dih).A2()) );
        } else {
          (*ss).thetas_.push_back((float)(SS[myType].psi * DEGRAD));
          (*ss).Rmasks_.push_back( MovingAtoms(*currentParm, Visited, (*dih).A1(), (*dih).A2()) );
        }
      }
    }
    mprintf("\n");
  } // END loop over SS types
  return Action::OK;
}

// Action_MakeStructure::DoAction()
Action::RetType Action_MakeStructure::DoAction(int frameNum, Frame* currentFrame, 
                                               Frame** frameAddress) 
{
  Matrix_3x3 rotationMatrix;
  for (std::vector<SecStructHolder>::iterator ss = secstruct_.begin();
                                              ss != secstruct_.end(); ++ss)
  {
    std::vector<float>::iterator theta = (*ss).thetas_.begin();
    std::vector<AtomMask>::iterator Rmask = (*ss).Rmasks_.begin();
    for (DihedralSearch::mask_it dih = (*ss).dihSearch_.begin();
                                 dih != (*ss).dihSearch_.end(); ++dih, ++theta, ++Rmask)
    {
      double theta_in_radians = (double)*theta;
      // Calculate current value of dihedral
      double torsion = Torsion( currentFrame->XYZ( (*dih).A0() ),
                                currentFrame->XYZ( (*dih).A1() ),
                                currentFrame->XYZ( (*dih).A2() ),
                                currentFrame->XYZ( (*dih).A3() ) );
      // Calculate delta needed to get to theta
      double delta = theta_in_radians - torsion;
      // Set axis of rotation
      Vec3 axisOfRotation = currentFrame->SetAxisOfRotation((*dih).A1(), (*dih).A2());
      // Calculate rotation matrix for delta 
      rotationMatrix.CalcRotationMatrix(axisOfRotation, delta);
//      if (debug_ > 0) {
//        std::string a0name = CurrentParm_->TruncResAtomName( (*dih).A0() );
//        std::string a1name = CurrentParm_->TruncResAtomName( (*dih).A1() );
//        std::string a2name = CurrentParm_->TruncResAtomName( (*dih).A2() );
//        std::string a3name = CurrentParm_->TruncResAtomName( (*dih).A3() );
//        mprintf("\tRotating Dih %s-%s-%s-%s (@%.2f) by %.2f deg to get to %.2f.\n",
//                 a0name.c_str(), a1name.c_str(), a2name.c_str(), a3name.c_str(),
          mprintf("\tRotating Dih %i:%s (%i-%i-%i-%i) (@%.2f) by %.2f deg to get to %.2f.\n",
                  (*dih).ResNum()+1, (*dih).Name().c_str(),
                  (*dih).A0() + 1, (*dih).A1() + 1, (*dih).A2() + 1, (*dih).A3() + 1, 
                  torsion*RADDEG, delta*RADDEG, theta_in_radians*RADDEG);
//      }
    // Rotate around axis
      currentFrame->Rotate(rotationMatrix, *Rmask);
    }
  }
  return Action::OK;
}
