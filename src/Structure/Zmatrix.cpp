#include <vector>
#include <algorithm> // std::sort, std::min, std::max
#include <stack>
#include <cmath> // cos
#include <utility> // std::pair
#include "Zmatrix.h"
#include "BuildAtom.h"
#include "Model.h"
#include "../Frame.h"
#include "../CpptrajStdio.h"
#include "../Constants.h"
#include "../DistRoutines.h"
#include "../Topology.h"
#include "../TorsionRoutines.h"

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Zmatrix::Zmatrix() :
  debug_(0),
  icseed0_(InternalCoords::NO_ATOM),
  icseed1_(InternalCoords::NO_ATOM),
  icseed2_(InternalCoords::NO_ATOM),
  seed0Pos_(0.0),
  seed1Pos_(0.0),
  seed2Pos_(0.0),
  seedAt0_(InternalCoords::NO_ATOM),
  seedAt1_(InternalCoords::NO_ATOM),
  seedAt2_(InternalCoords::NO_ATOM)
{}

/** Clear the zmatrix. */
void Zmatrix::clear() {
  IC_.clear();
  icseed0_ = InternalCoords::NO_ATOM;
  icseed1_ = InternalCoords::NO_ATOM;
  icseed2_ = InternalCoords::NO_ATOM;
  seed0Pos_ = Vec3(0.0);
  seed1Pos_ = Vec3(0.0);
  seed2Pos_ = Vec3(0.0);
  seedAt0_ = InternalCoords::NO_ATOM;
  seedAt1_ = InternalCoords::NO_ATOM;
  seedAt2_ = InternalCoords::NO_ATOM;
}

/// Error message for seed already set
static inline int seed_err(int iseed) {
  mprinterr("Internal Error: Internal coord seed %i is already set.\n", iseed);
  return 1;
}

/** Add internal coords. If any of the atoms are not set assume this is one
  * of the 3 seed atoms and determine which one.
  */
int Zmatrix::AddIC(InternalCoords const& ic) {
  if (ic.AtJ() == InternalCoords::NO_ATOM &&
      ic.AtK() == InternalCoords::NO_ATOM &&
      ic.AtL() == InternalCoords::NO_ATOM)
  { // Potential seed0
    if (icseed0_ != InternalCoords::NO_ATOM)
      return seed_err(0);
    icseed0_ = IC_.size();
  } else if (ic.AtK() == InternalCoords::NO_ATOM &&
             ic.AtL() == InternalCoords::NO_ATOM)
  { // Potential seed1
    if (icseed1_ != InternalCoords::NO_ATOM)
      return seed_err(1);
    icseed1_ = IC_.size();
  } else if (ic.AtL() == InternalCoords::NO_ATOM) {
    // Potential seed2
    if (icseed2_ != InternalCoords::NO_ATOM)
      return seed_err(2);
    icseed2_ = IC_.size();
  }

  IC_.push_back( ic );
  return 0;
}

/** Set IC at specified position. */ // TODO check if seed?
void Zmatrix::SetIC(unsigned int idx, InternalCoords const& ic) {
  IC_[idx] = ic;
}

/** \return Array containing indices of ICs with specified atom i. */
std::vector<int> Zmatrix::AtomI_indices(int atomi) const {
  std::vector<int> indices;
  for (unsigned int idx = 0; idx != IC_.size(); idx++) {
    if (IC_[idx].AtI() == atomi)
      indices.push_back( (int)idx );
  }
  return indices;
}

/** Print to stdout */
void Zmatrix::print(Topology* topIn) const {
  mprintf("%zu internal coords.\n", IC_.size());
  mprintf("Seed IC indices    : %i %i %i\n", icseed0_+1, icseed1_+1, icseed2_+1);
  mprintf("Seed Cart. indices : %i %i %i\n", seedAt0_+1, seedAt1_+1, seedAt2_+1);
  seed0Pos_.Print("Seed0");
  seed1Pos_.Print("Seed1");
  seed2Pos_.Print("Seed2");
  if (topIn == 0) {
    mprintf("%-8s %8s %8s %8s %8s %12s %12s %12s\n",
            "#Idx", "AtI", "AtJ", "AtK", "AtL", "Dist", "Theta", "Phi");
    for (ICarray::const_iterator it = IC_.begin(); it != IC_.end(); ++it)
      mprintf("%-8li %8i %8i %8i %8i %12.4f %12.4f %12.4f\n", it - IC_.begin()+1,
              it->AtI()+1, it->AtJ()+1, it->AtK()+1, it->AtL()+1,
              it->Dist(), it->Theta(), it->Phi());
  } else {
    mprintf("%-6s %6s %8s %6s %8s %6s %8s %6s %8s %6s %6s %6s\n",
            "#Idx", "AtI", "NmI", "AtJ", "NmJ", "AtK", "NmK", "AtL", "NmL", "Dist", "Theta", "Phi");
    for (ICarray::const_iterator it = IC_.begin(); it != IC_.end(); ++it)
      mprintf("%-6li %6i %8s %6i %8s %6i %8s %6i %8s %6.2f %6.2f %6.2f\n", it - IC_.begin()+1,
              it->AtI()+1, topIn->AtomMaskName(it->AtI()).c_str(),
              it->AtJ()+1, topIn->AtomMaskName(it->AtJ()).c_str(),
              it->AtK()+1, topIn->AtomMaskName(it->AtK()).c_str(),
              it->AtL()+1, topIn->AtomMaskName(it->AtL()).c_str(),
              it->Dist(), it->Theta(), it->Phi());
  }
}

/** \return True if all Cartesian seeds are set. */
bool Zmatrix::HasCartSeeds() const {
  bool has_cart_seed = (seedAt0_ != InternalCoords::NO_ATOM &&
                        seedAt1_ != InternalCoords::NO_ATOM &&
                        seedAt2_ != InternalCoords::NO_ATOM); 
  return has_cart_seed;
}

/** Set seeds from specified atoms. */
int Zmatrix::SetSeedPositions(Frame const& frameIn, Topology const& topIn, int a1, int a2, int a3)
{
  // a1 must be bonded to a2
  if (!topIn[a1].IsBondedTo(a2)) {
    mprinterr("Error: In topology '%s', seed 0 atom '%s' is not bonded to seed 1 atom '%s'\n",
              topIn.c_str(), topIn.AtomMaskName(a1).c_str(), topIn.AtomMaskName(a2).c_str());
    return 1;
  }
  // a2 must be bonded to a3
  if (!topIn[a2].IsBondedTo(a3)) {
    mprinterr("Error: In topology '%s', seed 1 atom '%s' is not bonded to seed 2 atom '%s'\n",
              topIn.c_str(), topIn.AtomMaskName(a2).c_str(), topIn.AtomMaskName(a3).c_str());
    return 1;
  }
  // Store seed positions
  seedAt0_  = a1;
  seed0Pos_ = Vec3(frameIn.XYZ(a1));
  seedAt1_  = a2;
  seed1Pos_ = Vec3(frameIn.XYZ(a2));
  seedAt2_  = a3;
  seed2Pos_ = Vec3(frameIn.XYZ(a3));
  if (debug_ > 0)
    mprintf("DEBUG: Seed atoms: %s - %s - %s\n",
            topIn.AtomMaskName(a1).c_str(),
            topIn.AtomMaskName(a2).c_str(),
            topIn.AtomMaskName(a3).c_str());

  return 0;
}

/// Mark an IC as used, update used count
static inline void MARK(int i, std::vector<bool>& isUsed, unsigned int& Nused) {
  if (!isUsed[i]) {
    isUsed[i] = true;
    Nused++;
  }
}

/** Calculate and add an internal coordinate given atom indices
  * and corresponding Cartesian coordinates.
  */
void Zmatrix::addIc(int at0, int at1, int at2, int at3,
                   const double* xyz0, const double* xyz1,
                   const double* xyz2, const double* xyz3)
{
  IC_.push_back( InternalCoords(at0, at1, at2, at3,
                                sqrt(DIST2_NoImage(xyz0, xyz1)),
                                CalcAngle(xyz0, xyz1, xyz2) * Constants::RADDEG,
                                Torsion(xyz0, xyz1, xyz2, xyz3) * Constants::RADDEG) );
}

/** Simple automatic setting of seeds for a molecule.
  * seed0 - seed1 - seed 2
  * Prefer that seed1 has only exactly 2 bonds. It cannot have 1.
  */
int Zmatrix::autoSetSeeds_simple(Frame const& frameIn, Topology const& topIn, Molecule const& mol)
{
  seedAt0_ = InternalCoords::NO_ATOM;
  seedAt1_ = InternalCoords::NO_ATOM;
  seedAt2_ = InternalCoords::NO_ATOM;

  // Handle special cases
  if (mol.NumAtoms() < 1) {
    mprinterr("Internal Error: Zmatrix::autoSetSeeds_simple() called with an empty molecule.\n");
    return 1;
  }
  if (mol.NumAtoms() == 1) {
    seedAt0_ = mol.MolUnit().Front();
    seed0Pos_ = Vec3(frameIn.XYZ(seedAt0_));
    return 0;
  } else if (mol.NumAtoms() == 2) {
    seedAt0_ = mol.MolUnit().Front();
    seed0Pos_ = Vec3(frameIn.XYZ(seedAt0_));
    if (topIn[seedAt0_].Nbonds() != 1) {
      mprinterr("Internal Error: Zmatrix::autoSetSeeds_simple(): 2 atoms but no bonds.\n");
      return 1;
    }
    seedAt1_ = topIn[seedAt0_].Bond(0);
    seed1Pos_ = Vec3(frameIn.XYZ(seedAt1_));
    return 0;
  }

  int potentialSeed1 = -1;
  int potentialSeedNbonds = -1;
  for (Unit::const_iterator seg = mol.MolUnit().segBegin();
                            seg != mol.MolUnit().segEnd(); ++seg)
  {
    for (int idx = seg->Begin(); idx != seg->End(); ++idx) {
      if (topIn[idx].Nbonds() > 1) {
        // New seed1 if this seed1 has fewer bonds
        if (potentialSeed1 == -1 || topIn[idx].Nbonds() < potentialSeedNbonds) {
          potentialSeed1 = idx;
          potentialSeedNbonds = topIn[idx].Nbonds();
          if (potentialSeedNbonds == 2) {
            seedAt1_ = idx;
            break;
          }
        }
      }
    } // END loop over segment atoms
    // Exit if a good seed was found
    if (seedAt1_ != InternalCoords::NO_ATOM) break;
  }
  if (seedAt1_ == InternalCoords::NO_ATOM) {
    if (debug_ > 0)
      mprintf("DEBUG: No seed1 with just 2 bonds found. Using seed1 with %i bonds.\n", potentialSeedNbonds);
    seedAt1_ = potentialSeed1;
  }
  if (seedAt1_ == InternalCoords::NO_ATOM) {
    mprinterr("Error: No seed1 could be found.\n");
    return 1;
  }

  // seed0 will be lowest index, seed 1 will be highest
  int lowestidx = -1;
  int highestidx = -1;
  for (Atom::bond_iterator bat = topIn[seedAt1_].bondbegin(); bat != topIn[seedAt1_].bondend(); ++bat)
  {
    if (lowestidx == -1) {
      lowestidx = *bat;
      highestidx = *bat;
    } else {
      lowestidx = std::min(lowestidx, *bat);
      highestidx = std::max(highestidx, *bat);
    }
  }
  seedAt0_ = lowestidx;
  seedAt2_ = highestidx;

  if (debug_ > 0) {
    mprintf("DEBUG: Potential seed 0: %i %s\n", seedAt0_+1, topIn.AtomMaskName(seedAt0_).c_str());
    mprintf("DEBUG: Potential seed 1: %i %s\n", seedAt1_+1, topIn.AtomMaskName(seedAt1_).c_str());
    mprintf("DEBUG: Potential seed 2: %i %s\n", seedAt2_+1, topIn.AtomMaskName(seedAt2_).c_str());
  }
    
  seed0Pos_ = Vec3(frameIn.XYZ(seedAt0_));
  seed1Pos_ = Vec3(frameIn.XYZ(seedAt1_));
  seed2Pos_ = Vec3(frameIn.XYZ(seedAt2_));

  return 0;
}

/// For DEBUG, print integer array
static inline void printIarray(std::vector<int> const& arr, const char* desc, Topology const& topIn) {
  mprintf("DEBUG:\t\t%s:", desc);
  for (std::vector<int>::const_iterator it = arr.begin(); it != arr.end(); ++it)
    //mprintf(" %i", *it);
    mprintf(" %s", topIn.AtomMaskName(*it).c_str());
  mprintf("\n");
}

/// Hold 4 atoms
class PHI {
  public:
    PHI(int ai, int aj, int ak, int al) : ai_(ai), aj_(aj), ak_(ak), al_(al) {}
    int ai_, aj_, ak_, al_;
};

/** Trace molecule starting from L - K - J - ??. */
int Zmatrix::traceMol(int atL0, int atK0, int atJ0,
                      Frame const& frameIn, Topology const& topIn, unsigned int maxnatom,
                      unsigned int& nHasIC, Barray& hasIC)
{
  // Will hold branches that need to be investigated
  std::stack<PHI> Branches;

  int atL = atL0;
  int atK = atK0;
  int atJ = atJ0;

  //int debug_it = 0; // DEBUG
  unsigned int maxStack = 0;

  while (nHasIC < maxnatom) {
    if (debug_ > 1)
      mprintf("\nDEBUG: nHasIC= %8u / %8u : %s - %s - %s -\n", nHasIC, maxnatom,
              topIn.AtomMaskName(atL).c_str(),
              topIn.AtomMaskName(atK).c_str(),
              topIn.AtomMaskName(atJ).c_str());
    // List all atoms bonded to J with the following priority:
    //   1) Atoms with 1 bond.
    //   2) Low index atoms.
    Iarray OneBondAtoms; // TODO allocate outside loop?
    Iarray OtherAtoms;
    for (Atom::bond_iterator bat = topIn[atJ].bondbegin(); bat != topIn[atJ].bondend(); ++bat)
    {
      if (!hasIC[*bat]) {
        if (topIn[*bat].Nbonds() == 1)
          OneBondAtoms.push_back( *bat );
        else if (*bat != atK)
          OtherAtoms.push_back( *bat );
      }
    }
    if (OtherAtoms.size() > 1)
      std::sort( OtherAtoms.begin(), OtherAtoms.end() );
    if (debug_ > 1) {
      printIarray( OneBondAtoms, "OneBondAtoms", topIn );
      printIarray( OtherAtoms, "OtherAtoms", topIn );
    }
    // Create ICs for 1 bond atoms
    for (Iarray::const_iterator atI = OneBondAtoms.begin(); atI != OneBondAtoms.end(); ++atI) {
      addIc(*atI, atJ, atK, atL, frameIn.XYZ(*atI), frameIn.XYZ(atJ), frameIn.XYZ(atK), frameIn.XYZ(atL));
      if (debug_ > 1) {
        mprintf("DEBUG: Added (1 atom) ");
        IC_.back().printIC(topIn);
      }
      MARK(*atI, hasIC, nHasIC);
    }
    // If nothing else, check the stack
    if (OtherAtoms.empty()) {
      // If no branches, done for now.
      if (Branches.empty()) {
        if (debug_ > 1)
          mprintf("DEBUG: No more branches. Exiting traceMol.\n");
        break;
      }
      // Add IC for the branch
      bool has_next = false;
      while (!Branches.empty()) {
        PHI const& p = Branches.top();
        if (debug_ > 1) {
          mprintf("DEBUG:\t\tPopped off stack: %s - %s - %s - %s\n",
                  topIn.AtomMaskName(p.ai_).c_str(),
                  topIn.AtomMaskName(p.aj_).c_str(),
                  topIn.AtomMaskName(p.ak_).c_str(),
                  topIn.AtomMaskName(p.al_).c_str());
        }
        if (!hasIC[p.ai_]) {
          addIc(p.ai_, p.aj_, p.ak_, p.al_, frameIn.XYZ(p.ai_), frameIn.XYZ(p.aj_), frameIn.XYZ(p.ak_), frameIn.XYZ(p.al_));
          if (debug_ > 1) {
            mprintf("DEBUG: Added (stack) ");
            IC_.back().printIC(topIn);
          }
          Branches.pop();
          MARK(p.ai_, hasIC, nHasIC);
          // Designate branch as next.
          atL = p.ak_;
          atK = p.aj_;
          atJ = p.ai_;
          has_next = true;
          break;
        } else {
          if (debug_ > 1) {
            mprintf("DEBUG:\t\t%s already has an IC.\n", topIn.AtomMaskName(p.ai_).c_str());
            std::vector<int> indices = AtomI_indices(p.ai_);
            // TODO check empty or size > 1
            IC_[indices.front()].printIC( topIn );
          }
          Branches.pop();
        }
      } // END while branches remain on stack
      if (!has_next) {
        if (debug_ > 1) mprintf("DEBUG:\t\tNothing left on the stack.\n");
        break;
      } 
    } else {
      int atI = OtherAtoms.front();
      if (debug_ > 1) {
        mprintf("DEBUG:\t\tNext: %s - %s - %s - %s\n",
                topIn.AtomMaskName(atI).c_str(),
                topIn.AtomMaskName(atJ).c_str(),
                topIn.AtomMaskName(atK).c_str(),
                topIn.AtomMaskName(atL).c_str());
      }
      // Add lowest index as IC
      if (!hasIC[atI]) {
        addIc(atI, atJ, atK, atL, frameIn.XYZ(atI), frameIn.XYZ(atJ), frameIn.XYZ(atK), frameIn.XYZ(atL));
        if (debug_ > 1) {
          mprintf("DEBUG: Added (next) ");
          IC_.back().printIC(topIn);
        }
        MARK(atI, hasIC, nHasIC);
      }
      // Place all above lowest index on the stack.
      for (unsigned int ii = 1; ii < OtherAtoms.size(); ii++) {
        Branches.push( PHI(OtherAtoms[ii], atJ, atK, atL) );
        maxStack = std::max(maxStack, (unsigned int)Branches.size());
        if (debug_ > 1) {
          PHI const& p = Branches.top();
          if (debug_ > 1) {
            mprintf("DEBUG:\t\tPlaced on stack: %s - %s - %s - %s (%zu)\n",
                    topIn.AtomMaskName(p.ai_).c_str(),
                    topIn.AtomMaskName(p.aj_).c_str(),
                    topIn.AtomMaskName(p.ak_).c_str(),
                    topIn.AtomMaskName(p.al_).c_str(), Branches.size());
          }
        }
      }
      // Designate lowest index as next
      atL = atK;
      atK = atJ;
      atJ = atI;
    }
    //if (debug_it == 1) break; // DEBUG
    //debug_it++; // DEBUG
  }
  if (debug_ > 0)
    mprintf("DEBUG: Max stack size (%s - %s - %s - ?)= %u\n",
            topIn.AtomMaskName(atL0).c_str(),
            topIn.AtomMaskName(atK0).c_str(),
            topIn.AtomMaskName(atJ0).c_str(),
            maxStack);
  return 0;
}

/** Set up Zmatrix from Cartesian coordinates and topology.
  * This algorithm attempts to "trace" the molecule in a manner that
  * should make internal coordinate assignment more "natural".
  */
int Zmatrix::SetFromFrame_Trace(Frame const& frameIn, Topology const& topIn, int molnum)
{
  if (molnum < 0) {
    mprinterr("Internal Error: Zmatrix::SetFromFrame(): Negative molecule index.\n");
    return 1;
  }
  if (topIn.Nmol() < 1) {
    mprinterr("Internal Error: Zmatrix::SetFromFrame(): No molecules.\n");
    return 1;
  }
  IC_.clear();
  Molecule const& currentMol = topIn.Mol(molnum);
  unsigned int maxnatom = currentMol.NumAtoms();

  // Keep track of which atoms are associated with an internal coordinate
  Barray hasIC( topIn.Natom(), false );
  unsigned int nHasIC = 0;

  // See if we need to assign seed atoms
  if (!HasCartSeeds()) {
    // First seed atom will just be first atom TODO lowest index heavy atom?
    if (autoSetSeeds_simple(frameIn, topIn, currentMol)) {
    //if (autoSetSeeds(frameIn, topIn, maxnatom, currentMol.MolUnit().Front())) {
      mprinterr("Error: Could not automatically determine seed atoms.\n");
      return 1;
    }

  } else {
    // Seed atoms already set
    if (debug_ > 0)
      mprintf("DEBUG: Cartesian Seed indices: %s %s %s\n",
              topIn.AtomMaskName(seedAt0_).c_str(),
              topIn.AtomMaskName(seedAt1_).c_str(),
              topIn.AtomMaskName(seedAt2_).c_str());
  }
  // Add IC seeds
  if (seedAt0_ != InternalCoords::NO_ATOM)
    AddIC( InternalCoords(seedAt0_, InternalCoords::NO_ATOM, InternalCoords::NO_ATOM, InternalCoords::NO_ATOM,
                              0, 0, 0) );
  if (seedAt1_ != InternalCoords::NO_ATOM)
    AddIC( InternalCoords(seedAt1_, seedAt0_, InternalCoords::NO_ATOM, InternalCoords::NO_ATOM,
                              sqrt(DIST2_NoImage(frameIn.XYZ(seedAt1_), frameIn.XYZ(seedAt0_))),
                              0, 0) );
  if (seedAt2_ != InternalCoords::NO_ATOM)
    AddIC( InternalCoords(seedAt2_, seedAt1_, seedAt0_, InternalCoords::NO_ATOM,
                              sqrt(DIST2_NoImage(frameIn.XYZ(seedAt2_), frameIn.XYZ(seedAt1_))),
                              CalcAngle(frameIn.XYZ(seedAt2_), frameIn.XYZ(seedAt1_), frameIn.XYZ(seedAt0_))*Constants::RADDEG,
                              0) );
  // If there are less than 4 atoms we are done
  if (maxnatom < 4) return 0;
  // Seeds are already done
  MARK(seedAt0_, hasIC, nHasIC);
  MARK(seedAt1_, hasIC, nHasIC);
  MARK(seedAt2_, hasIC, nHasIC);

  // Do the remaining atoms.
  // L - K - J - ?
  if (traceMol(seedAt0_, seedAt1_, seedAt2_, frameIn, topIn, maxnatom, nHasIC, hasIC)) return 1;
  // J - K - L - ?
  if (traceMol(seedAt2_, seedAt1_, seedAt0_, frameIn, topIn, maxnatom, nHasIC, hasIC)) return 1;
  if (topIn[seedAt1_].Nbonds() > 2) {
    unsigned int nUnfollowedBranches = 0;
    for (Atom::bond_iterator bat = topIn[seedAt1_].bondbegin(); bat != topIn[seedAt1_].bondend(); ++bat)
      if (!hasIC[*bat]) ++nUnfollowedBranches;
    if (nUnfollowedBranches > 0) {
      mprintf("Warning: Second seed atom %s has more than 2 bonds.\n", topIn.AtomMaskName(seedAt1_).c_str());
      // Potential tetrahedral or more.
      for (Atom::bond_iterator bat = topIn[seedAt1_].bondbegin(); bat != topIn[seedAt1_].bondend(); ++bat) {
        // Ignore seed atoms
        if (*bat != seedAt0_ && *bat != seedAt2_) {
          // Improper based off of seed
          //     I
          //     |
          // L - K - J
          addIc(*bat, seedAt2_, seedAt1_, seedAt0_, frameIn.XYZ(*bat), frameIn.XYZ(seedAt2_), frameIn.XYZ(seedAt1_), frameIn.XYZ(seedAt0_));
          MARK(*bat, hasIC, nHasIC);
          // Follow improper branch if needed: L - K - I - ?
          if (traceMol(seedAt0_, seedAt1_, *bat, frameIn, topIn, maxnatom, nHasIC, hasIC)) return 1;
        }
      }
    }
  }

  if (nHasIC < maxnatom) {
    mprintf("Warning: Not all atoms have an associated internal coordinate.\n");
  }

  return 0;
}

/** Setup Zmatrix from Cartesian coordinates/topology. */
int Zmatrix::SetFromFrame(Frame const& frameIn, Topology const& topIn) {
  return SetFromFrame_Trace(frameIn, topIn, 0);
}

/** Setup Zmatrix from Cartesian coordinates/topology. */
int Zmatrix::SetFromFrame(Frame const& frameIn, Topology const& topIn, int molnum)
{
  return SetFromFrame_Trace(frameIn, topIn, molnum);
}

/** Given two bonded atoms A and B, where B has a depth of at least 2
  * (i.e., it is possible to have B be atom J where we can form J-K-L),
  * set up a complete set of internal coordinates involving A and B in the
  * direction of atom A. This means all internal coordinates with A and B
  * as I and J (should be only 1), as J and K, and as K and L.
  */
int Zmatrix::SetupICsAroundBond(int atA, int atB, Frame const& frameIn, Topology const& topIn,
                                std::vector<bool> const& atomPositionKnown,
                                BuildAtom const& AtomA, BuildAtom const& AtomB)
{
  if (debug_ > 0)
    mprintf("DEBUG: SetupICsAroundBond: atA= %s (%i)  atB= %s (%i)\n",
            topIn.AtomMaskName(atA).c_str(), atA+1,
            topIn.AtomMaskName(atB).c_str(), atB+1);
  IC_.clear();

  Barray hasIC( topIn.Natom(), false );
  unsigned int nHasIC = 0;
  // Mark known atoms as already having IC
  for (std::vector<bool>::const_iterator it = atomPositionKnown.begin();
                                         it != atomPositionKnown.end(); ++it)
  {
    if (*it) MARK( it - atomPositionKnown.begin(), hasIC, nHasIC );
  }

  // First, make sure atom B as a bond depth of at least 2.
  // Choose K and L atoms given atA is I and atB is J.
  Atom const& AJ = topIn[atB];
  typedef std::pair<int,int> Apair;
  std::vector<Apair> KLpairs;
  for (Atom::bond_iterator kat = AJ.bondbegin(); kat != AJ.bondend(); ++kat)
  {
    if (*kat != atA) {
      //mprintf("DEBUG: kat= %s\n", topIn.AtomMaskName(*kat).c_str());
      Atom const& AK = topIn[*kat];
      for (Atom::bond_iterator lat = AK.bondbegin(); lat != AK.bondend(); ++lat)
      {
        if (*lat != atB && *lat != atA) {
          //mprintf("DEBUG: lat= %s\n", topIn.AtomMaskName(*lat).c_str());
          KLpairs.push_back( Apair(*kat, *lat) );
        }
      }
    }
  }
  if (debug_ > 0) {
    for (std::vector<Apair>::const_iterator it = KLpairs.begin();
                                            it != KLpairs.end(); ++it)
      mprintf("DEBUG:\t\tKL pair %s - %s\n", topIn.AtomMaskName(it->first).c_str(),
              topIn.AtomMaskName(it->second).c_str());
  }
  if (KLpairs.empty()) {
    mprinterr("Error: SetFromFrameAroundBond(): Could not find an atom pair bonded to atom %s\n",
              topIn.AtomMaskName(atB).c_str());
    return 1;
  }
  // TODO be smarter about how K and L are selected?
  double maxMass = topIn[KLpairs[0].first].Mass() + topIn[KLpairs[0].second].Mass();
  unsigned int maxIdx = 0;
  for (unsigned int idx = 1; idx < KLpairs.size(); idx++) {
    double sumMass = topIn[KLpairs[idx].first].Mass() + topIn[KLpairs[idx].second].Mass();
    if (sumMass > maxMass) {
      maxMass = sumMass;
      maxIdx = idx;
    }
  }
  int atk0 = KLpairs[maxIdx].first;
  int atl0 = KLpairs[maxIdx].second;
  int modelDebug = 0;
  if (debug_ > 0) {
    mprintf("DEBUG: Chosen KL pair: %s - %s\n",topIn.AtomMaskName(atk0).c_str(),
              topIn.AtomMaskName(atl0).c_str());
    modelDebug = debug_ - 1;
  }
  // ---- I J: Set dist, theta, phi for atA atB K L internal coord ---
  if (debug_ > 0)
    mprintf("DEBUG: IC (i j) %i - %i - %i - %i\n", atA+1, atB+1, atk0+1, atl0+1);
  double newDist = Atom::GetBondLength( topIn[atA].Element(), topIn[atB].Element() );
  if (debug_ > 0) mprintf("DEBUG:\t\tnewDist= %g\n", newDist);
  double newTheta = 0;
  if (Cpptraj::Structure::Model::AssignTheta(newTheta, atA, atB, atk0, topIn, frameIn, atomPositionKnown, modelDebug)) {
    mprinterr("Error: theta (i j) assignment failed.\n");
    return 1;
  }
  if (debug_ > 0) mprintf("DEBUG:\t\tnewTheta = %g\n", newTheta*Constants::RADDEG);
  double newPhi = 0;
  if (Cpptraj::Structure::Model::AssignPhi(newPhi, atA, atB, atk0, atl0, topIn, frameIn,
                                           atomPositionKnown, AtomB, modelDebug))
  {
    mprinterr("Error: phi (i j) assignment failed.\n");
    return 1;
  }
  if (debug_ > 0) mprintf("DEBUG:\t\tnewPhi = %g\n", newPhi*Constants::RADDEG);
  IC_.push_back(InternalCoords( atA, atB, atk0, atl0, newDist, newTheta*Constants::RADDEG, newPhi*Constants::RADDEG ));
  if (debug_ > 0) {
    mprintf("DEBUG: MODEL I J IC: ");
    IC_.back().printIC(topIn);
  }
  MARK( atA, hasIC, nHasIC );
  // ----- J K: Set up ICs for X atA atB K ---------------------------
  Atom const& AJ1 = topIn[atA];
  int ati = -1;
  //Atom const& AK1 = topIn[atB];
  //Atom const& AL1 = topIn[atk0];
  for (Atom::bond_iterator iat = AJ1.bondbegin(); iat != AJ1.bondend(); ++iat)
  {
    if (*iat != atB) {
      if (ati == -1) ati = *iat;
      // Use existing dist
      newDist = sqrt( DIST2_NoImage(frameIn.XYZ(*iat), frameIn.XYZ(atA)) );
      // Set theta for I atA atB
      newTheta = 0;
      if (Cpptraj::Structure::Model::AssignTheta(newTheta, *iat, atA, atB, topIn, frameIn, atomPositionKnown, modelDebug)) {
        mprinterr("Error: theta (j k) assignment failed.\n");
        return 1;
      }
      if (debug_ > 0)
        mprintf("DEBUG:\t\tnewTheta = %g\n", newTheta*Constants::RADDEG);
      // Set phi for I atA atB K
      newPhi = 0;
      if (Cpptraj::Structure::Model::AssignPhi(newPhi, *iat, atA, atB, atk0, topIn, frameIn,
                                               atomPositionKnown, AtomA, modelDebug))
      {
        mprinterr("Error: phi (j k) assignment failed.\n");
        return 1;
      }
      if (debug_ > 0)
        mprintf("DEBUG:\t\tnewPhi = %g\n", newPhi*Constants::RADDEG);
      IC_.push_back(InternalCoords( *iat, atA, atB, atk0, newDist, newTheta*Constants::RADDEG, newPhi*Constants::RADDEG ));
      if (debug_ > 0) {
        mprintf("DEBUG: MODEL J K IC: ");
        IC_.back().printIC(topIn);
      }
      MARK( *iat, hasIC, nHasIC );
      // ----- K L: Set up ICs for X iat atA atB ---------------------
      /*Atom const& AJ2 = topIn[*iat];
      for (Atom::bond_iterator i2at = AJ2.bondbegin(); i2at != AJ2.bondend(); ++i2at)
      {
        if (*i2at != atA && *i2at != atB && !hasIC[*i2at]) {
          // Use existing dist and theta
          newDist = sqrt( DIST2_NoImage(frameIn.XYZ(*i2at), frameIn.XYZ(*iat)) );
          newTheta = CalcAngle( frameIn.XYZ(*i2at), frameIn.XYZ(*iat), frameIn.XYZ(atA) );
          // Set phi for X iat atA atB
          newPhi = 0;
          if (Cpptraj::Structure::Model::AssignPhi(newPhi, *i2at, *iat, atA, atB, topIn, frameIn, atomPositionKnown, atomChirality)) {
            mprinterr("Error: phi (k l) assignment failed.\n");
            return 1;
          }
          mprintf("DEBUG:\t\tnewPhi = %g\n", newPhi*Constants::RADDEG);
          IC_.push_back(InternalCoords( *i2at, *iat, atA, atB, newDist, newTheta*Constants::RADDEG, newPhi*Constants::RADDEG ));
          mprintf("DEBUG: MODEL K L IC: ");
          IC_.back().printIC(topIn);
          MARK( *i2at, hasIC, nHasIC );
          // Trace from atA *iat *i2at outwards
          ToTrace.push_back(atA);
          ToTrace.push_back(*iat);
          ToTrace.push_back(*i2at);
          //if (traceMol(atA, *iat, *i2at, frameIn, topIn, topIn.Natom(), nHasIC, hasIC)) return 1;
        }
      } */
    }
  }
  // Handle remaining atoms.
  if (AJ1.Nbonds() > 1) {
    if (AJ1.Nbonds() == 2) {
      if (debug_ > 0) mprintf("DEBUG: 2 bonds to %s.\n", topIn.AtomMaskName(atA).c_str());
      if (traceMol(atB, atA, ati, frameIn, topIn, topIn.Natom(), nHasIC, hasIC))
        return 1;
    } else {
      // 3 or more bonds
      std::vector<int> const& priority = AtomA.Priority();
      int at0 = -1;
      int at1 = -1;
      std::vector<int> remainingAtoms;
      if (debug_ > 0) mprintf("DEBUG: %i bonds to %s\n", AJ1.Nbonds(), topIn.AtomMaskName(atA).c_str());
      for (std::vector<int>::const_iterator it = priority.begin(); it != priority.end(); ++it) {
        if (*it != atB) {
          if (debug_ > 0) mprintf("DEBUG:\t\t%s\n", topIn.AtomMaskName(*it).c_str());
          if (at0 == -1)
            at0 = *it;
          else if (at1 == -1)
            at1 = *it;
          else
            remainingAtoms.push_back( *it );
        }
      }
      // at0 atA at1
      if (traceMol(at1, atA, at0, frameIn, topIn, topIn.Natom(), nHasIC, hasIC))
        return 1;
      // at1 atA, at0
      if (traceMol(at0, atA, at1, frameIn, topIn, topIn.Natom(), nHasIC, hasIC))
        return 1;
      // Remaining atoms.
      for (std::vector<int>::const_iterator it = remainingAtoms.begin(); it != remainingAtoms.end(); ++it) {
        if (traceMol(at0, atA, *it, frameIn, topIn, topIn.Natom(), nHasIC, hasIC))
          return 1;
      }
    }
  }

  return 0;
}

// -----------------------------------------------------------------------------
/** \return Position of atom I from given internal coordinate. */
Vec3 Zmatrix::AtomIposition(InternalCoords const& ic, Frame const& frameOut)
{
    double rdist = ic.Dist();
    double theta = ic.Theta();
    double phi   = ic.Phi();

    double sinTheta = sin(theta * Constants::DEGRAD);
    double cosTheta = cos(theta * Constants::DEGRAD);
    double sinPhi   = sin(phi   * Constants::DEGRAD);
    double cosPhi   = cos(phi   * Constants::DEGRAD);

    // NOTE: Want -x
    Vec3 xyz( -(rdist * cosTheta),
                rdist * cosPhi * sinTheta,
                rdist * sinPhi * sinTheta );

    Vec3 posL = Vec3( frameOut.XYZ( ic.AtL()) );
    Vec3 posK = Vec3( frameOut.XYZ( ic.AtK()) );
    Vec3 posJ = Vec3( frameOut.XYZ( ic.AtJ()) );

    Vec3 LK = posK - posL;
    Vec3 KJ = posJ - posK;
    KJ.Normalize();
    Vec3 Norm = LK.Cross(KJ);
    Norm.Normalize();
    Vec3 NxKJ = Norm.Cross(KJ);

    Matrix_3x3 Rot( KJ[0], NxKJ[0], Norm[0],
                    KJ[1], NxKJ[1], Norm[1],
                    KJ[2], NxKJ[2], Norm[2] );

    Vec3 posI = (Rot * xyz) + posJ;

  return posI;
}

/** Set Cartesian coordinates in Frame from internal coordinates.
  * Assume none of the positions in frameOut can be used initially.
  */
int Zmatrix::SetToFrame(Frame& frameOut) const {
  // Track which atoms have Cartesian coords set
  Barray hasPosition( frameOut.Natom(), false );
  return SetToFrame(frameOut, hasPosition);
}

/** Set Cartesian coordinates in Frame from internal coordinates.
  * Any atom with hasPosition set to true is considered a "good"
  * position that other ICs can use.
  * The procedure used here is from:
  * Parsons et al., "Practical Conversion from Torsion Space to 
  * Cartesian Space for In Silico Protein Synthesis",
  * J Comput Chem 26: 1063–1068, 2005.
  */
int Zmatrix::SetToFrame(Frame& frameOut, Barray& hasPosition) const {
  // If any seed positions are defined, set them now
  if (seedAt0_ != InternalCoords::NO_ATOM) {
    frameOut.SetXYZ( seedAt0_, seed0Pos_ );
    hasPosition[ seedAt0_ ] = true;
  }
  if (seedAt1_ != InternalCoords::NO_ATOM) {
    frameOut.SetXYZ( seedAt1_, seed1Pos_ );
    hasPosition[ seedAt1_ ] = true;
  }
  if (seedAt2_ != InternalCoords::NO_ATOM) {
    frameOut.SetXYZ( seedAt2_, seed2Pos_ );
    hasPosition[ seedAt2_ ] = true;
  }
  // Track which ICs are used
  Barray isUsed( IC_.size(), false );
  unsigned int Nused = 0;
  if (HasCartSeeds()) {
    if (icseed0_ != InternalCoords::NO_ATOM) MARK(icseed0_, isUsed, Nused);
    if (icseed1_ != InternalCoords::NO_ATOM) MARK(icseed1_, isUsed, Nused);
    if (icseed2_ != InternalCoords::NO_ATOM) MARK(icseed2_, isUsed, Nused);
  } else {
    // Set positions of atoms from internal coordinate seeds. TODO check for clashes with seedAtX?
    if (icseed0_ != InternalCoords::NO_ATOM) {
      // First seed IC atom
      frameOut.SetXYZ(IC_[icseed0_].AtI(), Vec3(0.0));
      hasPosition[IC_[icseed0_].AtI()] = true;
      MARK(icseed0_, isUsed, Nused);
      // Set position of the second atom.
      if (icseed1_ != InternalCoords::NO_ATOM) {
        if (IC_[icseed1_].AtJ() != IC_[icseed0_].AtI()) {
          mprinterr("Internal Error: Atom j of seed 1 is not Atom i of seed 0.\n");
          return 1;
        }
        double r1 = IC_[icseed1_].Dist();
        frameOut.SetXYZ(IC_[icseed1_].AtI(), Vec3(r1, 0, 0));
        hasPosition[IC_[icseed1_].AtI()] = true;
        MARK(icseed1_, isUsed, Nused);
        // Set position of the third atom
        if (icseed2_ != InternalCoords::NO_ATOM) {
          if (IC_[icseed2_].AtJ() != IC_[icseed1_].AtI()) {
            mprinterr("Internal Error: Atom j of seed 2 is not Atom i of seed 1.\n");
            return 1;
          }
          if (IC_[icseed2_].AtK() != IC_[icseed0_].AtI()) {
            mprinterr("Internal Error: Atom k of seed 2 is not Atom i of seed 0.\n");
            return 1;
          }
          double r2 = IC_[icseed2_].Dist();
          double theta = IC_[icseed2_].Theta();

          double x = r2 * cos(180.0 - theta) * Constants::DEGRAD;
          double y = r2 * cos(180.0 - theta) * Constants::DEGRAD;

          frameOut.SetXYZ( IC_[icseed2_].AtI(), Vec3(r1 + x, y, 0) );
          hasPosition[IC_[icseed2_].AtI()] = true;
          MARK(icseed2_, isUsed, Nused);
        } // END seed atom 2
      } // END seed atom 1
    } // END seed atom 0
  } // END Does not have Cart. seeds
  // Find the lowest unused IC
  unsigned int lowestUnusedIC = 0;
  for (; lowestUnusedIC < IC_.size(); ++lowestUnusedIC)
    if (!isUsed[lowestUnusedIC]) break;
  if (debug_ > 0) mprintf("DEBUG: Lowest unused IC: %u\n", lowestUnusedIC+1);

  // Loop over remaining ICs 
  while (Nused < IC_.size()) {
    // Find the next IC that is not yet used.
    unsigned int idx = lowestUnusedIC;
    bool findNextIC = true;
    while (findNextIC) {
      while (idx < IC_.size() && isUsed[idx]) idx++;
      if (idx >= IC_.size()) {
        mprinterr("Error: Could not find next IC to use.\n");
        return 1;
      }
      // All 3 of the connecting atoms must be set
      if (hasPosition[ IC_[idx].AtJ() ] &&
          hasPosition[ IC_[idx].AtK() ] &&
          hasPosition[ IC_[idx].AtL() ])
      {
        findNextIC = false;
      } else
        idx++;
    } // END loop finding next atom to set
    if (debug_ > 0) mprintf("DEBUG: Next IC to use is %u\n", idx+1);

    InternalCoords const& ic = IC_[idx];
    Vec3 posI = AtomIposition(ic, frameOut);

    frameOut.SetXYZ( ic.AtI(), posI );
    hasPosition[ ic.AtI() ] = true;
    MARK(idx, isUsed, Nused);

    // Next lowest unused IC
    for (; lowestUnusedIC < IC_.size(); ++lowestUnusedIC)
      if (!isUsed[lowestUnusedIC]) break;


    //break; // DEBUG
  } // END loop over internal coordinates

  return 0;
}
