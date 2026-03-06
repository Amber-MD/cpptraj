#include <vector>
#include <algorithm> // std::sort, std::min, std::max, std::find, std::swap
#include <stack>
#include <cmath> // cos
#include "Zmatrix.h"
#include "GenerateConnectivityArrays.h"
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

/** Sort ICs. */
void Zmatrix::sort() {
  std::sort( IC_.begin(), IC_.end() );
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

/** Calculate and add internal coordinate from frame. */
int Zmatrix::AddIC(int ai, int aj, int ak, int al, Frame const& frameIn) {
  addIc(ai, aj, ak, al, frameIn);
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
void Zmatrix::print(Topology const* topIn) const {
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

/** Increment IC indices by given offset */
void Zmatrix::OffsetIcIndices(int offset) {
  for (ICarray::iterator it = IC_.begin(); it != IC_.end(); ++it)
    it->OffsetIndices(offset);
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

/** Calculate an internal coordinate given atom indices
  * and corresponding Cartesian coordinates.
  */
static inline InternalCoords calcIc(int at0, int at1, int at2, int at3,
                   const double* xyz0, const double* xyz1,
                   const double* xyz2, const double* xyz3)
{
  return InternalCoords(at0, at1, at2, at3,
                                sqrt(DIST2_NoImage(xyz0, xyz1)),
                                CalcAngle(xyz0, xyz1, xyz2) * Constants::RADDEG,
                                Torsion(xyz0, xyz1, xyz2, xyz3) * Constants::RADDEG);
}

/** Calculate and add an internal coordinate given atom indices
  * and corresponding Cartesian coordinates.
  */
void Zmatrix::addIc(int at0, int at1, int at2, int at3, Frame const& frameIn) {
  IC_.push_back( calcIc(at0, at1, at2, at3, frameIn.XYZ(at0), frameIn.XYZ(at1), frameIn.XYZ(at2), frameIn.XYZ(at3)) );
}

/** Set seeds as 3 consecutive atoms from mol 0 for which positions are known. */
int Zmatrix::AutoSetSeedsWithPositions(Frame const& frameIn, Topology const& topIn, int ires, Barray const& positionKnown)
{
  return autoSetSeeds_withPositions(frameIn, topIn, topIn.Res(ires).FirstAtom(), topIn.Res(ires).LastAtom(), positionKnown);
}

/** Set seeds as 3 consecutive atoms for which positions are known. */
int Zmatrix::autoSetSeeds_withPositions(Frame const& frameIn, Topology const& topIn, int startAtom, int endAtom, Barray const& positionKnown)
{
  seedAt0_ = InternalCoords::NO_ATOM;
  seedAt1_ = InternalCoords::NO_ATOM;
  seedAt2_ = InternalCoords::NO_ATOM;

  if (positionKnown.empty()) {
    mprinterr("InternalError: Zmatrix::autoSetSeeds_withPositions() called with an empty known position array.\n");
    return 1;
  }
  int numAtoms = endAtom - startAtom;
  if (numAtoms < 1) {
    mprinterr("Internal Error: Zmatrix::autoSetSeeds_withPositions() called with start <= end atom.\n");
    return 1;
  }
  mprintf("DEBUG: autoSetSeeds_withPositions from atoms %s to %s\n", topIn.AtomMaskName(startAtom).c_str(), topIn.AtomMaskName(endAtom-1).c_str());
  // Special cases
  if (numAtoms == 1) {
    seedAt0_ = startAtom;
    seed0Pos_ = Vec3(frameIn.XYZ(seedAt0_));
    return 0;
  } else if (numAtoms == 2) {
    seedAt0_ = startAtom;
    seed0Pos_ = Vec3(frameIn.XYZ(seedAt0_));
    if (topIn[seedAt0_].Nbonds() != 1) {
      mprinterr("Internal Error: Zmatrix::autoSetSeeds_simple(): 2 atoms but no bonds.\n");
      return 1;
    }
    seedAt1_ = topIn[seedAt0_].Bond(0);
    seed1Pos_ = Vec3(frameIn.XYZ(seedAt1_));
    return 0;
  }
  // Loop over atoms in the molecule
  int numS2bonds = -1;
  for (int at = startAtom; at < endAtom; at++) {
//  for (Unit::const_iterator seg = mol.MolUnit().segBegin();
//                            seg != mol.MolUnit().segEnd(); ++seg)
//  {
//    for (int at = seg->Begin(); at != seg->End(); ++at) {
      if (positionKnown[at]) {
        Atom const& AJ = topIn[at];
        if (AJ.Nbonds() > 1) {
          for (int bidx1 = 0; bidx1 < AJ.Nbonds(); bidx1++) {
            for (int bidx2 = bidx1 + 1; bidx2 < AJ.Nbonds(); bidx2++) {
              int bat1 = AJ.Bond(bidx1);
              int bat2 = AJ.Bond(bidx2);
              if (positionKnown[bat1] && positionKnown[bat2]) {
                int s1 = at;
                int s0, s2;
                // b1 - AJ - b2
                Atom const& b1 = topIn[bat1];
                Atom const& b2 = topIn[bat2];
                // The atom with more bonds should be AK (seed 2)
                if (b1.Nbonds() > b2.Nbonds()) {
                  s0 = bat2;
                  s2 = bat1;
                } else {
                  s0 = bat1;
                  s2 = bat2;
                }
                if (numS2bonds == -1 || topIn[s2].Nbonds() > numS2bonds) {
                  seedAt0_ = s0;
                  seedAt1_ = s1;
                  seedAt2_ = s2;
                  numS2bonds = topIn[s2].Nbonds();
                }
              }
            } // END inner loop over bonded atoms
          } // END outer loop over bonded atoms
        } // END AJ bonds > 1
      } // END position of AJ is known
    //} // END loop over segment atoms
  } // END loop over segments

  if (seedAt0_ == InternalCoords::NO_ATOM ||
      seedAt1_ == InternalCoords::NO_ATOM ||
      seedAt2_ == InternalCoords::NO_ATOM)
  {
    mprinterr("Error: No suitable seed atoms with known positions could be found.\n");
    return 1;
  }
  seed0Pos_ = Vec3(frameIn.XYZ(seedAt0_));
  seed1Pos_ = Vec3(frameIn.XYZ(seedAt1_));
  seed2Pos_ = Vec3(frameIn.XYZ(seedAt2_));
  return 0;
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

  if (mol.NumAtoms() < 1) {
    mprinterr("Internal Error: Zmatrix::autoSetSeeds_simple() called with an empty molecule.\n");
    return 1;
  }
  // Handle special cases
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
      addIc(*atI, atJ, atK, atL, frameIn);
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
          addIc(p.ai_, p.aj_, p.ak_, p.al_, frameIn);
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
        addIc(atI, atJ, atK, atL, frameIn);
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

/** Add IC for the given atom. */
int Zmatrix::addInternalCoordForAtom(int iat, Frame const& frameIn, Topology const& topIn)
{
  int maxScore = -1;
  int maxAtj = -1;
  int maxAtk = -1;
  int maxAtl = -1;
  Atom const& atomI = topIn[iat];
  for (Atom::bond_iterator jat = atomI.bondbegin(); jat != atomI.bondend(); ++jat) {
    Atom const& atomJ = topIn[*jat];
    for (Atom::bond_iterator kat = atomJ.bondbegin(); kat != atomJ.bondend(); ++kat) {
      if (*kat != iat) {
        Atom const& atomK = topIn[*kat];
        for (Atom::bond_iterator lat = atomK.bondbegin(); lat != atomK.bondend(); ++lat) {
          if (*lat != *jat && *lat != iat) {
            int nbondScore = topIn[*lat].Nbonds();
            if (maxScore == -1 || nbondScore > maxScore) {
              maxScore = nbondScore;
              maxAtj = *jat;
              maxAtk = *kat;
              maxAtl = *lat;
            }
            //mprintf("DEBUG: Potential IC for %s [ %s - %s - %s ] score= %i\n",
            //        topIn.AtomMaskName(iat).c_str(),
            //        topIn.AtomMaskName(*jat).c_str(),
            //        topIn.AtomMaskName(*kat).c_str(),
            //        topIn.AtomMaskName(*lat).c_str(),
            //        nbondScore);
          }
        }
      }
    }
  }
  if (maxScore == -1) {
    mprintf("Warning: Unable to define IC for atom %s\n", topIn.AtomMaskName(iat).c_str());
  } else {
      addIc( iat, maxAtj, maxAtk, maxAtl, frameIn );
  }
  return 0;
}

// -----------------------------------------------
/** \return the LEaP 'weight' of an atom.
  * Originally used to force the 'heaviest' atoms around a torsion trans to
  * each other. The 'weight' of an atom is defined as its element number,
  * unless the atom is CARBON, then it is 1000, making it the 'heaviest' atom.
  */
static int LeapAtomWeight(Atom const& At)
{
  if ( At.Element() == Atom::CARBON )
    return 1000;
  return At.AtomicNumber();
}

/** Order atoms bonded to the given atom in a manner similar to LEaP's
  * zModelOrderAtoms. In that routine, first atoms were sorted into
  * known position > unknown position. Then the heaviest atom in each
  * subgroup was swapped with the first element of that list. Since at this
  * point we assume all positions are known, we are just shifting the
  * heaviest atom to the front of the list.
  * The ignore atom is the index of the atom this atom is bonded to that
  * forms the torsion we are interested in.
  */
static inline std::vector<int> SortBondedAtomsLikeLeap(Atom const& At, Topology const& topIn, int ignoreAtom)
{
  std::vector<int> out;
  out.reserve( At.Nbonds() );
  // Find the index of the heaviest atom
  int iHighest = 0;
  int iPos = 0;
  for (int idx = 0; idx < At.Nbonds(); idx++) {
    int bat = At.Bond(idx);
    if (bat != ignoreAtom) {
      out.push_back( bat );
      int iWeight = LeapAtomWeight( topIn[bat] );
      if ( iHighest < iWeight ) {
        iHighest = iWeight;
        iPos = (int)out.size()-1;
      }
    }
  }
  // If highest weight atom not already in front, swap it there.
  if (iPos != 0) std::swap( out[0], out[iPos] );

  return out;
}

/** Set up Zmatrix from Cartesian coordinates and topology in the same
  * manner as LEaP's BuildInternalsForContainer/ModelAssignTorsionsAround.
  * Currently assumes all positions are known.
  */
int Zmatrix::GenerateInternals(Frame const& frameIn, Topology const& topIn)
{
  clear();
  // First generate the bond array
  BondArray bonds = GenerateBondArray( topIn.Residues(), topIn.Atoms() );
  // Loop over bonds
  for (BondArray::const_iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd)
  {
    Atom const& A2 = topIn[bnd->A1()];
    Atom const& A3 = topIn[bnd->A2()];
    if (A2.Nbonds() > 1 && A3.Nbonds() > 1) {
      //Residue const& R2 = topIn.Res(A2.ResNum());
      //Residue const& R3 = topIn.Res(A3.ResNum());
      mprintf("Building torsion INTERNALs around: %s - %s\n",
              topIn.LeapName(bnd->A1()).c_str(), topIn.LeapName(bnd->A2()).c_str());
      Iarray sorted_a2 = SortBondedAtomsLikeLeap(A2, topIn, bnd->A2());
      Iarray sorted_a3 = SortBondedAtomsLikeLeap(A3, topIn, bnd->A1());
      mprintf("Orientation around: %s {", *(A2.Name()));
      for (Atom::bond_iterator bat = A2.bondbegin(); bat != A2.bondend(); ++bat) mprintf(" %s", *(topIn[*bat].Name()));
      mprintf("}\n");
      for (Iarray::const_iterator it = sorted_a2.begin(); it != sorted_a2.end(); ++it)
        mprintf("Atom %li: %s\n", it - sorted_a2.begin(), *(topIn[*it].Name()));
      mprintf("Orientation around: %s {", *(A3.Name()));
      for (Atom::bond_iterator bat = A3.bondbegin(); bat != A3.bondend(); ++bat) mprintf(" %s", *(topIn[*bat].Name()));
      mprintf("}\n");
      for (Iarray::const_iterator it = sorted_a3.begin(); it != sorted_a3.end(); ++it)
        mprintf("Atom %li: %s\n", it - sorted_a3.begin(), *(topIn[*it].Name()));
      // Build the torsions
      int aj = bnd->A1();
      int ak = bnd->A2();
      for (Iarray::const_iterator ai = sorted_a2.begin(); ai != sorted_a2.end(); ++ai) {
        for (Iarray::const_iterator al = sorted_a3.begin(); al != sorted_a3.end(); ++al) {
          //double dval = Torsion(frameIn.XYZ(*ai), frameIn.XYZ(aj), frameIn.XYZ(ak), frameIn.XYZ(*al));
          addIc(*ai, aj, ak, *al, frameIn);
          mprintf("++++Torsion INTERNAL: %f to %s - %s - %s - %s\n",
                  IC_.back().Phi(),
                  topIn.LeapName(*ai).c_str(),
                  topIn.LeapName(aj).c_str(),
                  topIn.LeapName(ak).c_str(),
                  topIn.LeapName(*al).c_str());
          addIc(*al, ak, aj, *ai, frameIn);
        }
      }
    }
  }
  return 0;
}

// -----------------------------------------------
/** Set up Zmatrix from Cartesian coordinates and topology.
  * Use torsions based on connectivity to create a complete set of ICs.
  */
int Zmatrix::SetFromFrameAndConnect(Frame const& frameIn, Topology const& topIn) //, int molnum)
{
  clear();

  for (int iat1 = 0; iat1 < topIn.Natom(); iat1++)
  {
    Atom const& At1 = topIn[iat1];
    for (int bidx1 = 0; bidx1 < At1.Nbonds(); bidx1++) {
      int iat2 = At1.Bond(bidx1);
      Atom const& At2 = topIn[iat2];
      for (int bidx2 = 0; bidx2 < At2.Nbonds(); bidx2++) {
        int iat3 = At2.Bond(bidx2);
        if (iat3 != iat1) {
          Atom const& At3 = topIn[iat3];
          for (int bidx3 = 0; bidx3 < At3.Nbonds(); bidx3++) {
            int iat4 = At3.Bond(bidx3);
            if (iat4 != iat2 && iat1 < iat4) {
              //mprintf("DEBUG: DIHEDRAL  %i - %i - %i - %i (%i %i %i %i)\n", iat1+1, iat2+1, iat3+1, iat4+1, iat1*3, iat2*3, iat3*3, iat4*3);
              //out.push_back( DihedralType( iat1, iat2, iat3, iat4, -1 ) );
              addIc(iat1, iat2, iat3, iat4, frameIn);
              addIc(iat4, iat3, iat2, iat1, frameIn); // FIXME should the reverse one be put in?
            }
          }
        }
      }
    }
  }
  if (IC_.empty()) {
    // Either 4-5 atoms in a tetrahedral configuration or else
    // some other strange configuration.
    mprintf("Warning: No ICs created for %s\n", topIn.c_str());
/*
    // Create ICs using 1st 3 heaviest atoms.
    if (topIn.Natom() > 3) {
      typedef std::pair<double,int> MIpair;
      // Used to sort mass descending).
      struct MassCmp {
        inline bool operator()(MIpair const& lhs, MIpair const& rhs) const {
          if (lhs.first == rhs.first)
            return lhs.second < rhs.second;
          else
            return lhs.first > rhs.first;
        }
      };

      std::vector<MIpair> MassIndices;
      MassIndices.reserve( topIn.Natom() );
      for (int iat = 0; iat < topIn.Natom(); iat++)
        MassIndices.push_back( MIpair(topIn[iat].Mass(), iat) );
      // Sort by mass
      std::sort( MassIndices.begin(), MassIndices.end(), MassCmp() );
      mprintf("Sorted:");
      for (std::vector<MIpair>::const_iterator it = MassIndices.begin(); it != MassIndices.end(); ++it)
        mprintf(" %s", topIn.AtomMaskName(it->second).c_str());
      mprintf("\n");
      // Set up ICs
      for (int iat1 = 0; iat1 < topIn.Natom(); iat1++) {
        std::vector<int> iats;
        iats.reserve(3);
        for (std::vector<MIpair>::const_iterator it = MassIndices.begin(); it != MassIndices.end(); ++it)
        {
          if (it->second != iat1)
            iats.push_back( it->second );
          if (iats.size() == 3) break;
        }
        addIc(iat1, iats[0], iats[1], iats[2], frameIn);
      }
    }
*/
  }

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
    int seedErr = 0;
    //if (knownPositions.empty())
      // First seed atom will just be first atom TODO lowest index heavy atom?
      seedErr = autoSetSeeds_simple(frameIn, topIn, currentMol);
    //else
    //  seedErr = autoSetSeeds_withPositions(frameIn, topIn, currentMol, knownPositions);
    if (seedErr != 0) {
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
  // Add basic ICs for seeds
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
          addIc(*bat, seedAt2_, seedAt1_, seedAt0_, frameIn);
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
  // Add ICs for seeds
  /*if (maxnatom > 3) {
    //Atom const& SA0 = topIn[seedAt0_];
    //Atom const& SA1 = topIn[seedAt1_];
    //Atom const& SA2 = topIn[seedAt2_];
    // Seed 0; 0-1-2-X
    addInternalCoordForAtom(seedAt0_, frameIn, topIn);
    addInternalCoordForAtom(seedAt1_, frameIn, topIn);
    addInternalCoordForAtom(seedAt2_, frameIn, topIn);
  }*/
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

// -----------------------------------------------------------------------------
/*
 *  zvZMatrixCalculatePositionFromAngles()
 *
 *  Original Author: Christian Schafmeister (1991)
 *  Adapted by     : Dan Roe (2025) - any mistakes are mine!
 *
 *  Use NEWTON-RAPHSON method for finding the coordinate for the vector(vC)
 *  which is dAngleA from the vector vA (on the X axis) and dAngleB from a
 *  vector (vB) which lies in the XY plane dAngleC from the X axis.
 *
 *  The point is dBond from the origin.
 *
 */
Vec3 Zmatrix::calculatePositionFromAngles( double dAngleA, double dAngleB, 
                                           double dAngleC, double dBond )
{
  static const int MAXNEWTONSTEPS = 20;
  double dCosA = cos(dAngleA);
  double dSinA = sin(dAngleA);
  double dCosB = cos(dAngleB);
  double dCosC = cos(dAngleC);
  double dSinC = sin(dAngleC);

  // The idea is to minimize the function:
  // E = ( DOT(vC,vB) - cos(dAngleB) )^2
  // using NEWTONS method
  // The vector vC is constrained to make the angle
  // dAngleA with vA
  // The vector vC makes the angle (dX) with the XY plane
  // and the parameter that is optimized is dX

  // A reasonable starting point
  double dX = dAngleB;
  int iCount = 0;
  while ( iCount <MAXNEWTONSTEPS ) {

    double dCosX = cos(dX);
    double dSinX = sin(dX);

    double dF1 = -2*dSinA*dSinC*(-dCosB + dCosA*dCosC + dCosX*dSinA*dSinC)*dSinX;

    double dF2 = -2.0*dCosX*dSinA*dSinC* (-dCosB + dCosA*dCosC + dCosX*dSinA*dSinC) + 
                  2.0*dSinA*dSinA*dSinC*dSinC*dSinX*dSinX;

    //    MESSAGE( ( "Iteration %d dF1=%lf  dF2=%lf  dB=%lf\n",
    //                     iCount, dF1, dF2, dX ) );
    if ( fabs(dF1) < Constants::SMALL*10.0 ) break;
    if ( fabs(dF2) < Constants::SMALL ) {
      mprinterr( "Could not optimize! dF1 = %f, dF2 = %f dX = %f steps=%d", dF1, dF2, dX, iCount );
    }
    double dXNew = dX - dF1/dF2;
    if ( fabs(dXNew - dX) < Constants::SMALL ) break;
    dX = dXNew;
    iCount++;
  }
   
# ifdef CPPTRAJ_DEBUG_ZMATRIX
  if ( iCount > MAXNEWTONSTEPS ) 
    mprintf("DEBUG: Exceeded maximum number of Newton Raphson steps: %d\n", MAXNEWTONSTEPS);
# endif
                // Generate new coordinate
  return Vec3( dBond*cos(dAngleA), 
               dBond*sin(dAngleA)*cos(dX),
               dBond*sin(dAngleA)*sin(dX) );
}

/*
 * ZMatrixBondTwoAnglesOrientation()
 * Original Author: Christian Schafmeister (1991)
 * Adapted by     : Dan Roe (2025) - any mistakes are mine!
 *
 *  Build the external coordinate for the atom when the orientation,
 *  a bond length and two angles are supplied.
 *  The orientation is a positive or negative number which specifies
 *  the orientation of the new position.  It is calculated by:
 *    a=crossProduct( vPAtomA-vPCenter, vPAtomB-vPCenter );
 *    orientation = dotProduct( vPPos-vPCenter, a );
 *
 *  \param vPAtomC points to the position of the central atom.
 *  \return Calculated position
 */
Vec3 Zmatrix::PosFromBondTwoAnglesOrientation(
  Vec3 const& vPAtomC, Vec3 const& vPAtomA, Vec3 const& vPAtomB,
  double dBond, double dAngleA, double dAngleB, double dOrient )
{
  static const Vec3 vXAxis(1.0, 0.0, 0.0);
  static const Vec3 vYAxis(0.0, 1.0, 0.0);
  static const Vec3 vZAxis(0.0, 0.0, 1.0);
  //MATRIX          mT, mT1, mT2, mTX, mTY, mTZ, mTT;
  //double          dAngleX, dAngleY, dAngleZ;
  //double          dAngle;
  //VECTOR          vTrans, vTempAC, vTempBC, vTempXZ, vNew, vLab;
  // The procedure for finding the the coordinate is:
  // Translate vAtomC to the origin -> A'-C'
  // Find angle between PROJ((A'-C'),YZ plane) & Y axis
  // Rotate into XZ plane
  // Find angle between (A''-C'') and X axis
  // Rotate onto X axis
  // Find angle between PROJ((B'''-C'''),YZ plane) and Y axis
  // Rotate onto XY plane
  // Calculate coordinates in 3Space
  // Apply the reverse transformation to the new point
  // Actually, all that is done is the elements for the
  // forward transformations are calculated then used
  // to generate an inverse transform matrix.

  Vec3 vTrans = vPAtomC;
  Vec3 vTempAC = vPAtomA - vPAtomC; //vVectorSub( vPAtomA, vPAtomC );
  Vec3 vTempBC = vPAtomB - vPAtomC; //vVectorSub( vPAtomB, vPAtomC );
# ifdef CPPTRAJ_DEBUG_ZMATRIX
  mprintf("AC= %f, %f, %f\n", vTempAC[0], vTempAC[1], vTempAC[2]);
  mprintf("BC= %f, %f, %f\n", vTempBC[0], vTempBC[1], vTempBC[2]);
# endif
  Vec3 vTempXZ = vTempAC;
  vTempXZ[1] = 0.0; //  VectorSetY( &vTempXZ, 0.0 );
  double dAngleY;
  if (vTempXZ.Length() != 0.0) { //  if ( dVectorLen(&vTempXZ) != 0.0 )
    //dAngleY =     dAngleY = dVectorAbsAngle( &vTempXZ, &vXAxis, &vYAxis );
    vTempXZ.Normalize();
    dAngleY = vTempXZ.SignedAngle( vXAxis, vYAxis );
  } else
    dAngleY = 0.0;
# ifdef CPPTRAJ_DEBUG_ZMATRIX
  mprintf("dAngleY= %f\n", dAngleY);
# endif    
  Matrix_3x3 mT;
  mT.RotateAroundY( -dAngleY ); //  MatrixYRotate( mT, -dAngleY );
# ifdef CPPTRAJ_DEBUG_ZMATRIX
  mT.Print("mT (Y)");
# endif
  vTempAC = mT * vTempAC; //  MatrixTimesVector( vTempAC, mT, vTempAC );
  vTempBC = mT * vTempBC; //  MatrixTimesVector( vTempBC, mT, vTempBC );
# ifdef CPPTRAJ_DEBUG_ZMATRIX
  mprintf("Rotated around Y\n" );
  mprintf("New AC= %f, %f, %f\n", vTempAC[0], vTempAC[1], vTempAC[2]);
  mprintf("New BC= %f, %f, %f\n", vTempBC[0], vTempBC[1], vTempBC[2]);
# endif
  vTempAC.Normalize();
  double dAngleZ = vTempAC.SignedAngle( vXAxis, vZAxis ); //dAngleZ = dVectorAbsAngle( &vTempAC, &vXAxis, &vZAxis );
# ifdef CPPTRAJ_DEBUG_ZMATRIX
  mprintf("dAngleZ= %f\n", dAngleZ);
# endif
  mT.RotateAroundZ(-dAngleZ); //  MatrixZRotate( mT, -dAngleZ );
# ifdef CPPTRAJ_DEBUG_ZMATRIX
  mT.Print("mT (Z)");
# endif
  vTempBC = mT * vTempBC; //  MatrixTimesVector( vTempBC, mT, vTempBC );
//#ifdef DEBUG
//    MatrixTimesVector( vTempAC, mT, vTempAC );
//#endif
# ifdef CPPTRAJ_DEBUG_ZMATRIX
  mprintf("Rotated around Z\n" );
  mprintf("New AC= %f, %f, %f\n", vTempAC[0], vTempAC[1], vTempAC[2]);
  mprintf("New BC= %f, %f, %f\n", vTempBC[0], vTempBC[1], vTempBC[2]);
# endif        
  vTempBC[0] = 0; //  VectorSetX( &vTempBC, 0.0 );

  vTempBC.Normalize();
  double dAngleX = vTempBC.SignedAngle( vYAxis, vXAxis ); //   dAngleX = dVectorAbsAngle( &vTempBC, &vYAxis, &vXAxis );
# ifdef CPPTRAJ_DEBUG_ZMATRIX
  mprintf("dAngleX= %f\n", dAngleX);
# endif
                // Build the transformation matrix to convert from
                // lab coordinates to molecule coordinates in mT
//    MatrixXRotate( mTX, dAngleX );
//    MatrixZRotate( mTZ, dAngleZ );
//    MatrixYRotate( mTY, dAngleY );
//    MatrixTranslate( mTT, dVX(&vTrans), dVY(&vTrans), dVZ(&vTrans) );
//    MatrixMultiply( mT1, mTZ, mTX );
//    MatrixMultiply( mT2, mTY, mT1 );
//    MatrixMultiply( mT, mTT, mT2 );
  Matrix_3x3 mT2;
  mT2.CalcRotationMatrix( dAngleX, dAngleY, dAngleZ );
# ifdef CPPTRAJ_DEBUG_ZMATRIX
  mT2.Print("mT2");
# endif    
                // Calculate coordinates of new atom
  double dAngle = CalcAngle(vPAtomA.Dptr(), vPAtomC.Dptr(), vPAtomB.Dptr()); //dAngle = dVectorAtomAngle( vPAtomA, vPAtomC, vPAtomB );
# ifdef CPPTRAJ_DEBUG_ZMATRIX
  mprintf("dAngle= %f\n", dAngle);
# endif
  Vec3 vLab = calculatePositionFromAngles( dAngleA, dAngleB, dAngle, dBond );
# ifdef CPPTRAJ_DEBUG_ZMATRIX
  mprintf("vLab= %f %f %f\n", vLab[0], vLab[1], vLab[2]);
# endif
  if ( dOrient != 0.0 ) {
    vLab[2] = dOrient*vLab[2]; //  VectorSetZ( &vLab, dOrient*dVZ(&vLab) );
  }
                // If there is no chirality defined yet then just 
                // leave it the way it is 
        
  Vec3 vPPos = (mT2 * vLab) + vTrans; //MatrixTimesVector( vNew, mT, vLab );
# ifdef CPPTRAJ_DEBUG_ZMATRIX
  mprintf( "ZMatrix2Angle:  %f,%f,%f\n", vPPos[0], vPPos[1], vPPos[2]);
# endif
  return vPPos;
}

/** \return Position of atom I from positions of atoms J and K,
  *         the corresponding distance I-J (in Ang.) and angle I-J-K
  *         (in radians).
  */
Vec3 Zmatrix::AtomIposition(Vec3 const& posJ, Vec3 const& posK, double rdist, double theta_rad)
{
  // J-K vector
  Vec3 vTempX = posK - posJ;
  Vec3 vTempXZ( vTempX[0], 0.0, vTempX[2] );
  vTempXZ.Normalize();
  //
  static const Vec3 vXAxis(1.0, 0.0, 0.0);
  static const Vec3 vYAxis(0.0, 1.0, 0.0);
  static const Vec3 vZAxis(0.0, 0.0, 1.0);
  double dAngleY = vTempXZ.SignedAngle( vXAxis, vYAxis );
# ifdef CPPTRAJ_DEBUG_ZMATRIX
  mprintf( "Angle around Y=%f\n", dAngleY );
# endif
  Matrix_3x3 mT;
  mT.CalcRotationMatrix(vYAxis, dAngleY);
  vTempX = mT * vTempX;
# ifdef CPPTRAJ_DEBUG_ZMATRIX
  mprintf("Rotated around Y = %f, %f, %f\n", vTempX[0], vTempX[1], vTempX[2]);
# endif
  vTempX.Normalize();
  double dAngleZ = vTempX.SignedAngle( vXAxis, vZAxis );
# ifdef CPPTRAJ_DEBUG_ZMATRIX
  mprintf("Angle around Z=%f\n", dAngleZ );
# endif
  Vec3 vNew( rdist*cos(theta_rad), rdist*sin(theta_rad), 0.0 );
# ifdef CPPTRAJ_DEBUG_ZMATRIX
  mprintf("vNew: %f,%f,%f\n", vNew[0], vNew[1], vNew[2]);
# endif
  mT.CalcRotationMatrix(vZAxis, -dAngleZ);
  vNew = mT * vNew;
  mT.CalcRotationMatrix(vYAxis, -dAngleY);
  vNew = mT * vNew;
# ifdef CPPTRAJ_DEBUG_ZMATRIX
  mprintf("vNew before MatrixTranslate: %f,%f,%f\n", vNew[0], vNew[1], vNew[2]);
  mprintf("vTrans before MatrixTranslate: %f,%f,%f\n", posJ[0], posJ[1], posJ[2]);
# endif
  Vec3 posI = vNew + posJ;

  return posI;
}

/* \return Position of atom I from position of atom J, using the given I-J distance.
 * Atom I will be placed at a distance rdist from posJ along the X axis.
 *
 */
Vec3 Zmatrix::AtomIposition(Vec3 const& posJ, double rdist)
{
  Vec3 vNew( rdist, 0.0, 0.0 );
  vNew = vNew + posJ;
  return vNew;
}

/** \return Position of atom I from positions of atoms J, K, and L,
  *         and corresponding distance (in Ang) and angle/torsion
  *         values (in degrees).
  */
Vec3 Zmatrix::AtomIposition(Vec3 const& posJ, Vec3 const& posK, Vec3 const& posL,
                            double rdist, double theta, double phi)
{
//    double rdist = ic.Dist();
//    double theta = ic.Theta();
//    double phi   = ic.Phi();

    double sinTheta = sin(theta * Constants::DEGRAD);
    double cosTheta = cos(theta * Constants::DEGRAD);
    double sinPhi   = sin(phi   * Constants::DEGRAD);
    double cosPhi   = cos(phi   * Constants::DEGRAD);

    // NOTE: Want -x
    Vec3 xyz( -(rdist * cosTheta),
                rdist * cosPhi * sinTheta,
                rdist * sinPhi * sinTheta );

//    Vec3 posL = Vec3( frameOut.XYZ( ic.AtL()) );
//    Vec3 posK = Vec3( frameOut.XYZ( ic.AtK()) );
//    Vec3 posJ = Vec3( frameOut.XYZ( ic.AtJ()) );

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

/** \return Position of atom I from given internal coordinate and Frame. */
Vec3 Zmatrix::AtomIposition(InternalCoords const& ic, Frame const& frameOut)
{
    double rdist = ic.Dist();
    double theta = ic.Theta();
    double phi   = ic.Phi();

    Vec3 posJ = Vec3( frameOut.XYZ( ic.AtJ()) );
    Vec3 posK = Vec3( frameOut.XYZ( ic.AtK()) );
    Vec3 posL = Vec3( frameOut.XYZ( ic.AtL()) );

    return AtomIposition(posJ, posK, posL, rdist, theta, phi);
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
    mprintf("\tUsing Cartesian seeds.\n");
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
//  unsigned int lowestUnusedIC = 0;
//  for (; lowestUnusedIC < IC_.size(); ++lowestUnusedIC)
//    if (!isUsed[lowestUnusedIC]) break;
//  if (debug_ > 0) mprintf("DEBUG: Lowest unused IC: %u\n", lowestUnusedIC+1);

  // Out of the remaining ICs, count which ones do not have positions set.
  Iarray positionsToSet;
  //unsigned int remainingPositionsToSet = 0;
  for (unsigned int icIdx = 0; icIdx != IC_.size(); ++icIdx) {
    if (!isUsed[icIdx] && !hasPosition[IC_[icIdx].AtI()]) {
      Iarray::const_iterator it = std::find(positionsToSet.begin(), positionsToSet.end(), IC_[icIdx].AtI());
      if (it == positionsToSet.end()) {
        positionsToSet.push_back( IC_[icIdx].AtI() );
        if (debug_ > 1) mprintf("DEBUG:\t\tAtom %i needs its position set.\n", IC_[icIdx].AtI()+1);
      }
    }
  }
  if (debug_ > 0) mprintf("DEBUG: %zu positions to set.\n", positionsToSet.size());

  //while (remainingPositionsToSet > 0 && Nused < IC_.size()) {
  Iarray::const_iterator at_with_unset_pos = positionsToSet.begin();
  while (Nused < IC_.size()) {
    // Get next atom that needs its position set
    for (; at_with_unset_pos != positionsToSet.end(); ++at_with_unset_pos)
      if (!hasPosition[*at_with_unset_pos]) break;
    if (at_with_unset_pos == positionsToSet.end()) {
      if (debug_ > 0) mprintf("DEBUG: No more positions to set.\n");
      break;
    }
    if (debug_ > 1) mprintf("DEBUG: Setting position of atom %i\n", (*at_with_unset_pos) + 1);
    // Get IC that corresponds to this atom
    int icIdx = -1;
    for (unsigned int idx = 0; idx != IC_.size(); ++idx) {
      if (IC_[idx].AtI() == *at_with_unset_pos) {
        // All 3 of the connecting atoms must be set
        if (hasPosition[ IC_[idx].AtJ() ] &&
            hasPosition[ IC_[idx].AtK() ] &&
            hasPosition[ IC_[idx].AtL() ])
        {
          icIdx = (int)idx;
          break;
        }
      }
    }
/*
    // Get the next IC with a position to set
    int icIdx = -1;
    for (unsigned int idx = 0; idx != IC_.size(); ++idx) {
      if (!isUsed[idx] && !hasPosition[IC_[idx].AtI()]) {
        // All 3 of the connecting atoms must be set
        if (hasPosition[ IC_[idx].AtJ() ] &&
            hasPosition[ IC_[idx].AtK() ] &&
            hasPosition[ IC_[idx].AtL() ])
        {
          icIdx = (int)idx;
          break;
        } else {
          mprintf("DEBUG:\t\tIC %u is missing atoms.\n", idx+1);
        }
      }
    }*/
    if (icIdx < 0) {
      if (debug_ > 0) mprintf("DEBUG: Could not find complete IC yet.\n");
      ++at_with_unset_pos;
      if (at_with_unset_pos == positionsToSet.end()) {
        mprinterr("Error: Could not find next IC to use.\n");
        return 1;
      }
      continue; // FIXME
    }
    InternalCoords const& ic = IC_[icIdx];
    if (debug_ > 0) {
      mprintf("DEBUG: Next IC to use is %i : %i %i %i %i r=%g theta=%g phi=%g\n",
              icIdx+1, ic.AtI()+1, ic.AtJ()+1, ic.AtK()+1, ic.AtL()+1,
              ic.Dist(), ic.Theta(), ic.Phi());
      mprintf( "Torsion = %f\n", ic.Phi() );
      mprintf( "Angle   = %f\n", ic.Theta() );
      mprintf( "Bond    = %f\n", ic.Dist() );
    }
    Vec3 posI = AtomIposition(ic, frameOut);

    frameOut.SetXYZ( ic.AtI(), posI );
    hasPosition[ ic.AtI() ] = true;
    //remainingPositionsToSet--;
    MARK(icIdx, isUsed, Nused);
    // Go back to start, new ICs may now be enabled
    at_with_unset_pos = positionsToSet.begin();
  } // END loop over remaining positions
/*
  // Loop over remaining ICs 
  while (Nused < IC_.size()) {
    // Find the next IC that is not yet used.
//    unsigned int idx = lowestUnusedIC;
    unsigned int idx = 0;
    bool findNextIC = true;
    while (findNextIC) {
      while (idx < IC_.size() && isUsed[idx]) {
        mprintf("DEBUG:\t\tIC %u is used\n", idx+1);
        idx++;
      }
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
      } else {
        mprintf("DEBUG:\t\tIC %u is missing atoms.\n", idx+1);
        idx++;
      }
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
*/
  return 0;
}
