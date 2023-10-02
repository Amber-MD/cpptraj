#include <vector>
#include <algorithm> // std::sort, std::min, std::max
#include <stack>
#include <cmath> // cos
#include "Zmatrix.h"
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

/** Add internal coords as a IC seed. This is intended for use with systems
  * that have dummy atoms, such as those from Amber Prep files.
  */
/*int Zmatrix::AddICseed(InternalCoords const& ic) {
  if (icseed0_ == InternalCoords::NO_ATOM)
    icseed0_ = IC_.size();
  else if (icseed1_ == InternalCoords::NO_ATOM)
    icseed1_ = IC_.size();
  else if (icseed2_ == InternalCoords::NO_ATOM)
    icseed2_ = IC_.size();
  else {
    mprinterr("Error: Too many seed ICs.\n");
    return 1;
  }
  IC_.push_back( ic );
  return 0;
}*/

/** Print to stdout */
void Zmatrix::print() const {
  mprintf("%zu internal coords.\n", IC_.size());
  mprintf("Seed IC indices    : %i %i %i\n", icseed0_+1, icseed1_+1, icseed2_+1);
  mprintf("Seed Cart. indices : %i %i %i\n", seedAt0_+1, seedAt1_+1, seedAt2_+1);
  seed0Pos_.Print("Seed0");
  seed1Pos_.Print("Seed1");
  seed2Pos_.Print("Seed2");
  mprintf("%-8s %8s %8s %8s %8s %12s %12s %12s\n",
          "#Idx", "AtI", "AtJ", "AtK", "AtL", "Dist", "Theta", "Phi");
  for (ICarray::const_iterator it = IC_.begin(); it != IC_.end(); ++it)
    mprintf("%-8li %8i %8i %8i %8i %12.4f %12.4f %12.4f\n", it - IC_.begin()+1,
            it->AtI()+1, it->AtJ()+1, it->AtK()+1, it->AtL()+1,
            it->Dist(), it->Theta(), it->Phi());
}

// -------------------------------------
/// For bonded atoms, hold atom index, atom number, and # bonds.
/** Used to determine atom priority when looking for torsions to
  * base ICs off of.
  */
/*class AtnumNbonds {
  public:
    /// CONSTRUCTOR
    AtnumNbonds() : idx_(-1), priority_(-1), nbonds_(-1) {}
    /// CONSTRUCTOR
    AtnumNbonds(int idx, Atom const& atm) :
      idx_(idx), priority_(mainChainPriority(atm)), nbonds_(atm.Nbonds()) {}
    /// Sort by priority, # bonds, atom index
    bool operator<(const AtnumNbonds& rhs) const {
      if (idx_ < 0) return false;
      if (priority_ == rhs.priority_) {
        if (nbonds_ == rhs.nbonds_) {
          return (idx_ < rhs.idx_);
        } else {
          return (nbonds_ < rhs.nbonds_);
        }
      } else {
        return (priority_ < rhs.priority_);
      }
    }
    int Idx() const { return idx_; }
    int Priority() const { return priority_; }
    int Nbonds() const { return nbonds_; }
  private:
    /// Set priority based on how likely this is to be a main chain atom.
    static int mainChainPriority(Atom const& atm) {
      switch(atm.Element()) {
        case Atom::CARBON     : return 0; // highest priority
        case Atom::NITROGEN   :
        case Atom::BORON      : 
        case Atom::PHOSPHORUS : return 1;
        case Atom::OXYGEN     :
        case Atom::SULFUR     : return 2;
        // These atoms form only 1 bond and have lowest priority
        case Atom::HYDROGEN   :
        case Atom::FLUORINE   :
        case Atom::CHLORINE   :
        case Atom::BROMINE    :
        case Atom::LITHIUM    : return 4;
        default               : return 3; // 1 bond priority - 1
      }
      return 5; // should never get here
    }

    int idx_;      ///< Atom index
    int priority_; ///< Likelyhood of begin a main chain atom (0 most likely)
    int nbonds_;   ///< Number of bonded atoms
};

/// Get bonded atom priorities
static std::set<AtnumNbonds> getBondedAtomPriorities(int seed0, Topology const& topIn, int ignoreIdx0, int ignoreIdx1) {
  std::set<AtnumNbonds> bondedAtoms;
  for (Atom::bond_iterator bat = topIn[seed0].bondbegin();
                           bat != topIn[seed0].bondend(); ++bat)
  {
    if (*bat != ignoreIdx0 && *bat != ignoreIdx1)
      bondedAtoms.insert( AtnumNbonds(*bat, topIn[*bat]) );
  }
  for (std::set<AtnumNbonds>::const_iterator it = bondedAtoms.begin(); it != bondedAtoms.end(); ++it)
    mprintf("DEBUG: Atom %s bonded atom %s idx=%i priority= %i nbonds=%i\n",
            topIn.AtomMaskName(seed0).c_str(),
            topIn.AtomMaskName(it->Idx()).c_str(),
            it->Idx(), it->Priority(), it->Nbonds());
  return bondedAtoms;
}

/// \return Index of atom at front of the set (highest priority).
static inline int FrontIdx(std::set<AtnumNbonds> const& in) {
  std::set<AtnumNbonds>::const_iterator it = in.begin();
  return it->Idx();
}

/// \return Index of first atom that isSet (has coords), or atom at front of the set (highest priority).
static inline int FirstOrFrontIdx(std::set<AtnumNbonds> const& in, std::vector<bool> const& isSet)
{
  if (in.empty()) return -1;
  int firstIdx = FrontIdx( in );
  int firstSetIdx = -1;
  if (in.size() == 1)
    firstSetIdx = firstIdx;
  else {
    for (std::set<AtnumNbonds>::const_iterator it = in.begin();
                                               it != in.end(); ++it)
    {
      if (isSet[it->Idx()]) {
        firstSetIdx = it->Idx();
        break;
      }
    }
  }
  mprintf("DEBUG: First idx= %i  First set idx= %i\n", firstIdx, firstSetIdx);
  if (firstSetIdx > -1)
    return firstSetIdx;
  return firstIdx;
}*/
// -------------------------------------

/** \return True if all IC seeds are set. */
/*bool Zmatrix::HasICSeeds() const {
  bool has_ic_seed = (icseed0_ != InternalCoords::NO_ATOM &&
                      icseed1_ != InternalCoords::NO_ATOM &&
                      icseed2_ != InternalCoords::NO_ATOM   );
  return has_ic_seed;
}*/

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

/** Given a first seed, automatically determine remaining 2 seeds.
  * Set all seed indices and positions.
  */
/*
int Zmatrix::autoSetSeeds(Frame const& frameIn, Topology const& topIn, unsigned int maxnatom, int firstSeed)
{
  int at0 = firstSeed;
  if (maxnatom < 2) {
    seedAt0_ = at0;
    seed0Pos_ = Vec3(frameIn.XYZ(seedAt0_));
    return 0;
  }
  // Choose second seed atom as bonded atom with lowest index. Prefer heavy
  // atoms and atoms with more than 1 bond.
  std::set<AtnumNbonds> bondedTo0 = getBondedAtomPriorities(at0, topIn, -1, -1);
  if (bondedTo0.empty()) {
    mprinterr("Internal Error: Zmatrix::SetFromFrame(): could not get second seed atom.\n");
    return 1;
  }
  int at1 = FrontIdx( bondedTo0 );
  if (maxnatom < 3) {
    seedAt1_ = at1;
    seed1Pos_ = Vec3(frameIn.XYZ(seedAt1_));
    return 0;
  }
  // The third seed atom will either be bonded to seed 0 or seed 1.
  AtnumNbonds potential2From0, potential2From1;
  std::set<AtnumNbonds> bondedTo1 = getBondedAtomPriorities(at1, topIn, at0, -1);
  if (!bondedTo1.empty()) {
    std::set<AtnumNbonds>::const_iterator it = bondedTo1.begin();
    potential2From1 = *it;
  }
  if (bondedTo0.size() >= 2) {
    std::set<AtnumNbonds>::const_iterator it = bondedTo0.begin();
    ++it;
    potential2From0 = *it;
  }
  if (potential2From0 < potential2From1) {
    mprintf("DEBUG: 2 - 0 - 1\n");
    seedAt0_ = potential2From0.Idx();
    seedAt1_ = at0;
    seedAt2_ = at1;
    // D0 D1 D2 A2
    //addIc(at2, DUMMY2, DUMMY1, DUMMY0, frameIn.XYZ(at2), seed2Pos_.Dptr(), seed1Pos_.Dptr(), seed0Pos_.Dptr());
    // D1 D2 A2 A0
    //addIc(at0, at2, DUMMY2, DUMMY1, frameIn.XYZ(at0), frameIn.XYZ(at2), seed2Pos_.Dptr(), seed1Pos_.Dptr());
    // D2 A2 A0 A1
    //addIc(at1, at0, at2, DUMMY2, frameIn.XYZ(at1), frameIn.XYZ(at0), frameIn.XYZ(at2), seed2Pos_.Dptr());
  } else {
    mprintf("DEBUG: 0 - 1 - 2\n");
    seedAt0_ = at0;
    seedAt1_ = at1;
    seedAt2_ = potential2From1.Idx();
    // D0 D1 D2 A0
    //addIc(at0, DUMMY2, DUMMY1, DUMMY0, frameIn.XYZ(at0), seed2Pos_.Dptr(), seed1Pos_.Dptr(), seed0Pos_.Dptr());
    // D1 D2 A0 A1
    //addIc(at1, at0, DUMMY2, DUMMY1, frameIn.XYZ(at1), frameIn.XYZ(at0), seed2Pos_.Dptr(), seed1Pos_.Dptr());
    // D2 A0 A1 A2
    //addIc(at2, at1, at0, DUMMY2, frameIn.XYZ(at2), frameIn.XYZ(at1), frameIn.XYZ(at0), seed2Pos_.Dptr());
  }
  mprintf("DEBUG: Seed atoms: %s - %s - %s\n",
          topIn.AtomMaskName(seedAt0_).c_str(), 
          topIn.AtomMaskName(seedAt1_).c_str(), 
          topIn.AtomMaskName(seedAt2_).c_str());
  seed0Pos_ = Vec3(frameIn.XYZ(seedAt0_));
  seed1Pos_ = Vec3(frameIn.XYZ(seedAt1_));
  seed2Pos_ = Vec3(frameIn.XYZ(seedAt2_));
  return 0;
}*/

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
        else
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
      PHI const& p = Branches.top();
      if (debug_ > 1) {
        mprintf("DEBUG:\t\tPopped off stack: %s - %s - %s - %s\n",
                topIn.AtomMaskName(p.al_).c_str(),
                topIn.AtomMaskName(p.ak_).c_str(),
                topIn.AtomMaskName(p.aj_).c_str(),
                topIn.AtomMaskName(p.ai_).c_str());
      }
      addIc(p.ai_, p.aj_, p.ak_, p.al_, frameIn.XYZ(p.ai_), frameIn.XYZ(p.aj_), frameIn.XYZ(p.ak_), frameIn.XYZ(p.al_));
      Branches.pop();
      MARK(p.ai_, hasIC, nHasIC);
      // Designate branch as next.
      atL = p.ak_;
      atK = p.aj_;
      atJ = p.ai_;
    } else {
      int atI = OtherAtoms.front();
      if (debug_ > 1) {
        mprintf("DEBUG:\t\tNext: %s - %s - %s - %s\n",
                topIn.AtomMaskName(atL).c_str(),
                topIn.AtomMaskName(atK).c_str(),
                topIn.AtomMaskName(atJ).c_str(),
                topIn.AtomMaskName(atI).c_str());
      }
      // Add lowest index as IC
      addIc(atI, atJ, atK, atL, frameIn.XYZ(atI), frameIn.XYZ(atJ), frameIn.XYZ(atK), frameIn.XYZ(atL));
      MARK(atI, hasIC, nHasIC);
      // Place all above lowest index on the stack.
      for (unsigned int ii = 1; ii < OtherAtoms.size(); ii++) {
        Branches.push( PHI(OtherAtoms[ii], atJ, atK, atL) );
        maxStack = std::max(maxStack, (unsigned int)Branches.size());
        if (debug_ > 1) {
          PHI const& p = Branches.top();
          mprintf("DEBUG:\t\tPlaced on stack: %s - %s - %s - %s (%zu)\n",
                  topIn.AtomMaskName(p.al_).c_str(),
                  topIn.AtomMaskName(p.ak_).c_str(),
                  topIn.AtomMaskName(p.aj_).c_str(),
                  topIn.AtomMaskName(p.ai_).c_str(), Branches.size());
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
/*
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
  // Flatten the molecule array
  Iarray atomIndices;
  atomIndices.reserve( currentMol.NumAtoms() );
  for (Unit::const_iterator seg = currentMol.MolUnit().segBegin();
                            seg != currentMol.MolUnit().segEnd(); ++seg)
    for (int idx = seg->Begin(); idx != seg->End(); ++idx)
      atomIndices.push_back( idx );
  unsigned int maxnatom = atomIndices.size();

  // Keep track of which atoms are associated with an internal coordinate
  Barray hasIC( topIn.Natom(), false );
  unsigned int nHasIC = 0;

  // See if we need to assign seed atoms
  if (!HasCartSeeds()) {
    // First seed atom will just be first atom TODO lowest index heavy atom?
    if (autoSetSeeds(frameIn, topIn, maxnatom, atomIndices.front())) {
      mprinterr("Error: Could not automatically determine seed atoms.\n");
      return 1;
    }

  } else {
    // Seed atoms already set
    mprintf("DEBUG: Cartesian Seed indices: %s %s %s\n",
            topIn.AtomMaskName(seedAt0_).c_str(),
            topIn.AtomMaskName(seedAt1_).c_str(),
            topIn.AtomMaskName(seedAt2_).c_str());
  }
  // If there are less than 4 atoms we are done
  if (maxnatom < 4) return 0;
  // Seeds are already done
  MARK(seedAt0_, hasIC, nHasIC);
  MARK(seedAt1_, hasIC, nHasIC);
  MARK(seedAt2_, hasIC, nHasIC);
  // Do the remaining atoms.
  // Find the lowest unset atom.
  unsigned int lowestUnsetIdx = 0;
  for (; lowestUnsetIdx < maxnatom; ++lowestUnsetIdx)
    if (!hasIC[atomIndices[lowestUnsetIdx]]) break;
  mprintf("DEBUG: Lowest unset atom %i at index %u\n", atomIndices[lowestUnsetIdx]+1, lowestUnsetIdx);

  // Loop over remaining atoms
  while (nHasIC < maxnatom) {
    // Find the next atom that is not yet set.
    unsigned int idx = lowestUnsetIdx;
    bool findNextAtom = true;
    int ai = -1;
    int aj = -1;
    int ak = -1;
    int al = -1;
    while (findNextAtom) {
      while (idx < maxnatom && hasIC[atomIndices[idx]]) idx++;
      if (idx >= maxnatom) {
        mprinterr("Error: Could not find next atom to set.\n");
        return 1;
      }
      ai = atomIndices[idx];
      mprintf("DEBUG:\tAttempting to set atom i %i\n", ai+1);
      // Determine j k and l indices. Prefer atoms that have already been set.
      std::set<AtnumNbonds> bondedAtomsI = getBondedAtomPriorities(ai, topIn, -1, -1);
      aj = FirstOrFrontIdx( bondedAtomsI, hasIC );
      mprintf("DEBUG:\t\tAtom j = %i\n", aj+1);
      std::set<AtnumNbonds> bondedAtomsJ = getBondedAtomPriorities(aj, topIn, ai, -1); // TODO reuse bondedAtomsI?
      ak = FirstOrFrontIdx( bondedAtomsJ, hasIC );
      mprintf("DEBUG:\t\tAtom k = %i\n", ak+1);
      // FIXME check for cycle here?
      std::set<AtnumNbonds> bondedAtomsK = getBondedAtomPriorities(ak, topIn, aj, ai); // TODO reuse bondedAtomsI?
      al = FirstOrFrontIdx( bondedAtomsK, hasIC );
      mprintf("DEBUG:\t\tAtom l = %i\n", al+1);
      if (aj < 0 || ak < 0 || al < 0) {
        //mprinterr("Internal Error: Could not find torsion for atom %i %s\n", ai+1, topIn.AtomMaskName(ai).c_str());
        //return 1;
        mprintf("Warning: Could not find torsion for atom %i %s\n", ai+1, topIn.AtomMaskName(ai).c_str());
        mprintf("Warning: Creating pseudo torsion using seed atoms.\n");
        aj = seedAt2_;
        ak = seedAt1_;
        al = seedAt0_;
        // TODO mark such torsions as problematic in InternalCoords?
      }
      // All 3 of the connecting atoms must be set
      if (hasIC[aj] && hasIC[ak] && hasIC[al]) {
        findNextAtom = false;
      } else
        idx++;
    } // END loop over findNextAtom
    // Create IC for i j k l
    addIc(ai, aj, ak, al, frameIn.XYZ(ai), frameIn.XYZ(aj), frameIn.XYZ(ak), frameIn.XYZ(al));
    MARK(ai, hasIC, nHasIC);

    // Next lowest unset atom
    for (; lowestUnsetIdx < maxnatom; ++lowestUnsetIdx)
      if (!hasIC[atomIndices[lowestUnsetIdx]]) break;

    //break; //DEBUG
  } // END loop over remaining atoms

  return 0;*/
}

/** Set Cartesian coordinates in Frame from internal coordinates.
  * The procedure used here is from:
  * Parsons et al., "Practical Conversion from Torsion Space to 
  * Cartesian Space for In Silico Protein Synthesis",
  * J Comput Chem 26: 1063–1068, 2005.
  */
int Zmatrix::SetToFrame(Frame& frameOut) const {
  // Track which atoms have Cartesian coords set
  Barray hasPosition( frameOut.Natom(), false );
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
