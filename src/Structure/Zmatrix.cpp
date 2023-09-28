#include <vector>
#include <set>
#include "Zmatrix.h"
#include "../Frame.h"
#include "../CpptrajStdio.h"
#include "../Constants.h"
#include "../DistRoutines.h"
#include "../Topology.h"
#include "../TorsionRoutines.h"
#include <cmath> // cos

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Zmatrix::Zmatrix() :
  debug_(0),
  icseed0_(InternalCoords::NO_ATOM),
  icseed1_(InternalCoords::NO_ATOM),
  icseed2_(InternalCoords::NO_ATOM),
  seedAt0_(InternalCoords::NO_ATOM),
  seedAt1_(InternalCoords::NO_ATOM),
  seedAt2_(InternalCoords::NO_ATOM)
{}

/** Add internal coords */
void Zmatrix::AddIC(InternalCoords const& ic) {
  IC_.push_back( ic );
}

/** Add internal coords as a seed. */
int Zmatrix::AddICseed(InternalCoords const& ic) {
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
}

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

/// For bonded atoms, hold atom index, atom number, and # bonds.
class AtnumNbonds {
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
static std::set<AtnumNbonds> getBondedAtomPriorities(int seed0, Topology const& topIn, int ignoreIdx) {
  std::set<AtnumNbonds> bondedAtoms;
  for (Atom::bond_iterator bat = topIn[seed0].bondbegin();
                           bat != topIn[seed0].bondend(); ++bat)
  {
    if (*bat != ignoreIdx)
      bondedAtoms.insert( AtnumNbonds(*bat, topIn[*bat]) );
  }
  for (std::set<AtnumNbonds>::const_iterator it = bondedAtoms.begin(); it != bondedAtoms.end(); ++it)
    mprintf("DEBUG: Atom %s bonded atom %s idx=%i priority= %i nbonds=%i\n",
            topIn.AtomMaskName(seed0).c_str(),
            topIn.AtomMaskName(it->Idx()).c_str(),
            it->Idx(), it->Priority(), it->Nbonds());
  return bondedAtoms;
}

/// \return Index of atom at front of set
static inline int FrontIdx(std::set<AtnumNbonds> const& in) {
  std::set<AtnumNbonds>::const_iterator it = in.begin();
  return it->Idx();
}

/// \return Index of first set atom, or atom at front of set
/*static inline int FirstOrFrontIdx(std::set<AtnumNbonds> const& in, std::vector<bool> const& isSet)
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

/// Set j k and l indices for given atom i
/*static inline int SetJKL(int ai, Topology const& topIn, std::vector<bool> const& isSet)
{
  int aj, ak, al;

  std::set<AtnumNbonds> bondedAtoms = getBondedAtomPriorities(ai, topIn, -1);
  if (bondedAtoms.empty()) return 1;
  int firstIdx = FirstOrFrontIdx( bondedAtoms, isSet );

  return 0;
}*/

/** \return True if all IC seeds are set. */
bool Zmatrix::HasICSeeds() const {
  bool has_ic_seed = (icseed0_ != InternalCoords::NO_ATOM &&
                      icseed1_ != InternalCoords::NO_ATOM &&
                      icseed2_ != InternalCoords::NO_ATOM   );
  return has_ic_seed;
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
  mprintf("DEBUG: Seed atoms: %s - %s - %s\n",
          topIn.AtomMaskName(a1).c_str(),
          topIn.AtomMaskName(a2).c_str(),
          topIn.AtomMaskName(a3).c_str());

  return 0;
}

//const int Zmatrix::DUMMY0 = -2;
//const int Zmatrix::DUMMY1 = -3;
//const int Zmatrix::DUMMY2 = -4;

/// Mark an IC as used, update used count
static inline void MARK(int i, std::vector<bool>& isUsed, unsigned int& Nused) {
  if (!isUsed[i]) {
    isUsed[i] = true;
    Nused++;
  }
}

void Zmatrix::addIc(int at0, int at1, int at2, int at3,
                   const double* xyz0, const double* xyz1,
                   const double* xyz2, const double* xyz3)
{
  IC_.push_back( InternalCoords(at0, at1, at2, at3,
                                sqrt(DIST2_NoImage(xyz0, xyz1)),
                                CalcAngle(xyz0, xyz1, xyz2) * Constants::RADDEG,
                                Torsion(xyz0, xyz1, xyz2, xyz3) * Constants::RADDEG) );
}

/** Setup Zmatrix from Cartesian coordinates/topology. */
int Zmatrix::SetFromFrame(Frame const& frameIn, Topology const& topIn, int molnum)
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
  // Flatten the molecule array
  std::vector<int> atomIndices;
  atomIndices.reserve( currentMol.NumAtoms() );
  for (Unit::const_iterator seg = currentMol.MolUnit().segBegin();
                            seg != currentMol.MolUnit().segEnd(); ++seg)
    for (int idx = seg->Begin(); idx != seg->End(); ++idx)
      atomIndices.push_back( idx );
  unsigned int maxnatom = atomIndices.size();
  // Keep track of which atoms are associated with an internal coordinate
  std::vector<bool> hasIC( topIn.Natom(), false );
  unsigned int nHasIC = 0;
  // See if we need to assign seed atoms
  if (!HasCartSeeds()) {
    // First seed atom will just be first atom TODO lowest index heavy atom?
    int at0 = atomIndices.front();
    if (maxnatom < 2) {
      seedAt0_ = at0;
      seed0Pos_ = Vec3(frameIn.XYZ(seedAt0_));
      return 0;
    }
    // Choose second seed atom as bonded atom with lowest index. Prefer heavy
    // atoms and atoms with more than 1 bond.
    std::set<AtnumNbonds> bondedTo0 = getBondedAtomPriorities(at0, topIn, -1);
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
    std::set<AtnumNbonds> bondedTo1 = getBondedAtomPriorities(at1, topIn, at0);
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
    MARK(seedAt0_, hasIC, nHasIC);
    MARK(seedAt1_, hasIC, nHasIC);
    MARK(seedAt2_, hasIC, nHasIC);
    seed0Pos_ = Vec3(frameIn.XYZ(seedAt0_));
    seed1Pos_ = Vec3(frameIn.XYZ(seedAt1_));
    seed2Pos_ = Vec3(frameIn.XYZ(seedAt2_));
  } else {
    // Seed atoms already set
    mprintf("DEBUG: Cartesian Seed indices: %s %s %s\n",
            topIn.AtomMaskName(seedAt0_).c_str(),
            topIn.AtomMaskName(seedAt1_).c_str(),
            topIn.AtomMaskName(seedAt2_).c_str());
  }


/*
  // Do the remaining atoms
  unsigned int maxAtom = (unsigned int)topIn.Natom();
  std::vector<bool> isSet( maxAtom, false );
  isSet[seed0_] = true;
  isSet[seed1_] = true;
  isSet[seed2_] = true;
  unsigned int Nset = 3;

  // Find the lowest unset atom
  unsigned int lowestUnsetAtom = 0;
  for (; lowestUnsetAtom < maxAtom; ++lowestUnsetAtom)
    if (!isSet[lowestUnsetAtom]) break;
  //if (debug_ > 0)
    mprintf("DEBUG: Lowest unset atom: %u\n", lowestUnsetAtom+1);

  // Loop over remaining atoms
  while (Nset < maxAtom) {
   // Find the next atom that is not yet set.
    unsigned int idx = lowestUnsetAtom;
    bool findNextAtom = true;
    while (findNextAtom) {
      while (idx < maxAtom && isSet[idx]) idx++;
      if (idx >= maxAtom) {
        mprinterr("Error: Could not find next atom to set.\n");
        return 1;
      }
      SetJKL(idx, topIn, isSet);
      break; // DEBUG
      // All 3 of the connecting atoms must be set
      //if (isSet[ IC_[idx].AtJ() ] &&
      //    isSet[ IC_[idx].AtK() ] &&
      //    isSet[ IC_[idx].AtL() ])
      //{
      //  findNextAtom = false;
      //}
    } // END loop finding next atom to set
    //if (debug_ > 0)
    mprintf("DEBUG: Next atom to set is %u\n", idx+1);
    break; //DEBUG
  } // END loop over remaining atoms
   */

  return 0;
}

/** Set Cartesian coordinates in Frame from internal coordinates.
  * The procedure used here is from:
  * Parsons et al., "Practical Conversion from Torsion Space to 
  * Cartesian Space for In Silico Protein Synthesis",
  * J Comput Chem 26: 1063–1068, 2005.
  */
int Zmatrix::SetToFrame(Frame& frameOut) const {
  // Track which atoms have Cartesian coords set
  std::vector<bool> hasPosition( frameOut.Natom(), false );
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
  std::vector<bool> isUsed( IC_.size(), false );
  unsigned int Nused = 0;
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
      }
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
