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
  seed0_(InternalCoords::NO_ATOM),
  seed1_(InternalCoords::NO_ATOM),
  seed2_(InternalCoords::NO_ATOM)
{}

/** Add internal coords */
void Zmatrix::AddIC(InternalCoords const& ic, int topIdx) {
  IC_.push_back( ic );
  topIndices_.push_back( topIdx );
}

/** Add internal coords as a seed. */
int Zmatrix::AddICseed(InternalCoords const& ic, int topIdx) {
  if (seed0_ == InternalCoords::NO_ATOM)
    seed0_ = IC_.size();
  else if (seed1_ == InternalCoords::NO_ATOM)
    seed1_ = IC_.size();
  else if (seed2_ == InternalCoords::NO_ATOM)
    seed2_ = IC_.size();
  else {
    mprinterr("Error: Too many seed atoms.\n");
    return 1;
  }
  IC_.push_back( ic );
  topIndices_.push_back( topIdx );
  return 0;
}

/** Print to stdout */
void Zmatrix::print() const {
  mprintf("%zu internal coords.\n", IC_.size());
  mprintf("Seed indices: %i %i %i\n", seed0_+1, seed1_+1, seed2_+1);
  for (ICarray::const_iterator it = IC_.begin(); it != IC_.end(); ++it)
    mprintf("\t%8li %8i %8i %8i %12.4f %12.4f %12.4f\n", it - IC_.begin() + 1,
            it->AtJ()+1, it->AtK()+1, it->AtL()+1,
            it->Dist(), it->Theta(), it->Phi());
}

/// For bonded atoms, hold atom index, atom number, and # bonds.
class AtnumNbonds {
  public:
    /// CONSTRUCTOR
    AtnumNbonds(int idx, Atom const& atm) :
      idx_(idx), priority_(mainChainPriority(atm)), nbonds_(atm.Nbonds()) {}
    /// Sort by priority, # bonds, atom index
    bool operator<(const AtnumNbonds& rhs) const {
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
}

/// Set j k and l indices for given atom i
static inline int SetJKL(int ai, Topology const& topIn, std::vector<bool> const& isSet)
{
  int aj, ak, al;

  std::set<AtnumNbonds> bondedAtoms = getBondedAtomPriorities(ai, topIn, -1);
  if (bondedAtoms.empty()) return 1;
  int firstIdx = FirstOrFrontIdx( bondedAtoms, isSet );

  return 0;
}

/** \return True if no seeds are set. */
bool Zmatrix::NoSeeds() const {
  return (seed0_ == InternalCoords::NO_ATOM ||
          seed1_ == InternalCoords::NO_ATOM ||
          seed2_ == InternalCoords::NO_ATOM   );
}

/** Setup Zmatrix from Cartesian coordinates/topology. */
int Zmatrix::SetFromFrame(Frame const& frameIn, Topology const& topIn, int molnum)
{
  IC_.clear();
  topIndices_.clear();
  // See if we need to assign seed atoms
  if (NoSeeds()) {
    mprinterr("Internal Error: Automatic seed generation not yet implemented.\n");
/*    mprintf("DEBUG: Generating dummy seed atoms.\n");
    seed0_ = InternalCoords::NO_ATOM;
    seed1_ = InternalCoords::NO_ATOM;
    seed2_ = InternalCoords::NO_ATOM;
    // Seed 0 is at the origin
    IC_.push_back( InternalCoords() );
    seed0_ = 0;
    topIndices_.push_back( -1 );
    seed0Pos_ = Vec3(0.0);
    // Seed 1 is at X=1
    IC_.push_back( InternalCoords(0, InternalCoords::NO_ATOM, InternalCoords::NO_ATOM,
                                  1.0, 0, 0) );
    seed1_ = 1;
    topIndices_.push_back( -1 );
    seed1Pos_ = Vec3(1.0, 0.0, 0.0);
    // Seed 2 is at X=Y=1
    IC_.push_back( InternalCoords(1, 0, InternalCoords::NO_ATOM,
                                  1.0, 90.0, 0) );
    seed2_ = 2;
    topIndices_.push_back( -1 );
    seed2Pos_ = Vec3(1.0, 1.0, 0.0);*/
  } else {
    mprintf("DEBUG: Seed atoms: %s %s %s\n",
            topIn.AtomMaskName(topIndices_[seed0_]).c_str(),
            topIn.AtomMaskName(topIndices_[seed1_]).c_str(),
            topIn.AtomMaskName(topIndices_[seed2_]).c_str());
  }

/*
  // First seed is first atom. No bonds, angles, or torsions. TODO should be lowest heavy atom?
  IC_.push_back( InternalCoords() );
  seed0_ = 0;

  // Choose second seed as bonded atom with lowest index. Prefer heavy atoms
  // and atoms with more than 1 bond.
  std::set<AtnumNbonds> bondedAtoms = getBondedAtomPriorities(seed0_, topIn, -1);
  if (!bondedAtoms.empty()) {
    seed1_ = FrontIdx(bondedAtoms);
    IC_.push_back( InternalCoords(seed0_, InternalCoords::NO_ATOM, InternalCoords::NO_ATOM,
                                  sqrt(DIST2_NoImage( frameIn.XYZ(seed1_), frameIn.XYZ(seed0_) )),
                                  0, 0) );
    // Choose third seed, ignoring first seed.
    bondedAtoms = getBondedAtomPriorities(seed1_, topIn, seed0_);
    if (!bondedAtoms.empty()) {
      seed2_ = FrontIdx(bondedAtoms);
      IC_.push_back( InternalCoords(seed1_, seed0_, InternalCoords::NO_ATOM,
                                    sqrt(DIST2_NoImage( frameIn.XYZ(seed2_), frameIn.XYZ(seed1_) )),
                                    CalcAngle( frameIn.XYZ(seed2_),
                                               frameIn.XYZ(seed1_),
                                               frameIn.XYZ(seed0_) ) * Constants::RADDEG,
                                    0) );
    }
  }
  // If less than 4 atoms, all done.
  if (topIn.Natom() < 4) return 0;

  if (seed0_ == InternalCoords::NO_ATOM ||
      seed1_ == InternalCoords::NO_ATOM ||
      seed2_ == InternalCoords::NO_ATOM)
  {
    mprinterr("Internal Error: Zmatrix::SetFromFrame(): Not enough seeds.\n");
    return 1;
  }*/


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
   

  return 0;
}

static inline void atomIsSet(int i, std::vector<bool>& isSet, unsigned int& Nset) {
  if (!isSet[i]) {
    isSet[i] = true;
    Nset++;
  }
}

/** Set Cartesian coordinates in Frame from internal coordinates.
  * The procedure used here is from:
  * Parsons et al., "Practical Conversion from Torsion Space to 
  * Cartesian Space for In Silico Protein Synthesis",
  * J Comput Chem 26: 1063–1068, 2005.
  */
int Zmatrix::SetToFrame(Frame& frameOut) const {
  if ((unsigned int)frameOut.Natom() != IC_.size()) {
    mprinterr("Internal Error: Output frame size (%i) != # internal coords (%zu)\n",
              frameOut.Natom(), IC_.size());
    return 1;
  }

  std::vector<bool> isSet( IC_.size(), false );
  unsigned int Nset = 0;
  // Set position of the first atom.
  if (seed0_ != InternalCoords::NO_ATOM) {
    frameOut.SetXYZ(topIndices_[seed0_], Vec3(0.0));
    atomIsSet(seed0_, isSet, Nset);
    // Set position of the second atom.
    if (seed1_ != InternalCoords::NO_ATOM) {
      if (IC_[seed1_].AtJ() != seed0_) {
        mprinterr("Internal Error: Atom j of seed 1 is not seed 0.\n");
        return 1;
      }
      double r1 = IC_[seed1_].Dist();
      frameOut.SetXYZ(topIndices_[seed1_], Vec3(r1, 0, 0));
      atomIsSet(seed1_, isSet, Nset);
      // Set position of the third atom
      if (seed2_ != InternalCoords::NO_ATOM) {
        if (IC_[seed2_].AtJ() != seed1_) {
          mprinterr("Internal Error: Atom j of seed 2 is not seed 1.\n");
          return 1;
        }
        if (IC_[seed2_].AtK() != seed0_) {
          mprinterr("Internal Error: Atom k of seed 2 is not seed 0.\n");
          return 1;
        }
        double r2 = IC_[seed2_].Dist();
        double theta = IC_[seed2_].Theta();

        double x = r2 * cos(180.0 - theta) * Constants::DEGRAD;
        double y = r2 * cos(180.0 - theta) * Constants::DEGRAD;

        frameOut.SetXYZ( topIndices_[seed2_], Vec3(r1 + x, y, 0) );
        atomIsSet(seed2_, isSet, Nset);
      } // END seed atom 2
    } // END seed atom 1
  } // END seed atom 0

  // Find the lowest unset atom
  unsigned int lowestUnsetAtom = 0;
  for (; lowestUnsetAtom < IC_.size(); ++lowestUnsetAtom)
    if (!isSet[lowestUnsetAtom]) break;
  if (debug_ > 0) mprintf("DEBUG: Lowest unset atom: %u\n", lowestUnsetAtom+1);

  // Loop over remaining atoms
  while (Nset < IC_.size()) {
    // Find the next atom that is not yet set.
    unsigned int idx = lowestUnsetAtom;
    bool findNextAtom = true;
    while (findNextAtom) {
      while (idx < IC_.size() && isSet[idx]) idx++;
      if (idx >= IC_.size()) {
        mprinterr("Error: Could not find next atom to set.\n");
        return 1;
      }
      // All 3 of the connecting atoms must be set
      if (isSet[ IC_[idx].AtJ() ] &&
          isSet[ IC_[idx].AtK() ] &&
          isSet[ IC_[idx].AtL() ])
      {
        findNextAtom = false;
      }
    } // END loop finding next atom to set
    if (debug_ > 0) mprintf("DEBUG: Next atom to set is %u\n", idx+1);

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

    frameOut.SetXYZ( topIndices_[idx], posI );
    atomIsSet(idx, isSet, Nset);

    // Next lowest unset atom
    for (; lowestUnsetAtom < IC_.size(); ++lowestUnsetAtom)
      if (!isSet[lowestUnsetAtom]) break;


    //break; // DEBUG
  } // END loop over internal coordinates



  return 0;
}
