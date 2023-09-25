#include <vector>
#include "Zmatrix.h"
#include "../Frame.h"
#include "../CpptrajStdio.h"
#include "../Constants.h"
#include <cmath> // cos

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Zmatrix::Zmatrix() :
  seed0_(InternalCoords::NO_ATOM),
  seed1_(InternalCoords::NO_ATOM),
  seed2_(InternalCoords::NO_ATOM)
{}

/** Add internal coords */
void Zmatrix::AddIC(InternalCoords const& ic) {
  IC_.push_back( ic );
}

/** Add internal coords as a seed. */
int Zmatrix::AddICseed(InternalCoords const& ic) {
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
  return 0;
}

/** Print to stdout */
void Zmatrix::print() const {
  mprintf("%zu internal coords.\n", IC_.size());
  mprintf("Seed atoms: %i %i %i\n", seed0_+1, seed1_+1, seed2_+1);
  for (ICarray::const_iterator it = IC_.begin(); it != IC_.end(); ++it)
    mprintf("\t%8li %8i %8i %8i %12.4f %12.4f %12.4f\n", it - IC_.begin() + 1,
            it->AtJ()+1, it->AtK()+1, it->AtL()+1,
            it->Dist(), it->Theta(), it->Phi());
}

/** Setup Zmatrix from Cartesian coordinates/topology. */
int Zmatrix::SetFromFrame(Frame const& frameIn, Topology const& topIn)
{
  IC_.clear();
  // First seed is first atom. No bonds, angles, or torsions.
  IC_.push_back( InternalCoords() );
  seed0_ = 0;

  // Choose second seed as bonded atom with lowest index. Prefer heavy atoms

  return 0;
}

static inline void atomIsSet(int i, std::vector<bool>& isSet, unsigned int& Nset) {
  if (!isSet[i]) {
    isSet[i] = true;
    Nset++;
  }
}

/** Set Cartesian coordinates in Frame from internal coordinates. */
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
    //IC_[seed0_].ZeroXYZ();
    frameOut.SetXYZ(seed0_, Vec3(0.0));
    atomIsSet(seed0_, isSet, Nset);
    // Set position of the second atom.
    if (seed1_ != InternalCoords::NO_ATOM) {
      if (IC_[seed1_].AtJ() != seed0_) {
        mprinterr("Internal Error: Atom j of seed 1 is not seed 0.\n");
        return 1;
      }
      double r1 = IC_[seed1_].Dist();
      //IC_[seed1_].SetXYZ( Vec3(r1, 0, 0) );
      frameOut.SetXYZ(seed1_, Vec3(r1, 0, 0));
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

        //IC_[seed2].SetXYZ( Vec3(r1 + x, y, 0) );
        frameOut.SetXYZ( seed2_, Vec3(r1 + x, y, 0) );
        atomIsSet(seed2_, isSet, Nset);
      } // END seed atom 2
    } // END seed atom 1
  } // END seed atom 0

  // Find the lowest unset atom
  unsigned int lowestUnsetAtom = 0;
  for (; lowestUnsetAtom < IC_.size(); ++lowestUnsetAtom)
    if (!isSet[lowestUnsetAtom]) break;
  mprintf("Lowest unset atom: %u\n", lowestUnsetAtom+1);

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
    mprintf("DEBUG: Next atom to set is %u\n", idx+1);
    break; // DEBUG
  } // END loop over internal coordinates



  return 0;
}
