#include "Chirality.h"
#include "../Atom.h"
#include "../Constants.h"
#include "../CpptrajStdio.h"
#include "../Frame.h"
#include "../Topology.h"
#include "../TorsionRoutines.h"
#include "../Vec3.h"
#include <cmath> // fabs
#include <algorithm> // sort

/// \return Total priority (i.e. sum of atomic numbers) of atoms bonded to given atom.
static int totalPriority(Topology const& topIn, int atnum, int rnum,
                         int depth, int tgtdepth, std::vector<bool>& Visited)
{
  if (Visited[atnum] || depth == tgtdepth) return 0;
  Visited[atnum] = true;
  if (topIn[atnum].ResNum() != rnum) {
    if (depth == 0)
      // Atom in another residue bonded to chiral atom (e.g. O to C1)
      return topIn[atnum].AtomicNumber();
    else
      // Atom in another residue bonded to atom bonded to chiral atom
      return 0;
  }
  int sum = 0;
  Atom const& atom = topIn[atnum];
  for (Atom::bond_iterator bat = atom.bondbegin(); bat != atom.bondend(); ++bat)
    sum += topIn[*bat].AtomicNumber() + totalPriority(topIn, *bat, rnum, depth+1, tgtdepth, Visited);
  return sum;
}

/// Used to determine priority of moieties bonded to an atom
class priority_element {
  public:
    /// CONSTRUCT from atom number and initial priority
    priority_element(int a, int p1) : atnum_(a), priority1_(p1), priority2_(-1) {}
    /// Set priority 2
    void SetPriority2(int p2) { priority2_ = p2; }
    /// \return Atom number
    int AtNum() const { return atnum_; }
    /// \return Priority 1
    int Priority1() const { return priority1_; }
    /// \return Priority 2
    int Priority2() const { return priority2_; }
    /// Sort on priority 1, then priority 2
    bool operator<(const priority_element& rhs) const {
      if (*this != rhs) {
        if (priority1_ == rhs.priority1_) {
          return (priority2_ > rhs.priority2_);
        } else {
          return (priority1_ > rhs.priority1_);
        }
      } else
        return false;
    }
    /// \return true if priorities are identical
    bool operator==(const priority_element& rhs) const {
      return (priority1_ == rhs.priority1_) && (priority2_ == rhs.priority2_);
    }
    /// \return true if priorities are not equal
    bool operator!=(const priority_element& rhs) const {
      if (priority2_ != rhs.priority2_ ||
          priority1_ != rhs.priority1_) return true;
      return false;
    }
  private:
    int atnum_;
    int priority1_;
    int priority2_;
};


/** Given an atom that is a chiral center, attempt to calculate a
  * torsion that will help determine R vs S. Priorities will be 
  * assigned to bonded atoms as 1, 2, 3, and optionally 4. The
  * torsion will then be calculated as
  *   1-2-3-0
  * where 0 is the chiral center. Negative is S, positive is R.
  */
Cpptraj::Structure::ChiralType
  Cpptraj::Structure::DetermineChirality(double& tors, int* AtomIndices,
                                         int atnum,
                                         Topology const& topIn,
                                         Frame const& frameIn, int debugIn)
{
  tors = 0.0;
  Atom const& atom = topIn[atnum];
  if (atom.Nbonds() < 3) {
    mprinterr("Error: CalcChiralAtomTorsion called for atom %s with less than 3 bonds.\n",
              topIn.AtomMaskName(atnum).c_str());
    return CHIRALITY_ERR;
  }
  // Calculate a priority score for each bonded atom.
  // First just use the atomic number.
  if (debugIn > 0)
    mprintf("DEBUG: Determining priorities around atom %s\n", topIn.AtomMaskName(atnum).c_str());
  std::vector<priority_element> priority;
  for (int idx = 0; idx != atom.Nbonds(); idx++) {
    priority.push_back( priority_element(atom.Bond(idx), topIn[atom.Bond(idx)].AtomicNumber()) );
    if (debugIn > 0)
      mprintf("DEBUG:\t\t%i Priority for %s is %i\n", idx, topIn.AtomMaskName(atom.Bond(idx)).c_str(), priority.back().Priority1());
  }
  // For any identical priorities, need to check who they are bonded to.
  bool depth_limit_hit = false;
  for (int idx1 = 0; idx1 != atom.Nbonds(); idx1++) {
    for (int idx2 = idx1+1; idx2 != atom.Nbonds(); idx2++) {
      if (priority[idx1] == priority[idx2]) {
        bool identical_priorities = true;
        int depth = 1;
        while (identical_priorities) {
          if (debugIn > 0)
            mprintf("DEBUG: Priority of index %i == %i, depth %i\n", idx1, idx2, depth);
          std::vector<bool> Visited(topIn.Natom(), false);
          Visited[atnum] = true;
          priority[idx1].SetPriority2(totalPriority(topIn, atom.Bond(idx1), atom.ResNum(), 0, depth, Visited));
          if (debugIn > 0)
            mprintf("DEBUG:\tPriority2 of %i is %i\n", idx1, priority[idx1].Priority2());

          Visited.assign(topIn.Natom(), false);
          Visited[atnum] = true;
          priority[idx2].SetPriority2(totalPriority(topIn, atom.Bond(idx2), atom.ResNum(), 0, depth, Visited));
          if (debugIn > 0)
            mprintf("DEBUG:\tPriority2 of %i is %i\n", idx2, priority[idx2].Priority2());
          if (priority[idx1] != priority[idx2]) {
            identical_priorities = false;
            break;
          }
          if (depth == 10) {
            if (debugIn > -1)
              mprintf("Warning: Could not determine priority around '%s'\n",
                        topIn.AtomMaskName(atnum).c_str());
            depth_limit_hit = true;
            break;
          }
          depth++;
        } // END while identical priorities
      }
    }
  }
  std::sort(priority.begin(), priority.end());
  if (debugIn > 0) {
    mprintf("DEBUG: Sorted by priority:");
    for (std::vector<priority_element>::const_iterator it = priority.begin();
                                                       it != priority.end(); ++it)
      mprintf(" %s", topIn.AtomMaskName(it->AtNum()).c_str());
    mprintf("\n");
  }

  if (AtomIndices != 0) {
    for (unsigned int ip = 0; ip != priority.size(); ++ip)
      AtomIndices[ip] = priority[ip].AtNum();
  }
  if (depth_limit_hit) return IS_UNKNOWN_CHIRALITY;

  tors = Torsion( frameIn.XYZ(priority[0].AtNum()),
                  frameIn.XYZ(priority[1].AtNum()),
                  frameIn.XYZ(priority[2].AtNum()),
                  frameIn.XYZ(atnum) );
  if (debugIn > 0)
    mprintf("DEBUG: Torsion around '%s' is %f",  topIn.AtomMaskName(atnum).c_str(), tors*Constants::RADDEG);
  ChiralType ret;
  if (tors < 0) {
    ret = IS_S;
    if (debugIn > 0) mprintf(" (S)\n");
  } else {
    ret = IS_R;
    if (debugIn > 0) mprintf(" (R)\n");
  }
  return ret;
}

/** Determine chirality around specified atom. */
Cpptraj::Structure::ChiralType
  Cpptraj::Structure::DetermineChirality(int atnum,
                                         Topology const& topIn,
                                         Frame const& frameIn, int debugIn)
{
  double tors;
  return DetermineChirality(tors, 0, atnum, topIn, frameIn, debugIn);
}

/** Set priority around a specified atom. */
Cpptraj::Structure::ChiralType
  Cpptraj::Structure::SetPriority(std::vector<int>& priority,
                                  int atnum, Topology const& topIn,
                                  Frame const& frameIn, int debugIn)
{
  priority.resize( topIn[atnum].Nbonds() );
  double tors;
  return DetermineChirality(tors, &priority[0], atnum, topIn, frameIn, debugIn);
}

/** Set priority around a specified atom. Do not warn if priority cannot be set. */
Cpptraj::Structure::ChiralType
  Cpptraj::Structure::SetPriority_silent(std::vector<int>& priority,
                                         int atnum, Topology const& topIn,
                                         Frame const& frameIn)
{
  priority.resize( topIn[atnum].Nbonds() );
  double tors;
  return DetermineChirality(tors, &priority[0], atnum, topIn, frameIn, -1);
}

// ===== LEaP Chirality routines ===============================================
/** LEaP routine for determining atom chirality.
  * This is done by crossing A to B and then dotting the
  * result with C. TODO use Chirality in StructureEnum?
  * The chirality of the vectors is determined by the sign of
  * the result, which is determined by whether or not C has
  * a component in the direction AxB or in the opposite direction.
  */
double Cpptraj::Structure::Chirality::VectorAtomChirality(Vec3 const& Center,
                                                          Vec3 const& A, Vec3 const& B, Vec3 const& C)
{
  Vec3 vA = A - Center;
  Vec3 vB = B - Center;
  Vec3 vC = C - Center;
  Vec3 vCross = vA.Cross( vB );
  double dot = vCross * vC;
  if (dot > 0)
    return 1.0;
  else if (dot < 0)
    return -1.0;
  return 0.0;
}

/** LEaP routine for determining atom chirality when positions may or
  * may not be defined. Currently the criteria for chirality is 
  * absolute orientation of the vectors joining this atom to its neighbors.
  * The neighbors are passed as vPA, vPB, vPC, vPD and bA, bB, bC, bD
  * define whether or not the position is defined.
  *
  * This routine calculates the chirality using the defined vectors and
  * then flips the sign depending on which vectors were used to calculate
  * the chirality. If the ATOMs have (fKnown) set then their coordinate
  * is considered to be known.
  */
double Cpptraj::Structure::Chirality::VectorAtomNormalizedChirality(Vec3 const& Center,
                                                                    Vec3 const& vPA, bool bA,
                                                                    Vec3 const& vPB, bool bB,
                                                                    Vec3 const& vPC, bool bC,
                                                                    Vec3 const& vPD, bool bD)
{
  double dChi = 0;
  
  if (!bA) {
    // If A is not known then use B,C,D to calc chirality.
    // The chirality calculated will be negative w.r.t. the
    // correct chirality.
    if (!bB || !bC || !bD) return dChi;
    dChi = -VectorAtomChirality( Center, vPB, vPC, vPD );
    return dChi;
  }

  if (!bB) {
    // If B is not known then use A,C,D to calc chirality.
    // The chirality calculated will be correct.
    if (!bB || !bD) return dChi;
    dChi = VectorAtomChirality( Center, vPA, vPC, vPD );
    return dChi;
  }

  if (!bC) {
    // If C is not known then use A,B,D to calc chirality.
    // The chirality calculated will be negative w.r.t. the
    // correct chirality.
    if (!bD) return dChi;
    dChi = -VectorAtomChirality( Center, vPA, vPB, vPD );
    return dChi;
  }

  dChi = VectorAtomChirality( Center, vPA, vPB, vPC );

  return dChi;
}

/** \return index of atom less than all others but larger than aLast */
static inline int findLeastLargerThan(Atom const& aAtom, int aLast)
{
  int aSmall = -1;
  for (Atom::bond_iterator aCur = aAtom.bondbegin(); aCur != aAtom.bondend(); ++aCur)
  {
    if (aLast != -1) {
      if (aLast >= *aCur) continue;
    }
    if (aSmall == -1)
      aSmall = *aCur;
    else if ( *aCur < aSmall )
      aSmall = *aCur;
  }
  return aSmall;
}

/** Transform the orientation that has been measured with 
 *      respect to the ordering in aaOrig[4] to the ordering
 *      in aaNew[4].  Return the result.
 *
 *      The transformation is done by swapping ATOMs in aaOrig until
 *      the order matches that of aaNew, each time two ATOMs are swapped,
 *      flip the sign of the orientation.
 *
 *      SIDE EFFECT:   The order in aaOrig is changed.
  */
static inline void chiralityTransformOrientation(double dOrig, int* aaOrig, double& dPNew, const int* aaNew)
{
  dPNew = dOrig;
  for (int i=0; i<4; i++ ) {
    int j = i;
    for ( ; j<4; j++ ) {
      if ( aaOrig[j] == aaNew[i] )
        break;
    }
    if ( j >= 4 ) {
      mprinterr("Error: Comparing atoms %i %i %i and %i to atoms %i %i %i and %i.\n",
                aaOrig[0]+1, aaOrig[1]+1, aaOrig[2]+1, aaOrig[3]+1,
                aaNew[0]+1, aaNew[1]+1, aaNew[2]+1, aaNew[3]+1);
      mprinterr("Error: This error may be due to faulty Connection atoms.\n");
      // TODO fatal
    }
    // Swap elements and flip sign
    if ( j != i ) {
      std::swap( aaOrig[j], aaOrig[i] );
      dPNew = -dPNew;
    }
  }
}

/** Sort neighbors of given atom in the same manner as LEaP. */
void Cpptraj::Structure::Chirality::chiralityOrderNeighbors(Atom const& aAtom,
                                                            int& aPAtomA, int& aPAtomB,
                                                            int& aPAtomC, int& aPAtomD)
{
  aPAtomA = -1;
  aPAtomB = -1;
  aPAtomC = -1;
  aPAtomD = -1;

  if (aAtom.Nbonds() < 1) return;

  aPAtomA = findLeastLargerThan(aAtom, -1);
  if (aAtom.Nbonds() < 2) return;

  aPAtomB = findLeastLargerThan(aAtom, aPAtomA);
  if (aAtom.Nbonds() < 3) return;

  aPAtomC = findLeastLargerThan(aAtom, aPAtomB);
  if (aAtom.Nbonds() < 4) return;

  aPAtomD = findLeastLargerThan(aAtom, aPAtomC);
}

/** Transform the chirality which has been measured with
 *  respect to ATOM ID ordering to an arbitrary ordering.
 */
double Cpptraj::Structure::Chirality::chiralityToOrientation(double dChirality, Atom const& aCenter,
                                            int aAtomA, int aAtomB, int aAtomC, int aAtomD)
{
  if (fabs(dChirality) < Constants::SMALL) return 0.0;

  int aaOrig[4];
  chiralityOrderNeighbors( aCenter, aaOrig[0], aaOrig[1], aaOrig[2], aaOrig[3] );

  int aaNew[4];
  aaNew[0] = aAtomA;
  aaNew[1] = aAtomB;
  aaNew[2] = aAtomC;
  aaNew[3] = aAtomD;

  bool newNull = (aaNew[3] == -1);
  bool origNull = (aaOrig[3] == -1);
  if (newNull && !origNull) {
    for (int i = 0; i < 4; i++) {
      bool found = false;
      for (int j=0; j<3; j++) found |= (aaOrig[i] == aaNew[j]);
      if ( !found ) {
        aaNew[3] = aaOrig[i];
        break;
      }
    }
  } else if (!newNull && origNull) {
    mprinterr("Error: Only three neighbors around: aCenter, but orientation has 4\n");
  }

  double dOrient;
  chiralityTransformOrientation( dChirality, aaOrig, dOrient, aaNew );

  return dOrient;
}


