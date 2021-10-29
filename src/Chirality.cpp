#include "Chirality.h"
#include "Constants.h"
#include "CpptrajStdio.h"
#include "Frame.h"
#include "Topology.h"
#include "TorsionRoutines.h"
#include <vector>
#include <algorithm> // sort

using namespace Cpptraj;

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
Chirality::ChiralType Chirality::DetermineChirality(double& tors, int atnum,
                                                    Topology const& topIn,
                                                    Frame const& frameIn)
{
  tors = 0.0;
  Atom const& atom = topIn[atnum];
  if (atom.Nbonds() < 3) {
    mprinterr("Error: CalcChiralAtomTorsion called for atom %s with less than 3 bonds.\n",
              topIn.AtomMaskName(atnum).c_str());
    return ERR;
  }
  // Calculate a priority score for each bonded atom.
  // First just use the atomic number.
  mprintf("DEBUG: Determining priorities around atom %s\n", topIn.AtomMaskName(atnum).c_str());
  std::vector<priority_element> priority;
  for (int idx = 0; idx != atom.Nbonds(); idx++) {
    priority.push_back( priority_element(atom.Bond(idx), topIn[atom.Bond(idx)].AtomicNumber()) );
    mprintf("DEBUG:\t\t%i Priority for %s is %i\n", idx, topIn.AtomMaskName(atom.Bond(idx)).c_str(), priority.back().Priority1());
  }
  // For any identical priorities, need to check who they are bonded to.
  for (int idx1 = 0; idx1 != atom.Nbonds(); idx1++) {
    for (int idx2 = idx1+1; idx2 != atom.Nbonds(); idx2++) {
      if (priority[idx1] == priority[idx2]) {
        bool identical_priorities = true;
        int depth = 1;
        while (identical_priorities) {
          mprintf("DEBUG: Priority of index %i == %i, depth %i\n", idx1, idx2, depth);
          std::vector<bool> Visited(topIn.Natom(), false);
          Visited[atnum] = true;
          priority[idx1].SetPriority2(totalPriority(topIn, atom.Bond(idx1), atom.ResNum(), 0, depth, Visited));
          mprintf("DEBUG:\tPriority2 of %i is %i\n", idx1, priority[idx1].Priority2());

          Visited.assign(topIn.Natom(), false);
          Visited[atnum] = true;
          priority[idx2].SetPriority2(totalPriority(topIn, atom.Bond(idx2), atom.ResNum(), 0, depth, Visited));
          mprintf("DEBUG:\tPriority2 of %i is %i\n", idx2, priority[idx2].Priority2());
          if (priority[idx1] != priority[idx2]) {
            identical_priorities = false;
            break;
          }
          if (depth == 10) {
            mprinterr("Error: Could not determine priority around '%s'\n",
                      topIn.AtomMaskName(atnum).c_str());
            return ERR;
          }
          depth++;
        } // END while identical priorities
      }
    }
  }
  std::sort(priority.begin(), priority.end());
  mprintf("DEBUG: Sorted by priority:");
  for (std::vector<priority_element>::const_iterator it = priority.begin();
                                                     it != priority.end(); ++it)
    mprintf(" %s", topIn.AtomMaskName(it->AtNum()).c_str());
  mprintf("\n");

  tors = Torsion( frameIn.XYZ(priority[0].AtNum()),
                  frameIn.XYZ(priority[1].AtNum()),
                  frameIn.XYZ(priority[2].AtNum()),
                  frameIn.XYZ(atnum) );
  mprintf("DEBUG: Torsion around '%s' is %f",  topIn.AtomMaskName(atnum).c_str(), tors*Constants::RADDEG);
  ChiralType ret;
  if (tors < 0) {
    ret = IS_S;
    mprintf(" (S)\n");
  } else {
    ret = IS_R;
    mprintf(" (R)\n");
  }
  return ret;
}

/** Determine chirality around specified atom. */
Chirality::ChiralType Chirality::DetermineChirality(int atnum,
                                                    Topology const& topIn,
                                                    Frame const& frameIn)
{
  double tors;
  return DetermineChirality(tors, atnum, topIn, frameIn);
}
