#ifndef INC_STRUCTURECHECK_H
#define INC_STRUCTURECHECK_H
#include "PairList.h"
#include "ImagedAction.h"
/// Used to count potential structure problems.
class StructureCheck {
  public:
    StructureCheck();
    /// Options: imageOn, checkBonds, saveProblems, mask1, mask2, ovrlpCut, bndLenOffset, PListCut
    int SetOptions(bool, bool, bool, std::string const&, std::string const&,
                   double, double, double);
    /// Setup for given topology and box.
    int Setup(Topology const&, Box const&);
    /// \return Number of abnormal bonds.
    int CheckBonds(Frame const&);
    /// \return Number of atomic overlaps.
    int CheckOverlaps(Frame const&);

    AtomMask const& Mask1()     const { return Mask1_; }
    AtomMask const& Mask2()     const { return Mask2_; }
    ImagedAction const& Image() const { return image_; }
    bool CheckBonds()           const { return bondcheck_; }
    double BondOffset()         const { return bondoffset_; }
    double NonBondCut2()        const { return nonbondcut2_; }
    double PairListCut()        const { return plcut_; }
    unsigned int Nbonds()       const { return bondList_.size(); }
#   ifdef _OPENMP
    unsigned int Nthreads()     const { return thread_problemAtoms_.size(); }
#   endif
    // -------------------------------------------
    /// Store problems between atoms. Also used to cache bond parameters.
    class Problem {
      public:
        Problem() : dist_(0.0), atom1_(-1), atom2_(-1) {}
        Problem(int a1, int a2, double d) : dist_(d) {
          if (a1 < a2) {
            atom1_ = a1; atom2_ = a2;
          } else {
            atom1_ = a2; atom2_ = a1;
          }
        }
        bool operator<(Problem const& rhs) const {
          if (atom1_ == rhs.atom1_)
            return (atom2_ < rhs.atom2_);
          else
            return (atom1_ < rhs.atom1_);
        }
        double D()  const { return dist_;  }
        int A1()    const { return atom1_; }
        int A2()    const { return atom2_; }
#       ifdef MPI
        const double* Dptr() const { return &dist_; }
#       endif
      private:
        double dist_; ///< Distance / Bond cutoff (Req+bondoffset)^2
        int atom1_;   ///< First atom
        int atom2_;   ///< Second atom
    };
    typedef std::vector<Problem> Parray;
    // -------------------------------------------

    typedef Parray::const_iterator const_iterator;
    /// Iterator to beginning of list of problems from latest call to CheckX
    const_iterator begin() const { return problemAtoms_.begin(); }
    /// Iterator to end of list of problems from latest call to CheckX
    const_iterator end()   const { return problemAtoms_.end();   }
  private:
    /// Type of overlap check
    enum CheckType { NO_PL_1_MASK=0, NO_PL_2_MASKS, PL_1_MASK };
    /// Add selected bonds in BondArray to list to be checked.
    void ProcessBondArray(BondArray const&, BondParmArray const&, CharMask const&);
    /// Add selected bonds in topology to list to be checked.
    void SetupBondList(AtomMask const&, Topology const&);
    /// PairList version of CheckOverlap, 1 mask
    int PL1_CheckOverlap(Frame const&, Matrix_3x3 const&, Matrix_3x3 const&);
    /// Non-pairlist version of CheckOverlap, 1 mask
    int Mask1_CheckOverlap(Frame const&, Matrix_3x3 const&, Matrix_3x3 const&);
    /// Non-pairlist version of CheckOverlap, 2 masks
    int Mask2_CheckOverlap(Frame const&, Matrix_3x3 const&, Matrix_3x3 const&);
    /// Sort problem list; combine results from threads for OpenMP
    void ConsolidateProblems();
#   ifdef _OPENMP
    std::vector<Parray> thread_problemAtoms_;
#   endif
    Parray problemAtoms_;

    PairList pairList_;     ///< Atom pair list
    ImagedAction image_;    ///< Hold imaging routines and info.
    Parray bondList_;       ///< Array of bonds to check.
    AtomMask Mask1_;        ///< Mask of atoms to check.
    AtomMask Mask2_;        ///< Optional mask of atoms to check against atoms in Mask1
    AtomMask OuterMask_;    ///< Mask with the most atoms.
    AtomMask InnerMask_;    ///< Mask with fewer atoms.
    double bondoffset_;     ///< Report bonds larger than Req + bondoffset
    double nonbondcut2_;    ///< Report distance^2 less than nonbondcut2
    double plcut_;          ///< Pairlist cutoff
    CheckType checkType_;   ///< Type of atom overlap check
    bool bondcheck_;        ///< If true check bonds as well
    bool saveProblems_;     ///< If true save problems in problemAtoms_
};
#endif
