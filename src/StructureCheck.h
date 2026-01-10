#ifndef INC_STRUCTURECHECK_H
#define INC_STRUCTURECHECK_H
#include "PairList.h"
#include "ImageOption.h"
#include "ExclusionArray.h"
#include "AtomMask.h"
#include "ParameterTypes.h"
#include "Structure/RingFinder.h"
#include "Timer.h"
// Forward declares
class Topology;
class CharMask;
class CpptrajFile;
namespace Cpptraj {
namespace Structure {
class LeastSquaresPlane;
}
}
/// Used to count potential structure problems.
class StructureCheck {
  public:
    // -------------------------------------------
    /// Cache heavy atom bond indices for checking with rings
    class Btype {
      public:
        Btype(int a1, int a2) : a1_(a1), a2_(a2) {}
        int A1() const { return a1_; }
        int A2() const { return a2_; }
      private:
        int a1_;
        int a2_;
    };
    // -------------------------------------------
    /// CONSTRUCTOR
    StructureCheck();
    /// Options: imageOn, checkBonds, saveProblems, check extrapoints, debug, mask1, mask2, extra points exclude mask, ovrlpCut, bndLenOffset, minBndLenOffset, PListCut
    int SetOptions(bool, bool, bool, bool, int, std::string const&, std::string const&, std::string const&,
                   double, double, double, double, bool, double, double, double);
    /// Options (no extra points): imageOn, checkBonds, saveProblems, debug, mask1, mask2, ovrlpCut, bndLenOffset, minBndLenOffset, PListCut
    int SetOptions(bool, bool, bool, int, std::string const&, std::string const&,
                   double, double, double, double, bool, double, double, double);

    /// Setup for given topology and box.
    int Setup(Topology const&, Box const&);
    /// \return Number of abnormal bonds.
    int CheckBonds(Frame const&);
    /// \return Number of atomic overlaps.
    int CheckOverlaps(Frame const&);
    /// Check if any bonds are passing through rings.
    int CheckRings(Frame const&);
    /// Check if any of the given bonds are passing through given rings.
    int CheckRings(Frame const&, Cpptraj::Structure::RingFinder const&, std::vector<Btype> const&);
    /// Write existing problems to the given file
    void WriteProblemsToFile(CpptrajFile*, int, Topology const&) const;


    void PrintTiming(int, double) const;
    /// Format of problems currently stored in problemAtoms_
    enum FmtType { F_ATOM =0, F_BOND, F_RING };

    AtomMask const& Mask1()     const { return Mask1_; }
    AtomMask const& Mask2()     const { return Mask2_; }
    ImageOption const& ImageOpt() const { return imageOpt_; }
    bool CheckBonds()           const { return bondcheck_; }
    double BondOffset()         const { return bondoffset_; }
    double BondMinOffset()      const { return bondMinOffset_; }
    bool CheckRings()           const { return ringcheck_; }
    double RingShortDist() const;
    double RingDist() const;
    double RingAngleCut_Deg() const;
    double NonBondCut2()        const { return nonbondcut2_; }
    double PairListCut()        const { return plcut_; }
    unsigned int Nbonds()       const { return bondList_.size(); }
    int Debug()                 const { return debug_; }
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
        /// Set problem atoms and distance, do not sort atoms
        void SetProb(int a1, int a2, double d) { atom1_ = a1; atom2_ = a2; dist_ = d; }
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
    /// Format strings corresponding to FmtType
    static const char* Fmt_[];
    /// Type of overlap check
    enum CheckType { NO_PL_1_MASK=0, NO_PL_2_MASKS, PL_1_MASK };
    /// \return mask string based on if extra points are being excluded
    std::string checkMaskStr(std::string const&) const;
    /// Add selected bonds in BondArray to list to be checked.
    void ProcessBondArray(BondArray const&, BondParmArray const&, CharMask const&, std::vector<Atom> const&);
    /// Add selected bonds in topology to list to be checked.
    void SetupBondList(AtomMask const&, Topology const&);
    /// Check for intersection between bond and ring
    void ring_bond_check(int&, double, Btype const&, Vec3 const&, AtomMask const&,
                         Cpptraj::Structure::LeastSquaresPlane const&
#                        ifdef _OPENMP
                         , int
#                        endif
                        );
    /// \return True if given bond not in given ring mask
    static inline bool check_bond_not_in_ring(AtomMask const&, Btype const&);
    /// Check for intersections between bonds and rings, no pair list
    int checkRings_NoPL(Frame const&, Cpptraj::Structure::RingFinder const&,
                        std::vector<Btype> const&,
                        std::vector<Cpptraj::Structure::LeastSquaresPlane> const&);
    /// Check for intersections between bonds and rings, use pair list
    int checkRings_PL(Frame const&, Cpptraj::Structure::RingFinder const&,
                      std::vector<Btype> const&,
                      std::vector<Cpptraj::Structure::LeastSquaresPlane> const&);
    /// PairList version of CheckOverlap, 1 mask
    int PL1_CheckOverlap(Frame const&);
    /// Non-pairlist version of CheckOverlap, 1 mask
    int Mask1_CheckOverlap(Frame const&);
    /// Non-pairlist version of CheckOverlap, 2 masks
    int Mask2_CheckOverlap(Frame const&);
    /// Sort problem list; combine results from threads for OpenMP
    void ConsolidateProblems();
    /// Check for/record non-bonded interaction problem
    inline void DistanceCheck(Frame const&, int, int, Parray&, int&) const;

    std::vector<Btype> ringBonds_;

#   ifdef _OPENMP
    std::vector<Parray> thread_problemAtoms_;
#   endif
    Parray problemAtoms_;

    PairList pairList_;     ///< Atom pair list
    ExclusionArray Excluded_; ///< Hold excluded atoms for pair list.
    ImageOption imageOpt_;    ///< Used to determine if imaging should be used.
    Cpptraj::Structure::RingFinder rings_; ///< Used to find rings
    Parray bondList_;       ///< Array of bonds to check.
    AtomMask Mask1_;        ///< Mask of atoms to check.
    AtomMask Mask2_;        ///< Optional mask of atoms to check against atoms in Mask1
    AtomMask OuterMask_;    ///< Mask with the most atoms.
    AtomMask InnerMask_;    ///< Mask with fewer atoms.
    std::string XP_Exclude_Mask_;   ///< Expression to deselect extra points
    double bondoffset_;     ///< Report bonds larger than Req + bondoffset_
    double bondMinOffset_;  ///< Report bonds less than Req - bondMinOffset_
    double nonbondcut2_;    ///< Report distance^2 less than nonbondcut2_
    double plcut_;          ///< Pairlist cutoff
    double ring_shortd2_;   ///< Ring center to bond center short distance (Ang) cutoff squared.
    double ring_dcut2_;     ///< Ring center to bond center distance (Ang) cutoff squared.
    double ring_acut_;      ///< Ring perpendicular vector to bond vector angle cutoff (rad).
    CheckType checkType_;   ///< Type of atom overlap check
    int debug_;             ///< Debug level.
    bool bondcheck_;        ///< If true check bonds as well
    bool ringcheck_;        ///< If true check bond/ring intersections
    bool saveProblems_;     ///< If true save problems in problemAtoms_
    bool checkExtraPts_;    ///< If true check extra points.
    FmtType lastFmt_;       ///< Format of problems currently stored in problemAtoms_

    Timer t_setup_;
    Timer t_setup_bonds_;
    Timer t_setup_exclusion_;
    Timer t_setup_ringfinder_;
};
#endif
