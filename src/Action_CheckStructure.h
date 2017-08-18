#ifndef INC_ACTION_CHECKSTRUCTURE_H
#define INC_ACTION_CHECKSTRUCTURE_H
#include "Action.h"
#include "ImagedAction.h"
#include "PairList.h"
class Action_CheckStructure : public Action {
  public:
    Action_CheckStructure();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_CheckStructure(); }
    void Help() const;
    // Interface that can be used outside ActionList
    int SeparateInit(bool, std::string const&, std::string const&, std::string const&,
                     double, double, bool, DataFileList&);
    int SeparateSetup(Topology const&, Box::BoxType, bool);
    int CheckBonds(int, Frame const& currentFrame, Topology const&);
    int CheckOverlap(int, Frame const&, Topology const&);
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();
#   ifdef MPI
    int SyncAction();
#   endif

    enum FmtType { F_ATOM =0, F_BOND };
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

    int PL_CheckOverlap(int, Frame const&, Topology const&);
    void ProcessBondArray(BondArray const&, BondParmArray const&, CharMask const&);
    void SetupBondList(AtomMask const&, Topology const&);
    void WriteProblems(FmtType, int, Topology const&);

#   ifdef _OPENMP
    std::vector<Parray> thread_problemAtoms_;
#   endif
    Parray problemAtoms_;

    PairList pairList_;
    ImagedAction image_; ///< Hold imaging routines and info.
    Parray bondList_;    ///< Array of bonds to check.
    AtomMask Mask1_;     ///< Mask of atoms to check.
    AtomMask Mask2_;     ///< Optional mask of atoms to check against atoms in Mask1
    AtomMask OuterMask_; ///< Mask with the most atoms.
    AtomMask InnerMask_; ///< Mask with fewer atoms.
    double bondoffset_;  ///< Report bonds larger than Req + bondoffset
    double nonbondcut2_; ///< Report distance^2 less than nonbondcut2
    double plcut_;       ///< Pairlist cutoff
    CpptrajFile* outfile_;  ///< Report file.
    Topology* CurrentParm_; ///< Current topology.
    DataSet* num_problems_; ///< Save number of problems each frame
    bool silent_;           ///< If true suppress output
    bool skipBadFrames_;    ///< If true skip frames with problems
    bool bondcheck_;        ///< If true check bonds as well (default)
    bool usePairList_;
    static const char* Fmt_[];
#   ifdef MPI
    Parallel::Comm trajComm_;
    DataSet* ds_fn_; ///< Frame number
    DataSet* ds_pt_; ///< Problem type
    DataSet* ds_a1_; ///< Atom 1
    DataSet* ds_n1_; ///< Name 1
    DataSet* ds_a2_; ///< Atom 2
    DataSet* ds_n2_; ///< Name 2
    DataSet* ds_d_;  ///< Distance
    int idx_;        ///< Index into ds_X data sets.
#   endif
};
#endif
