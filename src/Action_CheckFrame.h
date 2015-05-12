#ifndef INC_ACTION_CHECKFRAME_H
#define INC_ACTION_CHECKFRAME_H
#include "Action.h"
#include "ImagedAction.h"
class Action_CheckFrame : public Action {
  public:
    Action_CheckFrame();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_CheckFrame(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    void ProcessBondArray(BondArray const&, BondParmArray const&, AtomMask const&);
    void SetupBondList(AtomMask const&, Topology const&);
    int CheckBonds(int, Frame const& currentFrame, Topology const&);
    int CheckOverlap(int, Frame const&, Topology const&);
    /// Used to cache bond parameters
    struct BondType {
      double Req_off2_; ///< Bond (Req+bondoffset)^2
      int a1_;          ///< First atom in bond
      int a2_;          ///< Second atom in bond
    };
    typedef std::vector<BondType> BondList;

//#   ifdef _OPENMP
//    typedef std::vector<BondList> ProblemList;
//    ProblemList Problems_;
//#   endif
    ImagedAction image_;
    BondList bondList_; ///< Array of bonds to check.
    AtomMask Mask1_; ///< Mask of atoms to check.
    AtomMask Mask2_; ///< Optional mask of atoms to check against atoms in Mask1
    AtomMask OuterMask_; ///< Mask with the most atoms.
    AtomMask InnerMask_; ///< Mask with fewer atoms.
    double bondoffset_; ///< Report bonds larger than Req + bondoffset
    double nonbondcut2_; ///< Report distance^2 less than nonbondcut2
    CpptrajFile* outfile_; ///< Report file.
    Topology* CurrentParm_; ///< Current topology.
    int debug_; ///< Debug level
    bool silent_; ///< If true suppress output
    bool skipBadFrames_; ///< If true skip frames with problems
    bool bondcheck_; ///< If true check bonds as well (default)
};
#endif
