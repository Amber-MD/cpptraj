#ifndef INC_ACTION_CHECKSTRUCTURE_H
#define INC_ACTION_CHECKSTRUCTURE_H
#include "Action.h"
#include "ImagedAction.h"
class Action_CheckStructure : public Action {
  public:
    Action_CheckStructure();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_CheckStructure(); }
    static void Help();
    // Interface that can be used outside ActionList
    int SeparateInit(bool, std::string const&, std::string const&, std::string const&,
                     double, double, bool, DataFileList&);
    int SeparateSetup(Topology const&, bool);
    int CheckBonds(int, Frame const& currentFrame, Topology const&);
    int CheckOverlap(int, Frame const&, Topology const&);
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    void ProcessBondArray(BondArray const&, BondParmArray const&, AtomMask const&);
    void SetupBondList(AtomMask const&, Topology const&);
    /// Used to cache bond parameters
    struct BondType {
      double Req_off2_; ///< Bond cutoff (Req+bondoffset)^2
      int a1_;          ///< First atom in bond
      int a2_;          ///< Second atom in bond
    };
    typedef std::vector<BondType> BondList;

    ImagedAction image_; ///< Hold imaging routines and info.
    BondList bondList_;  ///< Array of bonds to check.
    AtomMask Mask1_;     ///< Mask of atoms to check.
    AtomMask Mask2_;     ///< Optional mask of atoms to check against atoms in Mask1
    AtomMask OuterMask_; ///< Mask with the most atoms.
    AtomMask InnerMask_; ///< Mask with fewer atoms.
    double bondoffset_;  ///< Report bonds larger than Req + bondoffset
    double nonbondcut2_; ///< Report distance^2 less than nonbondcut2
    CpptrajFile* outfile_;  ///< Report file.
    Topology* CurrentParm_; ///< Current topology.
    bool silent_;           ///< If true suppress output
    bool skipBadFrames_;    ///< If true skip frames with problems
    bool bondcheck_;        ///< If true check bonds as well (default)
};
#endif
