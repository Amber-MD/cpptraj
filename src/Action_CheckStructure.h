#ifndef INC_ACTION_CHECKSTRUCTURE_H
#define INC_ACTION_CHECKSTRUCTURE_H
#include "Action.h"
#include "ImagedAction.h"
// Class: Action_CheckStructure 
/// Action to check bond lengths and bad overlaps between non-bonded atoms 
class Action_CheckStructure: public Action, ImagedAction {
  public:
    Action_CheckStructure();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_CheckStructure(); }
    static void Help();
    ~Action_CheckStructure();
    // These are made public for use in other actions (e.g. Action_DihedralScan)
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    int CheckFrame(int, Frame const&);
    void Print() {}
  private:
    Action::RetType DoAction(int, Frame*, Frame**);

    void SetupBondlist(BondArray const&, BondParmArray const&, AtomMask const&);
    /// Used to cache bond parameters
    struct bond_list {
      double req; ///< Hold (req + bondoffset)^2
#     ifdef _OPENMP
      double D2;  ///< Hold distance^2 for this pair.
      int problem;///< Hold detected problem type for this pair.
#     endif
      int atom1;
      int atom2;
    };
    typedef std::vector<bond_list> BondListType;
    BondListType bondL_;
    /// Sort first by atom1, then by atom2
    struct bond_list_cmp {
      inline bool operator()(bond_list const& first, bond_list const& second) const {
        if (first.atom1 < second.atom1) {
          return true;
        } else if (first.atom1 == second.atom1) {
          if (first.atom2 < second.atom2) return true;
        } 
        return false;
      }
    };

    AtomMask Mask1_;
    double bondoffset_;
    double nonbondcut2_;
    bool bondcheck_;
    bool silent_;
    bool skipBadFrames_;
    CpptrajFile outfile_;
    Topology* CurrentParm_;
    int debug_;
};
#endif  
