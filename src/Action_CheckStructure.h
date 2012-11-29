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

    /// Used to cache bond parameters
    struct bond_list {
      double req;
      int atom1;
      int atom2;
      int param;
    };
    std::vector<bond_list> bondL_;
    // Sort first by atom1, then by atom2
    struct bond_list_cmp {
      inline bool operator()(bond_list first, bond_list second) const {
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
    CpptrajFile outfile_;
    Topology* CurrentParm_;
    int debug_;

    void SetupBondlist(std::vector<int> const&);
};
#endif  
