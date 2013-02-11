#ifndef INC_ACTION_MAKESTRUCTURE_H
#define INC_ACTION_MAKESTRUCTURE_H
#include "Action.h"
#include "DihedralSearch.h"
class Action_MakeStructure : public Action {
  public:
    Action_MakeStructure();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_MakeStructure(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    enum ssType { ALPHA = 0, LEFT, PP2, EXT, 
                  TI, TII, TVIII, TIp, TIIp, TVIa1, TVIa2, TVIb, NSS };
    struct SecStructHolder {
      Range resRange;                ///< Residues to set phi/psi for.
      DihedralSearch dihSearch_;     ///< Used to search for specified dihedrals
      std::vector<AtomMask> Rmasks_; ///< Masks of atoms to hold fixed during rotation.
      std::vector<float> thetas_;    ///< theta for each dihedral.
      ssType type;                   ///< Type of SS.
      std::string name;              ///< SS type name.
    };
    std::vector<SecStructHolder> secstruct_;
};
#endif
