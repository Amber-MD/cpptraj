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
    class SS_TYPE {
      public: 
        SS_TYPE() {}
        SS_TYPE(double ph,double ps,double ph2,double ps2,int t,const char* n) :
          phi(ph), psi(ps), phi2(ph2), psi2(ps2), isTurn(t), name(n) {}
        double phi, psi, phi2, psi2;
        int isTurn;
        const char* name;
    };
    std::vector<SS_TYPE> SS;

    struct SecStructHolder {
      Range resRange;                ///< Residues to set phi/psi for.
      DihedralSearch dihSearch_;     ///< Used to search for specified dihedrals
      std::vector<AtomMask> Rmasks_; ///< Masks of atoms to hold fixed during rotation.
      std::vector<float> thetas_;    ///< theta for each dihedral.
      std::vector<SS_TYPE>::const_iterator type; ///< Pointer to SS_TYPE.
      std::string name;              ///< SS type name.
    };
    std::vector<SecStructHolder> secstruct_;
};
#endif
