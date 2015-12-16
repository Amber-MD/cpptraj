#ifndef INC_MULTIDIHEDRAL_H
#define INC_MULTIDIHEDRAL_H
#include "Action.h"
#include "DihedralSearch.h"
class Action_MultiDihedral : public Action {
  public:
    Action_MultiDihedral();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_MultiDihedral(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    double minTorsion_;           ///< Values less than this will be shifted +360
    int debug_;
    DihedralSearch dihSearch_;    ///< Used to search for specified dihedrals
    std::vector<DataSet*> data_;  ///< Output DataSets, 1 per dihedral.
    Range resRange_;              ///< Residues to search for dihedrals.
    std::string dsetname_;
    DataFile* outfile_;
    // TODO: Replace these with new DataSet type
    DataSetList* masterDSL_;
};
#endif 
