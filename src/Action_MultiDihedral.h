#ifndef INC_MULTIDIHEDRAL_H
#define INC_MULTIDIHEDRAL_H
#include "Action.h"
#include "Range.h"
class Action_MultiDihedral : public Action {
  public:
    Action_MultiDihedral();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_MultiDihedral(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    void FindDihedralAtoms(Topology*, int, int, NameType const&, NameType const&,
                           NameType const&, NameType const&, std::string const&);

    std::vector<int> maskAtoms_;  ///< Selected atoms, 4 per dihedral
    std::vector<DataSet*> data_;  ///< Output DataSets, 1 per dihedral
    bool useMass_;
    bool range360_;
    bool findPHI_;
    bool findPSI_;
    Range resRange_;
    std::string dsetname_;
    DataFile* outfile_;
    // TODO: Replace these with new DataSet type
    DataSetList* masterDSL_;
};
#endif 
