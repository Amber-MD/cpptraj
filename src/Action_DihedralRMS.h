#ifndef INC_ACTION_DIHEDRALRMS_H
#define INC_ACTION_DIHEDRALRMS_H
#include "Action.h"
#include "DihedralSearch.h"
#include "ReferenceAction.h"
/// Calculate dihedral RMSD to reference 
class Action_DihedralRMS : public Action {
  public:
    Action_DihedralRMS();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_DihedralRMS(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    typedef std::vector<double> Darray;

    Range GetActualRange(Topology const&, Range const&) const;
    int SetupRefDihedrals(Topology const&);
    int CalcRefDihedrals(Frame const&);

    DihedralSearch dihSearch_; ///< Used to select specified dihedrals in target
    DihedralSearch refSearch_; ///< Used to select specified dihedrals in reference
    ReferenceAction REF_;      ///< Hold reference structure
    Darray refVals_;           ///< Reference dihedral values in radians
    Range tgtRange_;
    Range refRange_;
    DataSet* dataOut_;         ///< Output data set
    int debug_;
};
#endif
