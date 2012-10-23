#ifndef INC_ACTION_WATERSHELL_H
#define INC_ACTION_WATERSHELL_H
#include "Action.h"
#include "ImagedAction.h"
/// Calculate number of solvent residues in 1st/2nd solvation shell.
class Action_Watershell : public Action, ImagedAction {
  public:
    Action_Watershell();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Watershell(); }
    static void Help();

  private:
    AtomMask soluteMask_;
    AtomMask solventMask_;
    std::string solventmaskexpr_;
    std::vector<int> activeResidues_;
    std::vector<int> lower_;
    std::vector<int> upper_;
    double lowerCutoff_;
    double upperCutoff_;
    std::string filename_;
    Topology* CurrentParm_;

    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();
};
#endif
