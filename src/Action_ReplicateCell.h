#ifndef INC_ACTION_REPLICATECELL_H
#define INC_ACTION_REPLICATECELL_H
#include "Action.h"
#include "ImagedAction.h"
#include "Trajout_Single.h"
#include "DataSet_Coords.h"
/// Action to replicate unit cell in specified directions. 
class Action_ReplicateCell: public Action {
  public:
    Action_ReplicateCell();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_ReplicateCell(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    ImagedAction image_;
    Matrix_3x3 ucell_, recip_;
    typedef std::vector<int> Iarray;
    Iarray directionArray_;
    std::string trajfilename_;
    std::string parmfilename_;
    Trajout_Single outtraj_;
    DataSet_Coords* coords_;
    AtomMask Mask1_;
    int ncopies_;
    Topology combinedTop_;
    Frame combinedFrame_;
    ArgList trajArgs_;
    int ensembleNum_;
};
#endif
