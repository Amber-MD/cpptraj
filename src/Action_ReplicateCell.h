#ifndef INC_ACTION_REPLICATECELL_H
#define INC_ACTION_REPLICATECELL_H
#include "Action.h"
#include "Trajout_Single.h"
#include "ActionTopWriter.h"
// Forward declares
class DataSet_Coords;
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

    typedef std::vector<int> Iarray;
    Iarray directionArray_;     ///< Array of directions to replicate (x,y,z)
    ActionTopWriter topWriter_; ///< Used to write replicated cell topology
    Trajout_Single outtraj_;    ///< Output combined cell traj
    DataSet_Coords* coords_;    ///< Combined cell COORDS DataSet
    AtomMask Mask1_;            ///< Mask of atoms to replicate
    int ncopies_;               ///< Total # of replications to make
    bool writeTraj_;            ///< If true, write output combined cell traj
    Topology combinedTop_;      ///< Combined cell topology
    Frame combinedFrame_;       ///< Combined cell frame
};
#endif
