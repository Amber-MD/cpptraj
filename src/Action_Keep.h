#ifndef INC_ACTION_KEEP_H
#define INC_ACTION_KEEP_H
#include "Action.h"
#include "ActionTopWriter.h"
class DataSet_string;
/// Keep only specified parts of the system 
class Action_Keep : public Action {
  public:
    Action_Keep();
    ~Action_Keep();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Keep(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    Action::RetType keepBridge(int, ActionFrame&);

    Topology* currentParm_;
    Topology* keepParm_;         ///< Topology for atoms to keep
    Frame keepFrame_;            ///< Frame for atoms to keep

    ActionTopWriter topWriter_;

    //typedef std::pair<int,AtomMask> IdxMaskPairType;
    //typedef std::vector<IdxMaskPairType> IdxMaskPairArray;

    typedef std::vector<char> StatArray;
    StatArray resStat_;          ///< Hold status of each array
    static const char STAT_NONE_;
    static const char STAT_BRIDGERES_;
    static const char STAT_NONBRIDGERES_;

    DataSet_string* bridgeData_; ///< Bridging resdiue ID data set
    int nbridge_;                ///< Number of bridging residues to keep
    std::string bridgeResName_;  ///< Bridging residues name
    int nNonBridgeAtoms_;        ///< Number of non-bridge atoms, for resizing atomsToKeep_
    //IdxMaskPairArray idxMaskPair_; ///< For each res to keep, index into keep mask and mask of atoms to keep

    AtomMask keepMask_;          ///< Mask of atoms to keep no matter what.

    AtomMask atomsToKeep_;       ///< Mask of all atoms to keep.
};
#endif
