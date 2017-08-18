#ifndef INC_ACTION_CHECKSTRUCTURE_H
#define INC_ACTION_CHECKSTRUCTURE_H
#include "Action.h"
#include "StructureCheck.h"
class Action_CheckStructure : public Action {
  public:
    Action_CheckStructure();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_CheckStructure(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();
#   ifdef MPI
    int SyncAction();
#   endif

    enum FmtType { F_ATOM =0, F_BOND };

    void WriteProblems(FmtType, int, Topology const&);

    StructureCheck check_;  ///< Structure checker
    CpptrajFile* outfile_;  ///< Report file.
    Topology* CurrentParm_; ///< Current topology.
    DataSet* num_problems_; ///< Save number of problems each frame
    bool silent_;           ///< If true suppress output
    bool skipBadFrames_;    ///< If true skip frames with problems
    static const char* Fmt_[]; ///< Output format strings
#   ifdef MPI
    Parallel::Comm trajComm_;
    DataSet* ds_fn_; ///< Frame number
    DataSet* ds_pt_; ///< Problem type
    DataSet* ds_a1_; ///< Atom 1
    DataSet* ds_n1_; ///< Name 1
    DataSet* ds_a2_; ///< Atom 2
    DataSet* ds_n2_; ///< Name 2
    DataSet* ds_d_;  ///< Distance
    int idx_;        ///< Index into ds_X data sets.
#   endif
};
#endif
