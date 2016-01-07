#ifndef INC_ACTION_GRID_H
#define INC_ACTION_GRID_H
#include "Action.h"
#include "DataSet_GridFlt.h"
#include "GridAction.h"
class Action_Grid : public Action, private GridAction {
  public:
    Action_Grid();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Grid(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
#   ifdef MPI
    int ParallelActionInit(Parallel::Comm const& c) { return ParallelGridInit(c, grid_); }
#   endif
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();

    void PrintPDB(double);

    enum NormType { NONE=0, TO_FRAME, TO_DENSITY };
    NormType normalize_;
    double density_;
    double max_;
    double madura_;
    double smooth_;
    unsigned int nframes_;
    int debug_;
    bool invert_;
    AtomMask mask_;
    CpptrajFile* pdbfile_;
    DataSet_GridFlt* grid_;
};
#endif
