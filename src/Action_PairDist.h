// -*- mode: c++; -*-

#ifndef INC_ACTION_PAIRDIST_H
#define INC_ACTION_PAIRDIST_H

#include "Action.h"
#include "ImagedAction.h"
#include "OnlineVarT.h"



/** \author Hannes H. Loeffler
  */

class Action_PairDist : public Action, ImagedAction {
 public:
  Action_PairDist();

    DispatchObject* Alloc() const { return (DispatchObject*)new Action_PairDist(); }

    void Help() const;

 private:
  Action::RetType Init(ArgList&, ActionInit&, int);
  Action::RetType Setup(ActionSetup&);
  Action::RetType DoAction(int, ActionFrame&);
# ifdef MPI
  int SyncAction();
  Parallel::Comm trajComm_;
# endif
  void Print();

  CpptrajFile* output_;

  DataSet *Pr_;    /// distance vs P(r)
  DataSet *std_;   /// distance vs std

  AtomMask mask1_;
  AtomMask mask2_;

  double delta_;		// resolution

  std::vector<Stats<double> > histogram_;
  unsigned long maxbin_;

  bool same_mask_;
  unsigned long ub1_;
  unsigned long ub2_;
};
#endif
