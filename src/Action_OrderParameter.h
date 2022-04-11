// -*- mode: c++; -*-

#ifndef INC_ACTION_ORDERPARAMETER_H
#define INC_ACTION_ORDERPARAMETER_H

#include "Action.h"
#include "OnlineVarT.h"



/** \author Hannes H. Loeffler
  */

class Action_OrderParameter : public Action {
public:
  Action_OrderParameter();

    DispatchObject* Alloc() const {  return (DispatchObject*)new Action_OrderParameter(); }

    void Help() const;

private:
  Action::RetType Init(ArgList&, ActionInit&, int);
  Action::RetType Setup(ActionSetup&);
  Action::RetType DoAction(int, ActionFrame&);
  void Print();

  static const std::string emptystring;
  static const double MAXBOND;
  static const double MAXBOND2;

  CpptrajFile* output_;
  CpptrajFile* taildist_;

  enum DirectionType {DX = 0, DY, DZ};
  DirectionType axis_;

  double delta_;

  AtomMask tailstart_mask_;
  AtomMask tailend_mask_;
  AtomMask unsat_mask_;

  // std::vector<std::pair<AtomMask, bool> > masks_;
  std::vector<AtomMask> masks_;
  std::vector<std::vector<int> > dbonds_;

  bool scd_;
  //ImageOption imageOpt_; ///< Used to decide if imaging should be used.

  unsigned long maxbin_;
  std::vector<Stats<double> > tailhist_;
  std::vector<std::vector<Stats<double> > > orderParameter_;
};
#endif
