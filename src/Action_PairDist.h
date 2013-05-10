#ifndef INC_ACTION_PAIRDIST_H
#define INC_ACTION_PAIRDIST_H

#include "Action.h"
#include "ImagedAction.h"


/** \author Hannes H. Loeffler
  */

class Action_PairDist : public Action, ImagedAction {
 public:
  Action_PairDist();

  static DispatchObject* Alloc() {
    return (DispatchObject*)new Action_PairDist();
  }
    
  static void Help();

 private:
  Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
		       DataFileList*, int);
  Action::RetType Setup(Topology*, Topology**);
  Action::RetType DoAction(int, Frame*, Frame**);
  void Print();

  DataSet* Pr_;

  AtomMask mask1_;
  AtomMask mask2_;

  double delta_;		// resolution

  std::vector<double> histogram_;
  long int maxbin_;
};
#endif    
