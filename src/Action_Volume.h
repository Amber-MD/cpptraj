#ifndef INC_ACTION_VOLUME_H
#define INC_ACTION_VOLUME_H
#include "Action.h"
#include "ImagedAction.h"
/// Calculate unit cell volume. 
class Action_Volume: public Action {
  public:
    Action_Volume();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Volume(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    ImagedAction image_;
    DataSet *vol_;
    double sum_;
    double sum2_;
    int nframes_; 
};
#endif
