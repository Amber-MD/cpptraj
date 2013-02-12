#ifndef INC_ACTION_GIST_H
#define INC_ACTION_GIST_H
#include "Action.h"
// Class: Action_GIST
/// Calculate water energy and entropy
//class Action_GIST: public Action, ImagedAction  {
class Action_GIST: public Action  {
  public:
    Action_GIST();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_GIST(); }
    static void Help();

  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    DataSet* gist_;  ///< Will hold DataSet of calculated gist values
    bool watermodel_;
    bool useTIP3P_;
    bool useTIP4P_;
    bool useTIP4PEW_;

    //    Box gridbox_;    
    Vec3 gridcntr_;    
    Vec3 griddim_; 
    double gridspacn_;
    
};
#endif
