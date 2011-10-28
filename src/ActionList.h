#ifndef INC_ACTIONLIST_H
#define INC_ACTIONLIST_H
/// Class: ActionList
/// Hold actions that will be performed every frame. Each time a new parm is
/// loaded the setup routin will be called.
#include "Action.h"
class ActionList {
    Action **actionlist;
    int Naction;
    int debug;

  public:

    ActionList();
    ~ActionList();

    void SetDebug(int);
    int AddAction(ArgList &);
    int Init(DataSetList *, FrameList *, DataFileList *, ParmFileList*);
    int Setup(AmberParm **);
    void DoActions(Frame **, int);
    void Print();
};
#endif
