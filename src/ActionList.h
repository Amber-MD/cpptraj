#ifndef INC_ACTIONLIST_H
#define INC_ACTIONLIST_H
#include "Action.h"
// Class: ActionList
/// Hold actions that will be performed every frame.
/** This class is responsible for holding all actions that will be performed
  * during the course of trajectory processing. When an argument list 
  * containing a recognized action keyword is passed to AddAction, the
  * corresponding action class is added to the list along with the argument
  * list. Actions are not initialized until all parm and trajectory 
  * information has been read in. Actions are then set up for parm files,
  * and re-set up each time the parm changes. ActionList needs to know about 
  * any new actions that are implemented.
  */
class ActionList {
  public:

    ActionList();
    ~ActionList();

    /// Set the debug level for actions
    void SetDebug(int);
    /// Add an action to the action list.
    int AddAction(ArgList &);
    /// Initialize actions
    int Init(DataSetList *, FrameList *, DataFileList *, TopologyList*,bool);
    /// Set up actions for the given parm
    int Setup(Topology **);
    /// Perform actions on the given frame
    bool DoActions(Frame **, int);
    /// Call print for each action
    void Print();
  private:
    /// List of actions
    std::vector<Action*> actionlist_;
    typedef std::vector<Action*>::iterator action_it;
    /// debug level for actions
    int debug_;
};
#endif
