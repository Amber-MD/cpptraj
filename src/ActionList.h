#ifndef INC_ACTIONLIST_H
#define INC_ACTIONLIST_H
#include "Action.h"
// Class: ActionList
/// Hold actions that will be performed every frame.
/** This class is responsible for holding all actions that will be performed
  * during the course of trajectory processing.
  */
class ActionList {
  public:
    ActionList();
    ~ActionList();
    /// Clear the list
    void Clear();
    /// Set the debug level for actions.
    void SetDebug(int);
    int Debug() const { return debug_; }
    /// Add given action to the action list
    int AddAction( Action* actIn, ArgList const& argIn) {
      actionlist_.push_back( actIn );
      actioncmd_.push_back( argIn.ArgLine() );
      actionstatus_.push_back( INIT );
      return 0;
    }
    /// Set up actions for the given parm.
    int SetupActions(Topology **);
    /// Perform actions on the given frame.
    bool DoActions(Frame **, int);
    /// Call print for each action.
    void Print();
    /// List all actions in the action list.
    void List() const;
    bool Empty() const { return actionlist_.empty(); }
  private:
    /// Action initialization and setup status.
    enum ActionStatusType { NO_INIT=0, INIT, SETUP, INACTIVE };
    typedef std::vector<Action*> Aarray;
    /// List of actions
    Aarray actionlist_;
    /// List of action commands
    std::vector<std::string> actioncmd_;
    /// List of action statuses
    std::vector<ActionStatusType> actionstatus_;
    /// Default debug level for actions
    int debug_;
};
#endif
