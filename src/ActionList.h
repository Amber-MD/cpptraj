#ifndef INC_ACTIONLIST_H
#define INC_ACTIONLIST_H
#include "Action.h"
/// Hold Actions that will be performed every frame.
/** This class is responsible for holding all Actions that will be performed
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
    /// Set whether to supress Action Init/Setup output.
    void SetSilent(bool b) { actionsAreSilent_ = b; }
    /// Add given action to the action list and initialize.
    int AddAction(Action*, ArgList&, ActionInit&);
    /// Set up Actions for the given Topology.
    int SetupActions(ActionSetup&, bool);
    /// Perform Actions on the given Frame.
    bool DoActions(int, ActionFrame&);
    /// Call print for each Action.
    void PrintActions();
#   ifdef MPI
    /// Call parallel init for each Action
    void ParallelInitActions(Parallel::Comm const&);
    /// Call sync for each Action (parallel only)
    void SyncActions(Parallel::Comm const&);
#   endif
    /// List all Actions in the list.
    void List() const;
    /// \return Current debug level.
    int Debug()                      const { return debug_;                  }
    /// \return True if no Actions in list.
    bool Empty()                     const { return actionList_.empty();     }
    // The functions below help set up actions when ensemble processing.
    /// \return the number of Actions in the list.
    int Naction()                    const { return (int)actionList_.size(); }
    /// \return Arguments for corresponding Action.
    ArgList const& ActionArgs(int i) const { return actionList_[i].args_;    }
    /// \return Uninitialized copy of existing Action (for ensemble).
    Action* ActionAlloc(int i)       const { return (Action*)actionList_[i].ptr_->Alloc(); }
  private:
    /// Action initialization and setup status.
    enum ActionStatusType { NO_INIT=0, INIT, SETUP, INACTIVE };
    struct ActHolder {
      Action* ptr_;             ///< Pointer to Action.
      ArgList args_;            ///< Arguments associated with Action.
      ActionStatusType status_; ///< Current Action status.
    };
    typedef std::vector<ActHolder> Aarray;
    Aarray actionList_;     ///< List of Actions
    int debug_;             ///< Default debug level for new Actions
    bool actionsAreSilent_; ///< If true suppress all Init/Setup output from Actions.
};
#endif
