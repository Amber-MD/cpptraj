// ActionList
#include "ActionList.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
ActionList::ActionList() : debug_(0), actionsAreSilent_(false) {}

// DESTRUCTOR
ActionList::~ActionList() { Clear(); }

void ActionList::Clear() {
  // No need to cast back to whatever action was allocd since Action destructor is virtual
  for (Aarray::const_iterator act = actionList_.begin(); act != actionList_.end(); ++act)
    delete act->ptr_;
  actionList_.clear();
}

// ActionList::SetDebug()
void ActionList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_ > 0)
    mprintf("ActionList DEBUG LEVEL SET TO %i\n",debug_);
}

// ActionList::AddAction()
int ActionList::AddAction(Action* actIn, ArgList& argIn, ActionInit& init)
{
  if (actIn == 0) {
    mprinterr("Internal Error: AddAction() called with null Action.\n");
    return 1;
  }
  int err = 0;
  if (actionsAreSilent_) SetWorldSilent( true );
  ActHolder act;
  act.ptr_ = actIn;
  act.args_ = argIn;
  // Attempt to initialize action
  if ( act.ptr_->Init( argIn, init, debug_ ) != Action::OK ) {
    mprinterr("Error: Could not initialize action [%s]\n", argIn.Command());
    delete act.ptr_;
    err = 1;
  } else {
    act.status_ = INIT;
    actionList_.push_back( act );
    if (argIn.CheckForMoreArgs()) err = 1;
  }
  if (actionsAreSilent_) SetWorldSilent( false );
  return err;
}

// ActionList::SetupActions()
/** Attempt to set up all actions in the action list with the given parm
  * If an action cannot be set up skip it.
  */
int ActionList::SetupActions(ActionSetup& setup, bool exitOnError) {
  if (actionList_.empty()) return 0;
  ActionSetup OriginalSetup = setup;
  mprintf(".....................................................\n");
  mprintf("ACTION SETUP FOR PARM '%s' (%zu actions):\n", setup.Top().c_str(), actionList_.size());
  if (actionsAreSilent_) SetWorldSilent( true );
  for (Aarray::iterator act = actionList_.begin(); act != actionList_.end(); ++act)
  {
    // Only attempt to set up action if active 
    if (act->status_ != INACTIVE) {
      mprintf("  %u: [%s]\n", act - actionList_.begin(), act->args_.ArgLine());
      act->status_ = SETUP;
      Action::RetType err = act->ptr_->Setup(setup);
      if (err == Action::ERR) {
        mprinterr("Error: Setup failed for [%s]\n", act->args_.Command());
        if (exitOnError) return 1;
        act->status_ = INIT;
      } else if (err == Action::SKIP) {
        mprintf("Warning: Setup incomplete for [%s]: Skipping\n", act->args_.Command());
        // Reset action status to INIT (pre-setup)
        act->status_ = INIT;
      } else if (err == Action::USE_ORIGINAL_FRAME) {
        setup = OriginalSetup;
      }
    }
  }
  if (actionsAreSilent_) SetWorldSilent( false );
  return 0;
}

// ActionList::DoActions()
/** Perform actions in the action list on the given Frame. Skip actions not 
  * initialized or not setup. 
  * \param frameNumIn The current frame number.
  * \param frm Contains current Frame.
  * \return true if coordinate output should be suppressed, false if coordinate
  *         output should be performed.
  */
bool ActionList::DoActions(int frameNumIn, ActionFrame& frm) {
  ActionFrame OriginalFrame = frm;
  //fprintf(stdout,"DEBUG: Performing %i actions on frame %i.\n",Naction,frameNumIn);
  for (Aarray::iterator act = actionList_.begin(); act != actionList_.end(); ++act) 
  {
    // Only do actions which were properly set up
    if (act->status_ == SETUP) { 
      // Perform action on frame
      Action::RetType err = act->ptr_->DoAction(frameNumIn, frm);
      // Check for action special conditions/errors
      if (err == Action::USE_ORIGINAL_FRAME) {
        // Return to original frame
        frm = OriginalFrame;
      } else if (err == Action::SUPPRESS_COORD_OUTPUT) {
        // Skip the rest of the actions and suppress output. Necessary when
        // e.g. performing a running average over coords.
        return true;
      } else if (err == Action::ERR) {
        // If here return type is ACTION_ERR.
        // Treat actions that fail as if they could not be set up
        mprintf("Warning: Action [%s] failed, frame %i.\n", act->args_.Command(), frameNumIn);
        act->status_ = INIT;
      }
    }
  }
  return false;
}

// ActionList::Print()
void ActionList::PrintActions() {
  for (Aarray::const_iterator act = actionList_.begin(); act != actionList_.end(); ++act)
  { // Skip deactivated actions
    if (act->status_ != INACTIVE)
      act->ptr_->Print();
  }
}
#ifdef MPI
int ActionList::ParallelInitActions(Parallel::Comm const& commIn) {
  int err = 0;
  for (Aarray::iterator act = actionList_.begin(); act != actionList_.end(); ++act)
  { // Skip deactivated actions
    if (act->status_ != INACTIVE) {
      //rprintf("DEBUG: Calling ParallelActionInit() for '%s'\n", act->args_.Command());
      if (act->ptr_->ParallelActionInit(commIn)) {
        rprintf("Warning: Parallel Init failed for Action '%s'\n", act->args_.Command());
        act->status_ = INACTIVE;
        err++;
      }
    }
  }
  return err;
}

int ActionList::NumPreviousFramesReqd() const {
  int nrequired = 0;
  for (Aarray::const_iterator act = actionList_.begin(); act != actionList_.end(); ++act)
  { // Skip deactivated actions
    if (act->status_ != INACTIVE)
      nrequired = std::max( nrequired, act->ptr_->ParallelPreviousFramesRequired() );
  }
  if (debug_ > 0) rprintf("DEBUG: Action(s) require %i previous frames.\n", nrequired);
  return nrequired;
}

/** Should not be called by master. */
int ActionList::ParallelProcessPreload(Action::FArray const& preload_frames) {
  int err = 0;
  for (Aarray::iterator act = actionList_.begin(); act != actionList_.end(); ++act)
  { // Skip deactivated actions
    if (act->status_ != INACTIVE) {
      //rprintf("DEBUG: Calling ParallelPreloadFrames() for '%s'\n", act->args_.Command());
      if (act->ptr_->ParallelPreloadFrames( preload_frames )) {
        rprintf("Warning: Parallel Preload failed for Action '%s'\n", act->args_.Command());
        act->status_ = INACTIVE;
        err++;
      }
    }
  }
  return err;
}

void ActionList::SyncActions(Parallel::Comm const& commIn) {
  for (Aarray::const_iterator act = actionList_.begin(); act != actionList_.end(); ++act)
  { // Skip deactivated actions
    if (act->status_ != INACTIVE) {
      //rprintf("DEBUG: Calling SyncAction() for '%s'\n", act->args_.Command());
      if (act->ptr_->SyncAction(commIn))
        rprintf("Warning: Sync failed for Action '%s'\n", act->args_.Command());
    }
  }
}
#endif
void ActionList::List() const {
  if (!actionList_.empty()) {
    mprintf("\nACTIONS (%zu total):\n", actionList_.size());
    for (Aarray::const_iterator act = actionList_.begin(); act != actionList_.end(); ++act)
      mprintf("  %u: [%s]\n", act - actionList_.begin(), act->args_.ArgLine());
  }
}
