#include "Action_Time.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
Action_Time::Action_Time() :
  time0_(0.0),
  dt_(0.001),
  mode_(ADD)
{}

// Action_Time::Help()
void Action_Time::Help() const {
  mprintf("\t{time0 <initial time> dt <step> [update] | notime}\n"
          "  Add, modify, or remove time information from each frame.\n"
          "    time0 <initial time> : Time of the first frame (ps).\n"
          "    dt <step>            : Time step between frames (ps).\n"
          "    update               : If specified, modify any existing time info.\n"
          "    notime               : Remove any time info from frame.\n");
}

// Action_Time::Init()
Action::RetType Action_Time::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  /// Make sure that something has been specified.
  if (actionArgs.hasKey("notime"))
    mode_ = REMOVE;
  else {
    if (!actionArgs.Contains("time0") &&
        !actionArgs.Contains("dt"))
    {
      mprinterr("Error: Both 'time0' and 'dt' must be specified.\n");
      return Action::ERR;
    }
    time0_ = actionArgs.getKeyDouble("time0", 0.0);
    dt_    = actionArgs.getKeyDouble("dt", 0.001);
    if (actionArgs.hasKey("update"))
      mode_ = MODIFY;
    else
      mode_ = ADD;
  }

  mprintf("    TIME:");
  if (mode_ == REMOVE)
    mprintf(" Removing all time information from frames.\n");
  else {
    if (mode_ == MODIFY)
      mprintf(" Updating time information in frames.\n");
    else
      mprintf(" Adding/overwriting time information in frames.\n");
    mprintf("\tInitial time = %g ps\n", time0_);
    mprintf("\tTime step    = %g ps\n", dt_);
  }
  return Action::OK;
}

// Action_Time::Setup()
Action::RetType Action_Time::Setup(ActionSetup& setup)
{
  cInfo_ = setup.CoordInfo();
  if (!cInfo_.HasTime()) {
    if (mode_ == MODIFY)
      mprintf("Warning: 'update' specified but no time info in frame. Adding time info.\n");
    else if (mode_ == REMOVE) {
      mprintf("Warning: 'notime' specified but no time info in frame. Skipping.\n");
      return Action::SKIP;
    } else if (mode_ == ADD)
      mprintf("\tAdding time information to frames.\n");
  } else {
    if (mode_ == MODIFY)
      mprintf("\tUpdating time information in frames.\n");
    else if (mode_ == REMOVE)
      mprintf("\tRemoving time information in frames.\n");
    else if (mode_ == ADD)
      mprintf("\tOverwriting time information in frames.\n");
  }
  if (mode_ == REMOVE)
    cInfo_.SetTime( false );
  else
    cInfo_.SetTime( true );
  setup.SetCoordInfo( &cInfo_ );
  return Action::OK;
}

// Action_Time::DoAction()
Action::RetType Action_Time::DoAction(int frameNum, ActionFrame& frm)
{
  double currTime = 0.0;
  double newTime  = 0.0;
  switch (mode_) {
    case REMOVE : frm.ModifyFrm().SetTime( currTime ); break;
    case MODIFY : currTime = frm.Frm().Time();
    case ADD    :
      newTime = currTime + time0_ + ((double)frameNum * dt_);
      frm.ModifyFrm().SetTime( newTime );
      break;
  }
  return Action::MODIFY_COORDS;
}
