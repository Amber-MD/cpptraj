#include "Action_ConvertToFrac.h"
#include "CpptrajStdio.h"

// Action_ConvertToFrac::Help()
void Action_ConvertToFrac::Help() const {
  mprintf("THIS ACTION IS CURRENTLY ONLY FOR DEBUGGING PURPOSES.\n"
          "    Converts all Cartesian coordinates to fractional space.\n");
}

// Action_ConvertToFrac::Init()
Action::RetType Action_ConvertToFrac::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  mprintf("    CONVERT TO FRAC\n");
  mprintf("\tConverting all coordinates to fractional space.\n");
  return Action::OK;
}

// Action_ConvertToFrac::Setup()
Action::RetType Action_ConvertToFrac::Setup(ActionSetup& setup)
{
  // Determine Box info
  if (!setup.CoordInfo().TrajBox().HasBox()) {
    mprintf("Warning: Topology %s does not contain box information.\n", setup.Top().c_str());
    return Action::SKIP;
  }
  // Allocate space in new frame
  newFrame_.SetupFrameV(setup.Top().Atoms(), setup.CoordInfo());

  return Action::OK;
}

// Action_ConvertToFrac::DoAction()
Action::RetType Action_ConvertToFrac::DoAction(int frameNum, ActionFrame& frm)
{
  Frame const& frameIn = frm.Frm();
  Box const& box = frameIn.BoxCrd();

  for (int at = 0; at != frameIn.Natom(); at++) {
    Vec3 frac = box.FracCell() * Vec3(frameIn.XYZ(at));
    newFrame_.SetXYZ( at, frac );
  }
  // Set frame
  frm.SetFrame( &newFrame_ );
  // TODO Modify unit cell info?
  return Action::MODIFY_COORDS;
}
