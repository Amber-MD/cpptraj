#include "Action_MinMaxDist.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
Action_MinMaxDist::Action_MinMaxDist() :
  mode_(NO_MODE),
  distType_(NO_DIST)
{}

const char* Action_MinMaxDist::modeStr_[] = {
  "atoms",
  "residues",
  "molecules",
  0
};

const char* Action_MinMaxDist::distTypeStr_[] = {
  "minimum",
  "maximum",
  "minimum and maximum",
  0
};

// Action_MinMaxDist::Help()
void Action_MinMaxDist::Help() const {
  mprintf("mask1 <mask1> [mask2 <mask2>]\n"
          "  Record the min/max distance between atoms/residues/molecules.\n"
         );
}

// Action_MinMaxDist::Init()
Action::RetType Action_MinMaxDist::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  std::string mask1str = actionArgs.GetStringKey("mask1");
  if (mask1str.empty()) {
    mprinterr("Error: Must specify at least 'mask1'\n");
    return Action::ERR;
  }
  if (mask1_.SetMaskString( mask1str )) {
    mprinterr("Error: Could not set mask1 '%s'\n", mask1str.c_str());
    return Action::ERR;
  }
  std::string mask2str = actionArgs.GetStringKey("mask2");
  if (!mask2str.empty()) {
    if (mask2_.SetMaskString( mask2str )) {
      mprinterr("Error: Could not set mask2 '%s'\n", mask2str.c_str());
      return Action::ERR;
    }
  }
  // Default mode and distance
  if (mode_ == NO_MODE)
    mode_ = BY_ATOM;
  if (distType_ == NO_DIST) {
    if (actionArgs[0] == "mindist")
      distType_ = MIN_DIST;
    else if (actionArgs[0] == "maxdist")
      distType_ = MAX_DIST;
    else {
      mprintf("Warning: No distance type specified and command name '%s' unrecognized. Using default.\n");
      distType_ = MIN_DIST;
    }
  }

  mprintf("    MINMAXDIST: Calculating %s distance for selected %s\n",
          distTypeStr_[distType_], modeStr_[mode_]);
  mprintf("\tMask1: %s\n", mask1_.MaskString());
  if (mask2_.MaskStringSet()) {
    mprintf("\tMask2: %s\n", mask2_.MaskString());
  }

  return Action::OK;
}

// Action_MinMaxDist::Setup()
Action::RetType Action_MinMaxDist::Setup(ActionSetup& setup)
{
  
  return Action::OK;
}

// Action_MinMaxDist::DoAction()
Action::RetType Action_MinMaxDist::DoAction(int frameNum, ActionFrame& frm)
{
  return Action::OK;
}
