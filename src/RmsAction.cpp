#include "RmsAction.h"
#include "CpptrajStdio.h"

RmsAction::RmsAction() :
  fit_(true),
  rotate_(true),
  useMass_(false)
{}

void RmsAction::GetRmsKeywords(ArgList &actionArgs) {
  fit_ = !actionArgs.hasKey("nofit");
  if (fit_)
    rotate_ = !actionArgs.hasKey("norotate");
  useMass_ = actionArgs.hasKey("mass");
}

/** \return Reference mask. */
std::string RmsAction::GetRmsMasks(ArgList& actionArgs) {
  std::string mask0 = actionArgs.GetMaskNext();
  tgtMask_.SetMaskString(mask0);
  // Get the RMS mask string for reference
  std::string mask1 = actionArgs.GetMaskNext();
  if (mask1.empty())
    mask1 = mask0;
  return mask1;
}

void RmsAction::PrintRmsStatus() {
  if (!fit_)
    mprintf(", no fitting");
  else {
    mprintf(", with fitting");
    if (!rotate_)
      mprintf(" (no rotation)");
  }
  if (useMass_)
    mprintf(", mass-weighted");
  mprintf(".\n");
}

int RmsAction::SetupRmsMask(Topology const& topIn, const char* call) {
  if ( topIn.SetupIntegerMask( tgtMask_ ) ) return 1;
  tgtMask_.MaskInfo();
  if ( tgtMask_.None() ) {
    mprintf("Warning: %s: No atoms in mask.\n",call);
    return 1;
  }
  // Allocate space for selected atoms in the frame. This will also put the
  // correct masses in based on the mask.
  tgtFrame_.SetupFrameFromMask(tgtMask_, topIn.Atoms());
  return 0;
}
