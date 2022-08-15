#include "GridMover.h"
#include "Topology.h"

using namespace Cpptraj;

/** CONSTRUCTOR. */
GridMover::GridMover() :
  gridMoveType_(NO_MOVE),
  firstFrame_(true),
  x_align_(false)
{}

/** Set GridMover options. */
int GridMover::MoverInit(MoveType typeIn, bool x_alignIn) {
  gridMoveType_ = typeIn;
  x_align_ = x_alignIn;
  firstFrame_ = true;
  return 0;
}

/** Setup GridMover */
int GridMover::MoverSetup(Topology const& topIn, AtomMask const& maskIn) {
  // Set up frames if needed
  if (gridMoveType_ == RMS_FIT) {
    tgt_.SetupFrameFromMask(maskIn, topIn.Atoms());
    ref_ = tgt_;
    //firstFrame_ = true;
  }
  return 0;
}
