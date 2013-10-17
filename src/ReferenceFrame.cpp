#include "ReferenceFrame.h"
#include "Trajin_Single.h"
#include "CpptrajStdio.h"

// DESTRUCTOR
ReferenceFrame::~ReferenceFrame() {
  if (strippedParm_ && parm_ != 0) delete parm_;
  //if (frame_ != 0) delete frame_; // TODO: Enable.
}

int ReferenceFrame::LoadRef(std::string const& fname, Topology* parmIn, int debugIn)
{
  ArgList argIn;
  return LoadRef( fname, argIn, parmIn, debugIn);
}

// ReferenceFrame::LoadRef()
int ReferenceFrame::LoadRef(std::string const& fname, ArgList& argIn, 
                            Topology* parmIn, int debugIn)
{
  Trajin_Single traj;
  traj.SetDebug(debugIn);
  num_ = -1; // This is -1 as long as ref is not fully set up.
  // Set topology; remove old topology if it was stripped.
  if (strippedParm_ && parm_ != 0) delete parm_;
  strippedParm_ = false;
  parm_ = parmIn;
  // Set up trajectory - false = do not modify box info
  if ( traj.SetupTrajRead( fname, argIn, parm_, false ) ) {
    mprinterr("Error: reference: Could not set up trajectory.\n");
    return 1;
  }
  // Check for mask expression
  // NOTE: This is done after SetupTrajRead because forward slash is a valid
  //       mask operand for element, so getNextMask will unfortunately pick up 
  //       filenames as well as masks.
  std::string maskexpr = argIn.GetMaskNext();
  // Check for tag - done after SetupTrajRead so traj can process args
  std::string tag_ = argIn.getNextTag();
  // Tell trajectory to only process the start frame.
  traj.SingleFrame();
  // Check number of frames to be read
  int trajFrames = traj.TotalReadFrames();
  if (trajFrames < 1) {
    mprinterr("Error: No frames could be read for reference '%s'\n", traj.TrajFilename().full());
    return 1;
  }
  // Start trajectory read
  if ( traj.BeginTraj(false) ) {
    mprinterr("Error: Could not open reference '%s'\n.", traj.TrajFilename().full());
    return 1;
  }
  mprintf("\t");
  traj.PrintInfo(1);
  // Set up input frame
  if (frame_ == 0) frame_ = new Frame();
  frame_->SetupFrameV( parm_->Atoms(), traj.HasVelocity(), traj.NreplicaDimension() );
  // Read reference frame
  traj.GetNextFrame( *frame_ );
  if ( frame_->CheckCoordsInvalid() )
    mprintf("Warning: reference frame coords 1 & 2 overlap at origin; may be corrupt.\n");
  // If a mask expression was specified, strip to match the expression.
  if (!maskexpr.empty()) {
    if ( StripRef( maskexpr ) ) {
      delete frame_;
      frame_ = 0;
      return 1;
    }
  }
  // Set name and number.
  name_ = traj.TrajFilename();
  num_ = traj.Start();
  // If the top currently has no reference coords, set them now
  if (parm_->NoRefCoords()) parm_->SetReferenceCoords( frame_ );
  return 0;
}

/** Strip reference based on mask expression. */
int ReferenceFrame::StripRef(std::string const& maskexpr) {
  AtomMask stripMask( maskexpr );
  if (parm_->SetupIntegerMask(stripMask, *frame_)) return 1;
  return StripRef( stripMask );
}

/** Strip reference to match mask */
int ReferenceFrame::StripRef(AtomMask const& stripMask) {
  if (stripMask.None()) {
    mprinterr("Error: No atoms kept for reference.\n");
    return 1;
  }
  // Create new stripped frame
  Frame *strippedRefFrame = new Frame( *frame_, stripMask );
  mprintf("\tKept %i atoms in reference.\n", strippedRefFrame->Natom());
  // Create new stripped parm
  Topology *strippedRefParm = parm_->modifyStateByMask( stripMask );
  if (strippedRefParm==0) {
    mprinterr("Error: Could not create stripped reference topology.\n");
    return 1;
  }
  strippedParm_ = true;
  strippedRefParm->Summary();
  delete frame_;
  frame_ = strippedRefFrame;
  // Free parm_ if it was previously stripped.
  if (strippedParm_ && parm_ != 0) delete parm_;
  parm_ = strippedRefParm;
  return 0;
}

void ReferenceFrame::RefInfo() const {
  if (!tag_.empty())
    mprintf(" %s", tag_.c_str());
  mprintf(" '%s', frame %i\n", name_.full(), num_);
}
