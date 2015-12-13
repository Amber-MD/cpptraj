#include "Action_Box.h"
#include "CpptrajStdio.h"

Action_Box::Action_Box() : nobox_(false) {}

void Action_Box::Help() const {
  mprintf("\t[x <xval>] [y <yval>] [z <zval>] [alpha <a>] [beta <b>] [gamma <g>]\n"
          "\t[nobox] [truncoct]\n"
          "  For each input frame, replace any box information with the information given.\n");
}

// Action_Box::Init()
Action::RetType Action_Box::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get keywords
  if ( actionArgs.hasKey("nobox") )
    nobox_ = true; 
  else {
    box_.SetX( actionArgs.getKeyDouble("x", 0.0) );
    box_.SetY( actionArgs.getKeyDouble("y", 0.0) );
    box_.SetZ( actionArgs.getKeyDouble("z", 0.0) );
    box_.SetAlpha( actionArgs.getKeyDouble("alpha", 0.0) );
    box_.SetBeta(  actionArgs.getKeyDouble("beta",  0.0) );
    box_.SetGamma( actionArgs.getKeyDouble("gamma", 0.0) );
    if (actionArgs.hasKey("truncoct")) box_.SetTruncOct();
  }

  mprintf("    BOX:");
  if (nobox_)
    mprintf(" Removing box information.\n");
  else {
    if (box_.BoxX() > 0) mprintf(" X=%.3f", box_.BoxX());
    if (box_.BoxY() > 0) mprintf(" Y=%.3f", box_.BoxY());
    if (box_.BoxZ() > 0) mprintf(" Z=%.3f", box_.BoxZ());
    if (box_.Alpha() > 0) mprintf(" A=%.3f", box_.Alpha());
    if (box_.Beta() > 0) mprintf(" B=%.3f", box_.Beta());
    if (box_.Gamma() > 0) mprintf(" G=%.3f", box_.Gamma());
    mprintf("\n");
  }
  return Action::OK;
}

// Action_Box::Setup()
Action::RetType Action_Box::Setup(ActionSetup& setup) {
  cInfo_ = setup.CoordInfo();
  if (nobox_) {
    mprintf("\tRemoving box info.\n");
    cInfo_.SetBox( Box() );
  } else {
    Box pbox( box_ );
    // Fill in missing box information from current box 
    pbox.SetMissingInfo( setup.CoordInfo().TrajBox() );
    mprintf("\tNew box type is %s\n", pbox.TypeName() );
    cInfo_.SetBox( pbox );
  }
  setup.SetCoordInfo( &cInfo_ );
  return Action::OK;
}

Action::RetType Action_Box::DoAction(int frameNum, ActionFrame& frm) {
  double* frame_box = frm.ModifyFrm().bAddress();
  if (nobox_) {
    frame_box[0] = 0.0;
    frame_box[1] = 0.0;
    frame_box[2] = 0.0;
    frame_box[3] = 0.0;
    frame_box[4] = 0.0;
    frame_box[5] = 0.0;
  } else {
    Box fbox( box_ );
    fbox.SetMissingInfo( Box( frame_box ) );
    frame_box[0] = fbox.BoxX();
    frame_box[1] = fbox.BoxY();
    frame_box[2] = fbox.BoxZ();
    frame_box[3] = fbox.Alpha();
    frame_box[4] = fbox.Beta();
    frame_box[5] = fbox.Gamma();
  }
  return Action::MODIFY_COORDS;
}
