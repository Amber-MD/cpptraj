#include "Action_Align.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Align::Action_Align() :
  debug_(0),
  useMass_(false)
{}

void Action_Align::Help() const {
  mprintf("\t[<name>] <mask> [<refmask>] [mass]\n\t%s\n", ReferenceAction::Help());
  mprintf("  Align structure using specified <mask> onto reference.\n");
}

// Action_Align::Init()
Action::RetType Action_Align::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  // Check for keywords
  useMass_ = actionArgs.hasKey("mass");
  // Reference keywords: always fitting
  if (REF_.InitRef(actionArgs, init.DSL(), true, useMass_ )) return Action::ERR;
  // Get the fit mask string for target
  std::string tMaskExpr = actionArgs.GetMaskNext();
  tgtMask_.SetMaskString(tMaskExpr);
  // Get the fit mask string for reference
  std::string rMaskExpr = actionArgs.GetMaskNext();
  if (rMaskExpr.empty())
    rMaskExpr = tMaskExpr;
  REF_.SetRefMask( rMaskExpr );

# ifdef MPI
  if (REF_.SetTrajComm( init.TrajComm() )) return Action::ERR;
# endif
  mprintf("    ALIGN: Aligning atoms selected by mask '%s'\n", tgtMask_.MaskString());
  mprintf("\tReference is %s\n", REF_.RefModeString().c_str());
  if (useMass_)
    mprintf("\tFit will be mass-weighted.\n");

  return Action::OK;
}

// Action_Align::Setup()
/** Called every time the trajectory changes. Set up masks for target
  * and reference selection.
  */
Action::RetType Action_Align::Setup(ActionSetup& setup) {
  // Target setup
  if ( setup.Top().SetupIntegerMask( tgtMask_ ) ) return Action::ERR;
  mprintf("\tTarget mask:");
  tgtMask_.BriefMaskInfo();
  mprintf("\n");
  if ( tgtMask_.None() ) {
    mprintf("Warning: No atoms in mask '%s'.\n", tgtMask_.MaskString());
    return Action::SKIP;
  }
  // Allocate space for selected atoms in the frame. This will also put the
  // correct masses in based on the mask.
  tgtFrame_.SetupFrameFromMask(tgtMask_, setup.Top().Atoms());
  // Reference setup
  if (REF_.SetupRef(setup.Top(), tgtMask_.Nselected()))
    return Action::SKIP;
 
  // Warn if PBC and rotating
  if (setup.CoordInfo().TrajBox().Type() != Box::NOBOX) {
    mprintf("Warning: Coordinates are being rotated and box coordinates are present.\n"
            "Warning: Unit cell vectors are NOT rotated; imaging will not be possible\n"
            "Warning:  after the alignment is performed.\n");
  }
  return Action::OK;
}

// Action_Align::DoAction()
Action::RetType Action_Align::DoAction(int frameNum, ActionFrame& frm) {
  // Perform any needed reference actions
  REF_.ActionRef( frm.TrajoutNum(), frm.Frm() );
  // Set selected frame atoms. Masses have already been set.
  tgtFrame_.SetCoordinates(frm.Frm(), tgtMask_);
  tgtFrame_.RMSD_CenteredRef(REF_.SelectedRef(), rot_, tgtTrans_, useMass_);
  frm.ModifyFrm().Trans_Rot_Trans(tgtTrans_, rot_, REF_.RefTrans());
  REF_.PreviousRef( frm.Frm() );
  return Action::MODIFY_COORDS;
}
