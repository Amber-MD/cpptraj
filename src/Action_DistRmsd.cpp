#include "Action_DistRmsd.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_DistRmsd::Action_DistRmsd() : drmsd_(0) {}

void Action_DistRmsd::Help() const {
  mprintf("\t[<name>] [<mask>] [<refmask>] [out filename]\n%s"
          "  Calculate distance RMSD (DME) for specified atoms.\n", ReferenceAction::Help());
}

// Action_DistRmsd::Init()
/** Called once before traj processing. Set up reference info. */
Action::RetType Action_DistRmsd::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Check for keywords
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  // Reference keywords
  REF_.InitRef(actionArgs, init.DSL(), false, false);
  // Get the RMS mask string for target 
  std::string tMaskExpr = actionArgs.GetMaskNext();
  TgtMask_.SetMaskString( tMaskExpr );
  // Get the RMS mask string for reference
  std::string rMaskExpr = actionArgs.GetMaskNext();
  if (rMaskExpr.empty())
    rMaskExpr = tMaskExpr;
  REF_.SetRefMask( tMaskExpr );
 
  // Set up the RMSD data set
  drmsd_ = init.DSL().AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"DRMSD");
  if (drmsd_==0) return Action::ERR;
  // Add dataset to data file list
  if (outfile != 0) outfile->AddDataSet( drmsd_ );
# ifdef MPI
  if (REF_.SetTrajComm( init.TrajComm() )) return Action::ERR;
# endif
  mprintf("    DISTRMSD: (%s), reference is %s\n",TgtMask_.MaskString(),
          REF_.RefModeString().c_str());

  return Action::OK;
}

// Action_DistRmsd::Setup()
/** Called every time the trajectory changes. Set up TgtMask for the new 
  * parmtop and allocate space for selected atoms from the Frame.
  */
Action::RetType Action_DistRmsd::Setup(ActionSetup& setup) {

  if ( setup.Top().SetupIntegerMask(TgtMask_) ) return Action::ERR;
  if ( TgtMask_.None() ) {
    mprintf("Warning: No atoms in mask.\n");
    return Action::SKIP;
  }
  // Allocate space for selected atoms in the frame. This will also put the
  // correct masses in based on the mask.
  SelectedTgt_.SetupFrameFromMask(TgtMask_, setup.Top().Atoms());
  // Reference setup
  if (REF_.SetupRef(setup.Top(), TgtMask_.Nselected()))
    return Action::ERR; 

  return Action::OK;
}

// Action_DistRmsd::DoAction()
/** Called every time a frame is read in. Calc distance RMSD. */
Action::RetType Action_DistRmsd::DoAction(int frameNum, ActionFrame& frm) {
  // Perform any needed reference actions
  REF_.ActionRef( frm.TrajoutNum(), frm.Frm() );
  // Set selected frame atoms. Masses have already been set.
  SelectedTgt_.SetCoordinates(frm.Frm(), TgtMask_);
  double DR = SelectedTgt_.DISTRMSD( REF_.SelectedRef() );
  drmsd_->Add(frameNum, &DR);
  REF_.PreviousRef( frm.Frm() );
  return Action::OK;
}
